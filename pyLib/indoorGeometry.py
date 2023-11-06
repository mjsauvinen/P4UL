import numpy as np
import scipy.ndimage as sn

def ibound(x, dx, N):
  i = np.round(x/dx).astype(int)
  i = np.minimum( i , N )
  i = np.maximum( i , 0 )
  #print('i = {}'.format(i))
  return i

#= = = = = = = = = = = = = = = = = = = = = # 
class Box:
  def __init__(self, x1, x2, y1, y2, z1, z2, w=0. ):
    self.xb = [ x1, x2 ]
    self.yb = [ y1, y2 ]
    self.zb = [ z1, z2 ]
    self.w  = w          # Wall thickness

#= = = = = = = = = = = = = = = = = = = = = # 
class Hole_with_plate:
  def __init__(self, xc, yc, zc, L , Loffset, normal):
    if( normal not in ['x','y','z','+x','+y','+z','-x','-y','-z'] ):
      sys.exit(' normal direction {} not correctly specified. Exiting ...'.format(normal)) 
    self.xc = xc
    self.yc = yc
    self.zc = zc
    self.L  = L
    self.normal = normal # 'x', 'y', 'z' 
    self.Loff = Loffset
    if('-' in normal): self.Loff = -np.abs(Loffset)
    
#= = = = = = = = = = = = = = = = = = = = = #

class Cylinder:
  def __init__(self, xc, yc, zc, radius, H, normal ):
    if( normal not in ['x','y','z','+x','+y','+z','-x','-y','-z'] ):
      sys.exit(' normal direction {} not correctly specified. Exiting ...'.format(normal))
    self.xc = xc
    self.yc = yc
    self.zc = zc
    self.radius = radius
    self.H = H         # height
    self.normal = normal
    

#= = = = = = = = = = = = = = = = = = = = = #
class Domain:
  def __init__(self, name, Lx, Ly, Lz, dx, dy, dz):
    self.name = str(name)
    self.Lx = Lx    # length (m)
    self.Ly = Ly    # width  (m)
    self.Lz = Lz    # height (m)
    self.dx = dx    # resolution
    self.dy = dy    # resolution
    self.dz = dz    # resolution
    self.Nx = np.round( Lx/dx ).astype(int) # Number of grid points
    self.Ny = np.round( Ly/dy ).astype(int)
    self.Nz = np.round( Lz/dz ).astype(int)
    
    if( (Lx/dx)/self.Nx != 1.0 ): 
      print(' {}: Warning! Rounding Lx/dx = {}'.format(self.name,Lx/dx))
    if( (Ly/dy)/self.Ny != 1.0 ):
      print(' {}: Warning! Rounding Ly/dy = {}'.format(self.name,Ly/dy))    
    if( (Lz/dz)/self.Nz != 1.0 ):
      print(' {}: Warning! Rounding Lz/dz = {}'.format(self.name,Lz/dz))
    
    self.Ntot = (self.Nx * self.Ny * self.Nz)
    # The 3d raster needs to be arranged as S[j,i,k] to maintain compatibility with 2d rasters.
    self.S  = np.zeros( (self.Ny, self.Nx, self.Nz), np.int8 )
    
    ostr = '''
    {}:
    nx={}, ny={}, nz={}
    dx={:2.6f}, dy={:2.6f}, dz={:2.6f}
    '''.format(name,self.Nx, self.Ny, self.Nz, self.dx, self.dy, self.dz)
    print(ostr)
  
  #= = = = = = = = = = = = = = = = = = = = = #
  
  def addCylindricalRing( self, Cx ):
    # X Normal
    if( 'x' in Cx.normal ):
      y = np.linspace(0., (self.Ny-1)*self.dy, self.Ny)
      z = np.linspace(0., (self.Nz-1)*self.dz, self.Nz)
      Y,Z = np.meshgrid(y,z, indexing='ij') 
      idc  =  ( np.sqrt( (Y-Cx.yc)**2 + (Z-Cx.zc)**2 ) < (Cx.radius+self.dy) )
      idc *= ~( np.sqrt( (Y-Cx.yc)**2 + (Z-Cx.zc)**2 ) < (Cx.radius-self.dy) )
      
      if( '-' in Cx.normal ):
        ic1  = ibound( Cx.xc-Cx.H, self.dx , self.Nx )
        ic2  = ibound( Cx.xc     , self.dx , self.Nx )
      else:
        ic1  = ibound( Cx.xc     , self.dx , self.Nx )
        ic2  = ibound( Cx.xc+Cx.H, self.dx , self.Nx )
      
      if( ic2 == ic1 ): ic2 += 1
      
      for ic in range(ic1,ic2):
        Stmp = self.S[:,ic,:]
        Stmp += idc.astype(int)
        Stmp[(Stmp > 1)] = 1
        self.S[:,ic,:] = Stmp
    
    # Y Normal     
    if( 'y' in Cx.normal ):
      x = np.linspace(0., (self.Nx-1)*self.dx, self.Nx)
      z = np.linspace(0., (self.Nz-1)*self.dz, self.Nz)
      X,Z = np.meshgrid(x,z, indexing='ij') 
      idc  =  ( np.sqrt( (X-Cx.xc)**2 + (Z-Cx.zc)**2 ) < (Cx.radius+self.dx) )
      idc *= ~( np.sqrt( (X-Cx.xc)**2 + (Z-Cx.zc)**2 ) < (Cx.radius-self.dx) )
      
      if( '-' in Cx.normal ):
        jc1  = ibound( Cx.yc-Cx.H, self.dy , self.Ny )
        jc2  = ibound( Cx.yc     , self.dy , self.Ny )
      else:
        jc1  = ibound( Cx.yc     , self.dy , self.Ny )
        jc2  = ibound( Cx.yc+Cx.H, self.dy , self.Ny )
      
      if( jc2 == jc1 ): jc2 += 1
      
      for jc in range(jc1,jc2):
        Stmp = self.S[jc,:,:]
        Stmp += idc.astype(int)
        Stmp[(Stmp > 1)] = 1
        self.S[jc,:,:] = Stmp
        
    # Z Normal 
    if( 'z' in Cx.normal ):
      x = np.linspace(0., (self.Nx-1)*self.dx, self.Nx)
      y = np.linspace(0., (self.Ny-1)*self.dy, self.Ny)
      Y,X = np.meshgrid(y,x, indexing='ij') 
      idc  =  ( np.sqrt( (X-Cx.xc)**2 + (Y-Cx.yc)**2 ) < (Cx.radius+self.dx) )
      idc *= ~( np.sqrt( (X-Cx.xc)**2 + (Z-Cx.yc)**2 ) < (Cx.radius-self.dx) )
      
      if( '-' in Cx.normal ):
        kc1  = ibound( Cx.zc-Cx.H, self.dz , self.Nz )
        kc2  = ibound( Cx.zc     , self.dz , self.Nz )
      else:
        kc1  = ibound( Cx.zc     , self.dz , self.Nz )
        kc2  = ibound( Cx.zc+Cx.H, self.dz , self.Nz )
      
      if( kc2 == kc1 ): kc2 += 1
      
      for kc in range(kc1,kc2):
        Stmp = self.S[:,:,kc]
        Stmp += idc.astype(int)
        Stmp[(Stmp > 1)] = 1
        self.S[:,:,kc] = Stmp 
  
  #= = = = = = = = = = = = = = = = = = = = = # 
  def addCylinder( self, Cx ):
    # X Normal
    if( 'x' in Cx.normal ):
      y = np.linspace(0., (self.Ny-1)*self.dy, self.Ny)
      z = np.linspace(0., (self.Nz-1)*self.dz, self.Nz)
      Y,Z = np.meshgrid(y,z, indexing='ij') 
      idc  =  ( np.sqrt( (Y-Cx.yc)**2 + (Z-Cx.zc)**2 ) < (Cx.radius) )
      
      if( '-' in Cx.normal ):
        ic1  = ibound( Cx.xc-Cx.H, self.dx , self.Nx )
        ic2  = ibound( Cx.xc     , self.dx , self.Nx )
      else:
        ic1  = ibound( Cx.xc     , self.dx , self.Nx )
        ic2  = ibound( Cx.xc+Cx.H, self.dx , self.Nx )
      
      if( ic2 == ic1 ): ic2 += 1
      
      for ic in range(ic1,ic2):
        Stmp = self.S[:,ic,:]
        Stmp += idc.astype(int)
        Stmp[(Stmp > 1)] = 1
        self.S[:,ic,:] = Stmp
    
    # Y Normal     
    if( 'y' in Cx.normal ):
      x = np.linspace(0., (self.Nx-1)*self.dx, self.Nx)
      z = np.linspace(0., (self.Nz-1)*self.dz, self.Nz)
      X,Z = np.meshgrid(x,z, indexing='ij') 
      idc  =  ( np.sqrt( (X-Cx.xc)**2 + (Z-Cx.zc)**2 ) < (Cx.radius) )
      
      if( '-' in Cx.normal ):
        jc1  = ibound( Cx.yc-Cx.H, self.dy , self.Ny )
        jc2  = ibound( Cx.yc     , self.dy , self.Ny )
      else:
        jc1  = ibound( Cx.yc     , self.dy , self.Ny )
        jc2  = ibound( Cx.yc+Cx.H, self.dy , self.Ny )
      
      if( jc2 == jc1 ): jc2 += 1
      
      for jc in range(jc1,jc2):
        Stmp = self.S[jc,:,:]
        Stmp += idc.astype(int)
        Stmp[(Stmp > 1)] = 1
        self.S[jc,:,:] = Stmp
        
    # Z Normal 
    if( 'z' in Cx.normal ):
      x = np.linspace(0., (self.Nx-1)*self.dx, self.Nx)
      y = np.linspace(0., (self.Ny-1)*self.dy, self.Ny)
      Y,X = np.meshgrid(y,x, indexing='ij') 
      idc  =  ( np.sqrt( (X-Cx.xc)**2 + (Y-Cx.yc)**2 ) < (Cx.radius) )
      
      if( '-' in Cx.normal ):
        kc1  = ibound( Cx.zc-Cx.H, self.dz , self.Nz )
        kc2  = ibound( Cx.zc     , self.dz , self.Nz )
      else:
        kc1  = ibound( Cx.zc     , self.dz , self.Nz )
        kc2  = ibound( Cx.zc+Cx.H, self.dz , self.Nz )
      
      if( kc2 == kc1 ): kc2 += 1
      
      for kc in range(kc1,kc2):
        Stmp = self.S[:,:,kc]
        Stmp += idc.astype(int)
        Stmp[(Stmp > 1)] = 1
        self.S[:,:,kc] = Stmp 
  
  #= = = = = = = = = = = = = = = = = = = = = # 
  def addBox( self, Bx ):
    dx = self.dx; dy = self.dy; dz = self.dz
    k1 = ibound( Bx.zb[0] , dz , self.Nz )
    k2 = ibound( Bx.zb[1] , dz , self.Nz )
    j1 = ibound( Bx.yb[0] , dy , self.Ny )
    j2 = ibound( Bx.yb[1] , dy , self.Ny )
    i1 = ibound( Bx.xb[0] , dx , self.Nx )
    i2 = ibound( Bx.xb[1] , dx , self.Nx )

    iw1= ibound( Bx.xb[0]-Bx.w[0]/2., dx, self.Nx )
    iw2= ibound( Bx.xb[0]+Bx.w[0]/2., dx, self.Nx )
    self.S[j1:j2,iw1:iw2,k1:k2] = 1  # left wall
    
    iw1= ibound( Bx.xb[1]-Bx.w[0]/2., dx, self.Nx )
    iw2= ibound( Bx.xb[1]+Bx.w[0]/2., dx, self.Nx )
    self.S[j1:j2,iw1:iw2,k1:k2] = 1  # right wall (back wall)  
    
    jw1= ibound( Bx.yb[0]-Bx.w[1]/2., dy, self.Ny )
    jw2= ibound( Bx.yb[0]+Bx.w[1]/2., dy, self.Ny )
    self.S[jw1:jw2,i1:i2,k1:k2] = 1  # side wall
    
    jw1= ibound( Bx.yb[1]-Bx.w[1]/2., dy, self.Ny )
    jw2= ibound( Bx.yb[1]+Bx.w[1]/2., dy, self.Ny )
    self.S[jw1:jw2,i1:i2,k1:k2] = 1  # side wall  
    
    kw1= ibound( Bx.zb[1]-Bx.w[2]/2., dz, self.Nz )
    kw2= ibound( Bx.zb[1]+Bx.w[2]/2., dz, self.Nz )
    self.S[j1:j2,i1:i2,kw1:kw2] = 1 # top wall     
    
    #print(' k1={}, k2={}, j1={}, j2={}, i1={}, i2={}'.format(k1,k2,j1,j2,i1,i2))
    #print(' Number of non-zeros {}'.format( np.count_nonzero( self.S[kw1:kw2,j1:j2,i1:i2] ) ))
  
  #= = = = = = = = = = = = = = = = = = = = = # 
  def addSolidBlock( self, Bx ):
    dx = self.dx; dy = self.dy; dz = self.dz
    k1 = ibound( Bx.zb[0] , dz , self.Nz )
    k2 = ibound( Bx.zb[1] , dz , self.Nz )
    j1 = ibound( Bx.yb[0] , dy , self.Ny )
    j2 = ibound( Bx.yb[1] , dy , self.Ny )
    i1 = ibound( Bx.xb[0] , dx , self.Nx )
    i2 = ibound( Bx.xb[1] , dx , self.Nx )
    
    self.S[j1:j2,i1:i2,k1:k2] = 1 # fill in the whole box 

  #= = = = = = = = = = = = = = = = = = = = = #
  def addEmptyBlock( self, Bx ):
    dx = self.dx; dy = self.dy; dz = self.dz
    k1 = ibound( Bx.zb[0] , dz , self.Nz )
    k2 = ibound( Bx.zb[1] , dz , self.Nz )
    j1 = ibound( Bx.yb[0] , dy , self.Ny )
    j2 = ibound( Bx.yb[1] , dy , self.Ny )
    i1 = ibound( Bx.xb[0] , dx , self.Nx )
    i2 = ibound( Bx.xb[1] , dx , self.Nx )
    
    self.S[j1:j2,i1:i2,k1:k2] = 0 # fill in the whole box  

  #= = = = = = = = = = = = = = = = = = = = = # 
  
  def makeHoleWithPlate( self, Hx, nb=0 ):
    dx = self.dx; dy = self.dy; dz = self.dz
    k1 = ibound( Hx.zc-Hx.L/2 , dz , self.Nz )
    k2 = ibound( Hx.zc+Hx.L/2 , dz , self.Nz )
    j1 = ibound( Hx.yc-Hx.L/2 , dy , self.Ny )
    j2 = ibound( Hx.yc+Hx.L/2 , dy , self.Ny )
    i1 = ibound( Hx.xc-Hx.L/2 , dx , self.Nx )
    i2 = ibound( Hx.xc+Hx.L/2 , dx , self.Nx )
    
    if('z' in Hx.normal):
      koff = np.round( Hx.Loff/dz ).astype(int)
      ko1 = k1 + koff  #; print(' k1 = {}, koff = {}'.format(k1, koff))
      ko2 = k2 + koff  #; print(' k2 = {}, koff = {}'.format(k2, koff))
      self.S[j1:j2, i1:i2,  ko1:ko2] = self.S[j1:j2, i1:i2, k1:k2].copy()
      self.S[j1:j2, i1:i2,  k1:k2] = 0
      
      for kl in np.arange(ko1,ko2):
        if( np.count_nonzero(self.S[j1:j2, i1:i2, kl]) > 0 and (nb > 0) ):
          self.S[j1-nb:j2+nb, i1-nb:i2+nb, kl] \
          = sn.binary_dilation( self.S[j1-nb:j2+nb, i1-nb:i2+nb, kl] , iterations=nb )
      
      ostr = '''
      {}:
      !- z-normal
      i_sbl = {},
      i_sbr = {},
      j_sbs = {},
      j_sbn = {},
      k_sbb = {},
      k_sbt = {},
      '''.format(self.name,i1,i2,self.Ny-j2,self.Ny-j1,k1,k2)
      print(ostr)
      
    if('y' in Hx.normal):
      joff = np.round( Hx.Loff/dy ).astype(int)
      jo1 = j1 + joff
      jo2 = j2 + joff 
      self.S[jo1:jo2, i1:i2, k1:k2] = self.S[j1:j2, i1:i2,  k1:k2].copy()
      self.S[ j1:j2 , i1:i2, k1:k2] = 0
      
      for jl in np.arange(jo1,jo2):
        if( np.count_nonzero(self.S[jl, i1:i2, k1:k2]) > 0 and (nb > 0) ):
          self.S[jl, i1-nb:i2+nb, k1-nb:k2+nb] \
          = sn.binary_dilation( self.S[jl, i1-nb:i2+nb, k1-nb:k2+nb] , iterations=nb )
      
      ostr = '''
      {}:
      !- y-normal
      i_sbl = {},
      i_sbr = {},
      j_sbs = {},
      j_sbn = {},
      k_sbb = {},
      k_sbt = {},
      '''.format(self.name,i1,i2,self.Ny-j2,self.Ny-j1,k1,k2)
      print(ostr)
    
    if('x' in Hx.normal):
      ioff = np.round( Hx.Loff/dx ).astype(int)
      io1 = i1 + ioff
      io2 = i2 + ioff 
      self.S[j1:j2, io1:io2, k1:k2] = self.S[j1:j2, i1:i2, k1:k2].copy()
      self.S[j1:j2, i1:i2, k1:k2] = 0
      
      for il in np.arange(io1,io2):
        if( np.count_nonzero(self.S[j1:j2, il, k1:k2]) > 0 and (nb > 0) ):
          self.S[j1-nb:j2+nb, il, k1-nb:k2+nb] \
          = sn.binary_dilation( self.S[j1-nb:j2+nb, il, k1-nb:k2+nb] , iterations=nb )
      
      ostr = '''
      {}:
      !- x-normal
      i_sbl = {},
      i_sbr = {},
      j_sbs = {},
      j_sbn = {},
      k_sbb = {},
      k_sbt = {},
      '''.format(self.name,i1,i2,self.Ny-j2,self.Ny-j1,k1,k2)
      print(ostr)
      
  #= = = = = = = = = = = = = = = = = = = = = #
  
  def makeHoleOnly( self, Hx ):
    dx = self.dx; dy = self.dy; dz = self.dz
    k1 = ibound( Hx.zc-Hx.L/2 , dz , self.Nz )
    k2 = ibound( Hx.zc+Hx.L/2 , dz , self.Nz )
    j1 = ibound( Hx.yc-Hx.L/2 , dy , self.Ny )
    j2 = ibound( Hx.yc+Hx.L/2 , dy , self.Ny )
    i1 = ibound( Hx.xc-Hx.L/2 , dx , self.Nx )
    i2 = ibound( Hx.xc+Hx.L/2 , dx , self.Nx )
    
    if('z' in Hx.normal):
      self.S[j1:j2, i1:i2,  k1:k2] = 0
      ostr = '''
      {}:
      !- z-normal
      i_sbl = {},
      i_sbr = {},
      j_sbs = {},
      j_sbn = {},
      k_sbb = {},
      k_sbt = {},
      '''.format(self.name,i1,i2,self.Ny-j2,self.Ny-j1,k1,k2)
      print(ostr)
    
    if('y' in Hx.normal):
      self.S[ j1:j2 , i1:i2, k1:k2] = 0
      ostr = '''
      {}:
      !- y-normal
      i_sbl = {},
      i_sbr = {},
      j_sbs = {},
      j_sbn = {},
      k_sbb = {},
      k_sbt = {},
      '''.format(self.name,i1,i2,self.Ny-j2,self.Ny-j1,k1,k2)
      print(ostr)
    
    if('x' in Hx.normal):
      self.S[j1:j2, i1:i2, k1:k2] = 0
      ostr = '''
      {}:
      !- x-normal
      i_sbl = {},
      i_sbr = {},
      j_sbs = {},
      j_sbn = {},
      k_sbb = {},
      k_sbt = {},
      '''.format(self.name,i1,i2,self.Ny-j2,self.Ny-j1,k1,k2)
      print(ostr)
  
  #= = = = = = = = = = = = = = = = = = = = = #
  
  def save( self ):
    
    filename = self.name+'.npy'
    with open( filename, 'wb') as fout:
      pickle.dump(self, fout, pickle.HIGHEST_PROTOCOL)
  
  #= = = = = = = = = = = = = = = = = = = = = #
  
  def load( self ):
  
    filename = self.name+'.npy'
    with open( filename, 'rb') as fin:
      Dobj = pickle.load(fin)
      
    return Dobj
  
  #= = = = = = = = = = = = = = = = = = = = = #
  
  def writeout( self, fileout ):
    Sdict = dict()
    Sdict['S'] = self.S  # [j,i,k] format
    Sdict['dPx'] = np.array([ self.dy, self.dx, self.dz ])
    Sdict['GlobOrig'] = np.zeros( 3 ) # N, E, Z
    fileout = fileout.strip('.npz')+'.npz'
    np.savez_compressed(fileout, **Sdict )
    
  #= = = = = = = = = = = = = = = = = = = = = #
  
  def writeWallTemp( self, Temp, fileout ):
    T = self.S.copy().astype(np.float16)
    ids = (T > 0.1) 
    T[ids] = Temp  # 273. + 20.
    idnan = (T<(Temp-2.))
    T[idnan] = np.nan
    
    Tdict = dict()
    Tdict['S'] = T
    Tdict['dPx'] = np.array([ self.dy, self.dx, self.dz ])
    Tdict['GlobOrig'] = np.zeros( 3 )
    fileout = fileout.strip('.npz')+'_temp.npz'
    np.savez_compressed(fileout, **Tdict )
    
  #= = = = = = = = = = = = = = = = = = = = = #
  
  def coordsAsPalmIJK( self, strId, xc, yc, zc ):
    i = ibound( xc , self.dx , self.Nx )
    j = ibound( yc , self.dy , self.Ny )
    k = ibound( zc , self.dz , self.Nz )
    j = (self.Ny-j)  # Palm origin is bottom left
    
    ostr = ''' ID: {}
    {}: x={:4.1f}, y={:4.1f}, z={:4.1f}
      --> i={:4d}, j={:4d}, k={:4d}
    '''.format(strId,self.name, xc,yc,zc,i,j,k)
    print(ostr)
     
  #= = = = = = = = = = = = = = = = = = = = = #
