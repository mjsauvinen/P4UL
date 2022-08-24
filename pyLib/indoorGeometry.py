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
  def __init__(self, x1, x2, y1, y2, z1, z2, w ):
    self.xb = [ x1, x2 ]
    self.yb = [ y1, y2 ]
    self.zb = [ z1, z2 ]
    self.w  = w          # Wall thickness

#= = = = = = = = = = = = = = = = = = = = = # 
class Hole_with_plate:
  def __init__(self, xc, yc, zc, L , Loffset, normal):
    if( normal not in ['x','y','z'] ):
      sys.exit(' normal direction {} not correctly specified. Exiting ...'.format(normal)) 
    self.xc = xc
    self.yc = yc
    self.zc = zc
    self.L  = L
    self.Loff   = Loffset
    self.normal = normal # 'x', 'y', 'z' 


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
    
    self.Ntot = (self.Nx * self.Ny * self.Nz)
    # The 3d raster needs to be arranged as S[j,i,k] to maintain compatibility with 2d rasters.
    self.S  = np.zeros( (self.Ny, self.Nx, self.Nz), np.int8 )
    ostr = '''
    {}:
    nx={}, ny={}, nz={}
    dx={:2.5f}, dy={:2.5f}, dz={:2.5f}
    '''.format(name,self.Nx, self.Ny, self.Nz, self.dx, self.dy, self.dz)
    print(ostr)
    
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

  def makeHoleWithPlate( self, Hx, nb=0 ):
    dx = self.dx; dy = self.dy; dz = self.dz
    k1 = ibound( Hx.zc-Hx.L/2 , dz , self.Nz )
    k2 = ibound( Hx.zc+Hx.L/2 , dz , self.Nz )
    j1 = ibound( Hx.yc-Hx.L/2 , dy , self.Ny )
    j2 = ibound( Hx.yc+Hx.L/2 , dy , self.Ny )
    i1 = ibound( Hx.xc-Hx.L/2 , dx , self.Nx )
    i2 = ibound( Hx.xc+Hx.L/2 , dx , self.Nx )
    
    if( Hx.normal == 'z'):
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
      
    if( Hx.normal == 'y'):
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
    
    if( Hx.normal == 'x'):
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
    
    if( Hx.normal == 'z'):
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
    
    if( Hx.normal == 'y'):
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
    
    if( Hx.normal == 'x'):
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

#= = = = = = = = = = = = = = = = = = = = = #
