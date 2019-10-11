import operator
import numpy as np
import sys
''' 
Description:


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def writeNumpyZFootprintRaw( filename, arr ):
  fstr = filename.strip('.npz')
  print(' Writing raw footprint data to file {}.npz ...'.format(fstr))
  dims = np.shape(arr)
  if( dims[1] != 9 ):
    sys.exit(" Error: dims[1] does not equal 9. Exiting ...")
    
  np.savez_compressed(fstr, \
    xO=arr[:,0], yO=arr[:,1], zO=arr[:,2], \
      xt=arr[:,3], yt=arr[:,4], zt=arr[:,5], \
        ut=arr[:,6], vt=arr[:,7], wt=arr[:,8] )
  
  print(' ... done! ')
  
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def writeNumpyZFootprintIJK(fn, xO, yO, zO, xt, yt, zt, ut, vt, wt, dxyz):
  fstr = fn.split('.npz')[0]
  np.savez_compressed( fstr, \
    xO=xO, yO=yO, zO=zO, xt=xt, yt=yt, zt=zt, ut=ut, vt=vt, wt=wt, dxyz=dxyz )
  print(' {}.npz saved successfully!'.format(fstr) )


# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def readNumpyZFootprintRaw( filename ):
  '''
  The saved .npz file contains 
  '''
  print(' Read raw footprint file {} ...'.format(filename))
  try: dat = np.load(filename)
  except: sys.exit(' Cannot read file {}. Exiting ...'.format(filename))
  
  xO = dat['xO']; yO = dat['yO']; zO = dat['zO']
  xt = dat['xt']; yt = dat['yt']; zt = dat['zt']
  ut = dat['ut']; vt = dat['vt']; wt = dat['wt'] 
  
  dat.close()
  
  return xO, yO, zO, xt, yt, zt, ut, vt, wt

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def writeNumpyZFootprint(filename, F, X, Y, Z, C, Ids=None ):
  fstr = filename.split('.npz')[0]
  if( Ids != None ):
    np.savez_compressed( fstr , F=F, X=X, Y=Y, Z=Z, C=C, Ids=Ids )
  else:
    np.savez_compressed( fstr , F=F, X=X, Y=Y, Z=Z, C=C )
  print(' {}.npz saved successfully!'.format(fstr))
  

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def readNumpyZFootprint( filename, IdsOn=False ):
  print(' Read footprint file {} ...'.format(filename))
  try: dat = np.load(filename)
  except: sys.exit(' Cannot read file {}. Exiting ...'.format(filename))
  
  F = dat['F']; X = dat['X']; Y = dat['Y']; Z = dat['Z']; C = dat['C']
  
  if( IdsOn ):
    try:    Ids = dat['Ids'].item()  # .item() returns the dict inside 0-array.
    except: Ids = None
    
  dat.close()
  
  if( IdsOn ):
    return F, X, Y, Z, C, Ids
  else:
    return F, X, Y, Z, C

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def fp2mshIJ(pxO, pyO, pzO, xG, yG, dx, dy ):  # IJ as in indecies.
  # Elegant and much faster. Use this!
  # pxO: particle x-origin, xG: x-grid coordinates.
  # First, Create meshgrid from the grid coords.
  X, Y = np.meshgrid( xG, yG )

  T  = np.zeros( np.shape(X) )  # Target
  Z  = np.zeros( np.shape(X) )  # Heights
  
  ix = ( pxO / dx ).astype(int); iy = ( pyO / dy ).astype(int)

  # The loop must be explicitly written open because 
  # the repeated additions to cells are not accounted properly.
  for i in range(len(ix)):
    T[iy[i],ix[i]] += 1.
  
  Z[iy[:],ix[:]] = pzO[:]
  
  return T, X, Y, Z
  

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


def fp2mshBM( pxO, pyO, pzO, xG, yG, dx, dy ):  # BM as in boolean matrix. 
  # Elegant routine, but awfully slow. Don't use this!
  # pxO: particle x-origin, xG: x-grid coordinates.
  # First, Create meshgrid from the grid coords.
  X, Y = np.meshgrid( xG, yG )
  
  # Then a mesh variable for storing the hits and the topography height.
  T = np.zeros( np.shape(X) )
  Z = np.zeros( np.shape(X) )
  
  for xi in xG:
    print(' Grid x-coord = {} '.format(xi))
    x1 = xi-dx/2.; x2 = xi+dx/2.
    PXb = ((x1 <= pxO) * (pxO < x2))
    if( PXb.any() ):
      for yi in yG:
        y1 = yi-dy/2.; y2 = yi+dy/2.
        # Utilizing the seeding coordinates (origin coords), extract how many hits each cell gets.
        PXYb = PXb * ((y1 <= pyO) * (pyO < y2))
        if( PXYb.any()):
          # Create a boolean matrix which isolates (with True value) the desired grid cell.
          MXb = ((x1 <= X)   * (X   < x2)  )
          MXYb = MXb * ((y1 <= Y)   * (Y < y2) )
        
          Z[MXYb] = np.mean( pzO[ PXYb ] )
          T += np.sum( PXYb.astype(int) ) * MXYb.astype(int)
  
  PXb = None; MXb = None; MXYb = None  # Clear
  return T, X, Y, Z

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def coordsFootprintGrid( NxG, dxG, pxO, pyO, verbose=False ):
  # Max values.
  xG_max = NxG[0]*dxG[0]  # Max dimensions.
  yG_max = NxG[1]*dxG[1]

  '''
  Note: At this point we assume non-cyclic boundary cond. for 
  the south/north boundaries. Therefore, the particles will be 
  absorbed if they come in contact with the y-normal boundaries.
  The footprint-grid will have to be extended backward only in
  the x-direction.
  '''

  # Smallest and largest x/y-value recorded:
  x_min = np.min( pxO ); y_min = np.min( pyO )
  x_max = np.max( pxO ); y_max = np.max( pyO )
  
  if(verbose):
    print( ' min(xO) = {}, max(xO) = {}'.format(x_min, x_max))
    print( ' min(yO) = {}, max(yO) = {}'.format(y_min, y_max))

  # Define an integer factor for domain multiplication/extension.
  fx = 0.
  if( x_min < 0. ):
    fx = int( abs(x_min) / xG_max ) + 1.
    
  # Coordinates for the extended footprint grid. Cell-centers.
  xD = np.linspace(-fx*xG_max+dxG[0]/2., xG_max-dxG[0]/2., (fx*NxG[0]+NxG[0])) # Last term: resolution.
  yD = np.linspace(dxG[1]/2.           , yG_max-dxG[1]/2.,  NxG[1]           )

  return xD, yD

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def idAppendices(fstring, ijkOn=False):
  if( ijkOn ):
    fileId = fstring.strip('.npz')                  # file ID string.
    fileId = fileId[-13:]
    varId  = fileId[-8:]; varId = varId.replace('.','_') # variable ID string.
  else:
    fileId = str()
    varId  = str(fn)

  return fileId, varId

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def percentileFootprintIds( F , p ):
  # 50, 75, 90
  p = p/100.
  Fsum = np.sum(F)
  Fpsum= p*Fsum
  
  fmax = np.max(F)  # maximum value.
  fv   = 0.5*fmax
  df   = fmax/350.  # values to increment.
  tol  = Fsum/2000.
  
  ic = 0
  while 1:
    ic += 1
    fv -= df
    idx = (F>fv)
    Fchecksum = np.sum(F[idx])
    
    if( (Fpsum-Fchecksum) < tol ):
      print(' i={}) TARGET vs. CURRENT: {} vs. {}'.format(ic,Fpsum,Fchecksum))
      break
    
  return idx

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def writeCrossWindSum( F , X, fname, idx=None ):
  import scipy.ndimage as sn # contains the filters
  
  nY, nX = np.shape(F)  # nRows, nCols
  Fm = np.zeros( nX )
  
  if( idx != None): Fx = F*idx
  else:             Fx = F.copy()
    
  for i in range( nX ):
    Fm[i] = np.sum(Fx[:,i])

  Fx  = None
  idx = (np.abs(Fm) > 0.)   # Select only non-zero entries
  Fm[idx] = sn.gaussian_filter( Fm[idx], sigma=2.5 )
  
  if( fname ):
    np.savetxt(fname+'_ysum.dat', np.c_[X[0,:],Fm] )   # x,y,z equal sized 1D arrays
    
  return Fm

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*= BEGIN KORMANN & MEIXNER =*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def kormann_and_meixner_fpr(z_0, z_m, u, sigma_v, L, X, Y, x_off=0., y_off=0. ):
  from scipy.optimize import fsolve
  from scipy.special import gamma
  
  Kappa = 0.41   # Von Karman const.
  
  # Bounds of integration for Eq. 42 to 46, defined on p.218
  z_1 = 3.*z_0
  z_2 = (1.+Kappa)*z_m
  
  # Data tuple for passing information to fsolve.
  data =(L, z_0, z_1, z_2, z_m)      
  
  
  # Final roots for m and n
  m0 = 0.5
  m  = fsolve( feqn_m, m0, args=data )[0]
  n  = fsolve( feqn_n, m0, args=data )[0]
  
  # Inversion of Eq 31
  u_star = u * Kappa / (np.log(z_m/z_0) + fopt1(L, z_m, z_m))
  
  # Eq (41), part 1
  U = u_star/Kappa * ( Iz_n(m   , L, z_0/z_m, z_1, z_2, z_m, 2 ) + \
      +                Iz_n(m   , L, z_0, z_1, z_2, z_m, 4, fopt1) ) \
      /              ( Iz_n(2.*m, L, z_0, z_1, z_2, z_m, 1 ) * z_m**m )
    
  # Eq (41), part 2
  K = Kappa*u_star * Iz_n(n,    L, z_0, z_1, z_2, z_m, 4, fopt2)\
    /              ( Iz_n(2.*n, L, z_0, z_1, z_2, z_m, 1 ) * z_m**(n-1.))
  
  # r is defined at the top of p.213, mu after Eq. 18
  r  = 2.+m-n
  mu = (1.+m)/r
  
  # Eq. 19
  xsi = U * z_m**r /( r**2 * K )
  
  # Eq. 21
  Xm    = np.abs(X-x_off)
  Ym    = np.abs(Y-y_off)
  Idm   = (X-x_off)>0.
  phi_x = ( gamma(mu)**(-1) * xsi**(mu)/( Xm**(1.+mu) ) * np.exp(-xsi/np.max(Xm,1e-10)) )* Idm
  
  # Cross wind diffusion
  # Eq. 18
  u_bar = gamma(mu)/gamma(1./r) * (r**2*K/U)**(m/r)*U*Xm**(m/r)
  
  # Eq. 9, definition of sig right after it
  sig = sigma_v*Xm/u_bar
  D_y = (np.sqrt(2.*np.pi)*sig)**(-1) * np.exp(-Ym**2./(2.*sig**2))
  
  phi = D_y * phi_x
  
  return phi[:,::-1]

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def fopt1(L, z, z_m=None):
  # This is used in eq 39 with J1 (Iz_4) and J2 (Iz_5). 
  if( L>0 ): 
    psi_m =  5.*z/L
  else:
    zeta = (1. - 16.*z/L)**(0.25)
    psi_m = (-2.)*np.log((1.+zeta)/2.) - np.log((1.+zeta**2)/2.) + 2.*np.arctan(zeta) - np.pi/2.
  
  return psi_m

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def fopt2(L, z, z_m):
  # This is used in eq 40 with J1 (Iz_4) and J2 (Iz_5).
  if( L>0 ): phi_c = 1. + 5.*z/L
  else:      phi_c = (1. - (16. * z/L))**(-0.5)
  
  rz = z/(phi_c * z_m)

  return rz

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
'''Following integrals (Eq 42-46) are solved numerically.
 They're all bundled within the same function to form a unified
 interface. This reduces code duplication. '''

def Iz_n(P, L, z_0, z_1, z_2, z_m, opt=1, fuser=None):

  az1 = (z_1/z_m); az2 = (z_2/z_m)
  dz = (az2-az1)/1000.
  az = np.arange(az1, az2, dz) + dz/2.

  
  if( opt == 1 ):    # I_1
    c = az**P * dz
  elif( opt == 2 ):  # I_2
    c = az**P * np.log(az/z_0) *dz
  elif( opt == 3 ):  # I_3
    c = az**P * np.log(az) * np.log(az/z_0) *dz
  elif( opt == 4 ):  # J_1
    c = az**P * fuser(L, az*z_m, z_m) * dz
  elif( opt == 5 ):  # J_2
    c = az**P * fuser(L, az*z_m, z_m)*np.log(az) * dz
    
  return np.sum(c)

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
def feqn_m( M, *data ):
  L, z_0, z_1, z_2, z_m = data
  A = Iz_n(2*M, L, z_0    , z_1, z_2, z_m, 1 ) * \
    ( Iz_n(  M, L, z_0/z_m, z_1, z_2, z_m, 3 ) + Iz_n(M, L, z_0, z_1, z_2, z_m, 5, fopt1) )
  B = Iz_n(2*M, L, 1      , z_1, z_2, z_m, 2 ) * \
    ( Iz_n(  M, L, z_0/z_m, z_1, z_2, z_m, 2 ) + Iz_n(M, L, z_0, z_1, z_2, z_m, 4, fopt1) )

  return (B - A)

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
def feqn_n( N, *data ):
  L, z_0, z_1, z_2, z_m  = data
  A = Iz_n(2*N, L, z_0, z_1, z_2, z_m, 1 ) * Iz_n(N, L, z_0, z_1, z_2, z_m, 5, fopt2)
  B = Iz_n(2*N, L, 1  , z_1, z_2, z_m, 2 ) * Iz_n(N, L, z_0, z_1, z_2, z_m, 4, fopt2)
  
  return (B - A)

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=  END KORMANN & MEIXNER  =*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=  BEGIN KLJUN  =*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def kljun_fpr(z_0, z_m, u_mean, sigma_v, L, Xt, Yt, z_i, us, x_off, y_off, nx=4000):
  rs=[1.]; wd = 0.; crop = True
  fdict = FFP(z_m, z_0, u_mean, z_i, L, sigma_v, us, None, rs, wd, nx, crop)
  fp = fdict['f_2d'].copy(); Xp = fdict['x_2d']; Yp = fdict['y_2d']
  
  dXt = Xt[0,2]-Xt[0,1]; dYt = Yt[2,0]-Yt[1,0]
  dXp = Xp[0,2]-Xt[0,1]; dYp = Yp[2,0]-Yt[1,0]
  
  ipt = (Xp[0,:]/dXt).astype(int); jpt = (Yp[:,0]/dYt).astype(int)
  print(' Xp = {} '.format(Xp[0,:]))
  print(' ipt = {} '.format(ipt))
  
  fdict = None
  
  # To be finalized ...
  
  return None  # Do not use yet


# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=  END KLJUN  =*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
