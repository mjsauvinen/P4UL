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
  print ' Read raw footprint file {} ...'.format(filename)
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
  print ' Read footprint file {} ...'.format(filename)
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
  for i in xrange(len(ix)):
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
    
  for i in xrange( nX ):
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
def FFP(zm=None, z0=None, umean=None, h=None, ol=None, sigmav=None, ustar=None, 
        wind_dir=None, rs=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], rslayer=0,
        nx=1000, crop=False):
    """
    Derive a flux footprint estimate based on the simple parameterisation FFP
    See Kljun, N., P. Calanca, M.W. Rotach, H.P. Schmid, 2015: 
    The simple two-dimensional parameterisation for Flux Footprint Predictions FFP.
    Geosci. Model Dev. 8, 3695-3713, doi:10.5194/gmd-8-3695-2015, for details.
    contact: n.kljun@swansea.ac.uk

    FFP Input
    zm     = Measurement height above displacement height (i.e. z-d) [m]
    z0     = Roughness length [m]; enter None if not known 
    umean  = Mean wind speed at zm [m/s]; enter None if not known 
             Either z0 or umean is required. If both are given,
             z0 is selected to calculate the footprint
    h      = Boundary layer height [m]
    ol     = Obukhov length [m]
    sigmav = standard deviation of lateral velocity fluctuations [ms-1]
    ustar  = friction velocity [ms-1]

    optional inputs:
    wind_dir = wind direction in degrees (of 360) for rotation of the footprint    
    rs       = Percentage of source area for which to provide contours, must be between 10% and 90%. 
               Can be either a single value (e.g., "80") or a list of values (e.g., "[10, 20, 30]")
               Expressed either in percentages ("80") or as fractions of 1 ("0.8"). 
               Default is [10:10:80]. Set to "None" for no output of percentages
    nx       = Integer scalar defining the number of grid elements of the scaled footprint.
               Large nx results in higher spatial resolution and higher computing time.
               Default is 1000, nx must be >=600.
    rslayer  = Calculate footprint even if zm within roughness sublayer: set rslayer = 1
               Note that this only gives a rough estimate of the footprint as the model is not 
               valid within the roughness sublayer. Default is 0 (i.e. no footprint for within RS).
               z0 is needed for estimation of the RS.
    crop     = Crop output area to size of the 80% footprint or the largest r given if crop=1
 
    FFP output
    x_ci_max = x location of footprint peak (distance from measurement) [m]
    x_ci         = x array of crosswind integrated footprint [m]
    f_ci         = array with footprint function values of crosswind integrated footprint [m-1] 
    x_2d         = x-grid of 2-dimensional footprint [m], rotated if wind_dir is provided
    y_2d         = y-grid of 2-dimensional footprint [m], rotated if wind_dir is provided
    f_2d         = footprint function values of 2-dimensional footprint [m-2]
    rs       = percentage of footprint as in input, if provided
    fr       = footprint value at r, if r is provided
    xr       = x-array for contour line of r, if r is provided
    yr       = y-array for contour line of r, if r is provided
    flag_err = 0 if no error, 1 in case of error

    created: 15 April 2015 natascha kljun
    translated to python, December 2015 Gerardo Fratini, LI-COR Biosciences Inc.
    version: 1.22
    last change: 20/08/2016 natascha kljun
    Copyright (C) 2015, Natascha Kljun
    """
    import numpy as np
    import sys
    import numbers

    #===========================================================================
    ## Input check
    flag_err = 0
        
    ## Check existence of required input pars
    if None in [zm, h, ol, sigmav, ustar] or (z0 is None and umean is None):
        raise_ffp_exception(1)

    # Define rslayer if not passed
    if rslayer == None: rslayer == 0

    # Check passed values
    if zm <= 0.: raise_ffp_exception(2)
    if z0 is not None and umean is None and z0 <= 0.: raise_ffp_exception(3)
    if h <= 10.: raise_ffp_exception(4)
    if zm > h: raise_ffp_exception(5)        
    if z0 is not None and umean is None and zm <= 12.5*z0:
        if rslayer is 1: raise_ffp_exception(6)
        else: raise_ffp_exception(12)
    if float(zm)/ol <= -15.5: raise_ffp_exception(7)
    if sigmav <= 0: raise_ffp_exception(8)
    if ustar <= 0.1: raise_ffp_exception(9)
    if wind_dir is not None:
        if wind_dir> 360 or wind_dir < 0: raise_ffp_exception(10)
    if nx < 600: raise_ffp_exception(11)

    # Resolve ambiguity if both z0 and umean are passed (defaults to using z0)
    if None not in [z0, umean]: raise_ffp_exception(13)

    #===========================================================================
    # Handle rs
    if rs is not None:

        # Check that rs is a list, otherwise make it a list
        if isinstance(rs, numbers.Number): 
            if 0.9 < rs <= 1 or 90 < rs <= 100: rs = 0.9
            rs = [rs]
        if not isinstance(rs, list): raise_ffp_exception(14)

        # If rs is passed as percentages, normalize to fractions of one
        if np.max(rs) >= 1: rs = [x/100. for x in rs]

        # Eliminate any values beyond 0.9 (90%) and inform user
        if np.max(rs) > 0.9:
            raise_ffp_exception(15)
            rs = [item for item in rs if item <= 0.9]

        # Sort levels in ascending order
        rs = list(np.sort(rs))


    #===========================================================================
    # Model parameters
    a = 1.4524
    b = -1.9914
    c = 1.4622
    d = 0.1359
    ac = 2.17 
    bc = 1.66
    cc = 20.0

    #xstar_end = 30
    xstar_end  = 20  # Modified by Mikko to reduce the size of the domain.
    oln = 5000 #limit to L for neutral scaling
    k = 0.4 #von Karman

    #===========================================================================
    # Scaled X* for crosswind integrated footprint
    xstar_ci_param = np.linspace(d, xstar_end, nx+2)
    xstar_ci_param = xstar_ci_param[1:]

    # Crosswind integrated scaled F* 
    fstar_ci_param = a * (xstar_ci_param-d)**b * np.exp(-c/ (xstar_ci_param-d))
    ind_notnan     = ~np.isnan(fstar_ci_param)
    fstar_ci_param = fstar_ci_param[ind_notnan]
    xstar_ci_param = xstar_ci_param[ind_notnan]

    # Scaled sig_y*
    sigystar_param = ac * np.sqrt(bc * xstar_ci_param**2 / (1 + cc * xstar_ci_param))

    #===========================================================================
    # Real scale x and f_ci
    if z0 is not None:
        # Use z0
        if ol <= 0 or ol >= oln:
            xx  = (1. - 19.0 * zm/ol)**0.25
            psi_f = np.log((1. + xx**2)/2.) + 2.*np.log((1 + xx)/2.) - 2.*np.arctan(xx) + np.pi/2.
        elif ol > 0 and ol < oln:
            psi_f = -5.3 * zm / ol

        x = xstar_ci_param * zm / (1. - (zm/h)) * (np.log(zm/z0) - psi_f)
        if np.log(zm/z0) - psi_f > 0:
            x_ci = x
            f_ci = fstar_ci_param/zm * (1. - (zm/h)) / (np.log(zm/z0) - psi_f)
        else:
            x_ci_max, x_ci, f_ci, x_2d, y_2d, f_2d = None
            flag_err = 1
    else:
        # Use umean if z0 not available
        x = xstar_ci_param * zm/(1. - zm/h) * (umean/ustar * k)
        if umean / ustar > 0:
            x_ci = x
            f_ci = fstar_ci_param/zm * (1. - zm/h) / (umean/ustar * k)
        else:
            x_ci_max, x_ci, f_ci, x_2d, y_2d, f_2d = None
            flag_err = 1
                        
    #Maximum location of influence (peak location)
    xstarmax = -c / b + d
    if z0 is not None:
        x_ci_max = xstarmax * zm/(1. - (zm/h)) * (np.log(zm/z0) - psi_f)
    else:
        x_ci_max = xstarmax * zm/(1. - (zm/h)) * (umean/ustar * k)

    #Real scale sig_y
    if abs(ol) > oln:
        ol = -1E6
    if ol <= 0:   #convective
        scale_const = 1E-5 * abs(zm / ol)**(-1) + 0.80
    elif ol > 0:  #stable
        scale_const = 1E-5 * abs(zm / ol)**(-1) + 0.55
    if scale_const > 1:
        scale_const = 1.0
    sigy = sigystar_param/scale_const * zm * sigmav/ustar
    sigy[sigy < 0] = np.nan

    #Real scale f(x,y)
    dx = x_ci[2] - x_ci[1]
    y_pos = np.arange(0, (len(x_ci)/2.) * dx * 1.5, dx)
    #f_pos = np.full((len(f_ci), len(y_pos)), np.nan)
    f_pos = np.empty((len(f_ci), len(y_pos)))
    f_pos[:] = np.nan
    for ix in range(len(f_ci)):
        f_pos[ix,:] = f_ci[ix] * \
          1./(np.sqrt(2 * np.pi) * sigy[ix]) * np.exp(-y_pos**2/(2.* sigy[ix]**2))

    #Complete footprint for negative y (symmetrical)
    y_neg = - np.fliplr(y_pos[None, :])[0]
    f_neg = np.fliplr(f_pos)
    y = np.concatenate((y_neg[0:-1], y_pos))
    f = np.concatenate((f_neg[:, :-1].T, f_pos.T)).T

    #Matrices for output
    x_2d = np.tile(x[:,None], (1,len(y)))
    y_2d = np.tile(y.T,(len(x),1))
    f_2d = f
        

    #===========================================================================
    # Derive footprint ellipsoid incorporating R% of the flux, if requested,
    # starting at peak value.
    dy = dx
    if rs is not None:
        clevs = get_contour_levels(f_2d, dx, dy, rs)
        frs = [item[2] for item in clevs]
        xrs = []
        yrs = []
        for ix, fr in enumerate(frs):
            xr,yr = get_contour_vertices(x_2d, y_2d, f_2d, fr)
            if xr is None:
                frs[ix] = None
            xrs.append(xr)
            yrs.append(yr)
    else:
        if crop:
            rs_dummy = 0.8 #crop to 80%
            clevs = get_contour_levels(f_2d, dx, dy, rs_dummy)
            xrs = []
            yrs = []
            xrs,yrs = get_contour_vertices(x_2d, y_2d, f_2d, clevs[0][2])
                
    #===========================================================================
    # Crop domain and footprint to the largest rs value
    if crop:
        xrs_crop = [x for x in xrs if x is not None]
        yrs_crop = [x for x in yrs if x is not None]
        if rs is not None:
            dminx = np.floor(min(xrs_crop[-1]))
            dmaxx = np.ceil(max(xrs_crop[-1]))
            dminy = np.floor(min(yrs_crop[-1]))
            dmaxy = np.ceil(max(yrs_crop[-1]))
        else:
            dminx = np.floor(min(xrs_crop))
            dmaxx = np.ceil(max(xrs_crop))
            dminy = np.floor(min(yrs_crop))
            dmaxy = np.ceil(max(yrs_crop))
        jrange = np.where((y_2d[0] >= dminy) & (y_2d[0] <= dmaxy))[0]
        jrange = np.concatenate(([jrange[0]-1], jrange, [jrange[-1]+1]))
        jrange = jrange[np.where((jrange>=0) & (jrange<=y_2d.shape[0]-1))[0]]
        irange = np.where((x_2d[:,0] >= dminx) & (x_2d[:,0] <= dmaxx))[0]
        irange = np.concatenate(([irange[0]-1], irange, [irange[-1]+1]))
        irange = irange[np.where((irange>=0) & (irange<=x_2d.shape[1]-1))[0]]
        jrange = [[it] for it in jrange]
        x_2d = x_2d[irange,jrange]
        y_2d = y_2d[irange,jrange]
        f_2d = f_2d[irange,jrange]

    #===========================================================================
        #Rotate 3d footprint if requested
    if wind_dir is not None:                    
        wind_dir = wind_dir * np.pi / 180.
        dist = np.sqrt(x_2d**2 + y_2d**2)
        angle = np.arctan2(y_2d, x_2d)
        x_2d = dist * np.sin(wind_dir - angle)
        y_2d = dist * np.cos(wind_dir - angle)

        if rs is not None:
            for ix, r in enumerate(rs):
                xr_lev = np.array([x for x in xrs[ix] if x is not None])    
                yr_lev = np.array([x for x in yrs[ix] if x is not None])    
                dist = np.sqrt(xr_lev**2 + yr_lev**2)
                angle = np.arctan2(yr_lev,xr_lev)
                xr = dist * np.sin(wind_dir - angle)
                yr = dist * np.cos(wind_dir - angle)
                xrs[ix] = list(xr) 
                yrs[ix] = list(yr) 
                 
    #===========================================================================
    # Fill output structure
    if rs is not None:
        return {'x_ci_max': x_ci_max, 'x_ci': x_ci, 'f_ci': f_ci,
                'x_2d': x_2d, 'y_2d': y_2d, 'f_2d': f_2d,
                'rs': rs, 'fr': frs, 'xr': xrs, 'yr': yrs, 'flag_err':flag_err}
    else:
        return {'x_ci_max': x_ci_max, 'x_ci': x_ci, 'f_ci': f_ci,
                'x_2d': x_2d, 'y_2d': y_2d, 'f_2d': f_2d, 'flag_err':flag_err}

#===============================================================================
#===============================================================================
def get_contour_levels(f, dx, dy, rs=None):
    '''Contour levels of f at percentages of f-integral given by rs'''

    import numpy as np
    from numpy import ma
    import sys

    #Check input and resolve to default levels in needed
    if not isinstance(rs, (int, float, list)):
        rs = list(np.linspace(0.10, 0.90, 9))
    if isinstance(rs, (int, float)): rs = [rs]

    #Levels
    pclevs = np.empty(len(rs))
    pclevs[:] = np.nan
    ars = np.empty(len(rs))
    ars[:] = np.nan

    sf = np.sort(f, axis=None)[::-1]
    msf = ma.masked_array(sf, mask=(np.isnan(sf) | np.isinf(sf))) #Masked array for handling potential nan
        
    csf = msf.cumsum().filled(np.nan)*dx*dy
    for ix, r in enumerate(rs):
        dcsf = np.abs(csf - r)
        pclevs[ix] = sf[np.nanargmin(dcsf)]
        ars[ix] = csf[np.nanargmin(dcsf)]

    return [(round(r, 3), ar, pclev) for r, ar, pclev in zip(rs, ars, pclevs)]

#===============================================================================
def get_contour_vertices(x, y, f, lev):
    import matplotlib._cntr as cntr
        
    c = cntr.Cntr(x, y, f)
    nlist = c.trace(lev, lev, 0)
    segs = nlist[:len(nlist)//2]
    N = len(segs[0][:, 0])
    xr = [segs[0][ix, 0] for ix in range(N)] 
    yr = [segs[0][ix, 1] for ix in range(N)] 

    return [xr, yr]   # x,y coords of contour points.   

#===============================================================================
#===============================================================================
exTypes = {'message': 'Message',
           'alert': 'Alert',
           'error': 'Error',
           'fatal': 'Fatal error'}

exceptions = [
    {'code': 1,
     'type': exTypes['fatal'],
     'msg': 'At least one required parameter is missing. Please enter all '
            'required inputs. Check documentation for details.'},                                
    {'code': 2,
     'type': exTypes['fatal'],
     'msg': 'zm (measurement height) must be larger than zero.'},                                
    {'code': 3,
     'type': exTypes['fatal'],
     'msg': 'z0 (roughness length) must be larger than zero.'},                          
    {'code': 4,
     'type': exTypes['fatal'],
     'msg': 'h (BPL height) must be larger than 10m.'},                          
    {'code': 5,
     'type': exTypes['fatal'],
     'msg': 'zm (measurement height) must be smaller than h (PBL height).'},                             
    {'code': 6,
     'type': exTypes['alert'],
     'msg': 'zm (measurement height) should be above the roughness sub-layer (12.5*z0).'},
    {'code': 7,
     'type': exTypes['fatal'],
     'msg': 'zm/ol (measurement height to Obukhov length ratio) must be equal or larger than -15.5'},
    {'code': 8,
     'type': exTypes['fatal'],
     'msg': 'sigmav (standard deviation of crosswind) must be larger than zero'},
    {'code': 9,
     'type': exTypes['error'],
     'msg': 'ustar (friction velocity) must be >=0.1.'},
    {'code': 10,
     'type': exTypes['fatal'],
     'msg': 'wind_dir (wind direction) must be >=0 and <=360.'},
    {'code': 11,
     'type': exTypes['error'],
     'msg': 'nx must be >=600.'},
    {'code': 12,
     'type': exTypes['alert'],
     'msg': 'Using z0, ignoring umean.'},
    {'code': 13,
     'type': exTypes['error'],
     'msg': 'zm (measurement height) must be above roughness sub-layer (12.5*z0).'},
    {'code': 14,
     'type': exTypes['fatal'],
     'msg': 'if provided, rs must be in the form of a number or a list of numbers.'},
    {'code': 15,
     'type': exTypes['alert'],
     'msg': 'rs value(s) larger than 90% were found and eliminated.'},
        ]

def raise_ffp_exception(code):
    '''Raise exception or prints message according to specified code'''
        
    ex = [it for it in exceptions if it['code'] == code][0]
    string = ex['type'] + '(' + str(ex['code']).zfill(4) + '):\n '+ ex['msg'] 

    print('')
    if ex['type'] == exTypes['fatal']:
        string = string + '\n FFP_fixed_domain execution aborted.'
        raise Exception(string)
    else:
        print(string)



# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=  END KLJUN  =*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
