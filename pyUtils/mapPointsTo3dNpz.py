#!/usr/bin/env python3
import sys
import argparse
import numpy as np
from numTools import rotation_by_euler_angles
from netcdfTools import *
from scipy.ndimage.morphology import binary_closing

#==========================================================#
def centerCoords(xi,yi,zi, printOn=False):
  x2, x1 = np.max(x), np.min(x)
  y2, y1 = np.max(y), np.min(y)
  z2, z1 = np.max(z), np.min(z)
  if( printOn ):
    print('RANGES x:{0:.3f}/{1:.3f},  y:{2:.3f}/{3:.3f}, z:{4:.3f}/{5:.3f}'.format(x2,x1,y2,y1,z2,z1))
  
  xc, yc, zc = (x2+x1)*0.5, (y2+y1)*0.5, (z2+z1)*0.5
  
  if( printOn ):
    print('Center Coords: x,y,z = {0:.3f}, {1:.3f}, {2:.3f}'.format(xc,yc,zc))
  
  return xc, yc, zc

#==========================================================#

def coordsToIJK(x, y, z, Ix, Iy, Iz, dx, dy, dz, kZOn, d=None):
  
  if( d is None ): d = np.zeros(4)
  
  jt = (Ny-1)-(np.round(((y+d[2])/dy), decimals=0).astype(int) + Iy)
  it = np.round(((x+d[1])/dx), decimals=0).astype(int) + Ix
  kt = np.round(((z+d[3])/dz), decimals=0).astype(int)
  if( kZOn ): iCk = Iz - np.min(kt)
  kt += iCk
  
  return it, jt, kt 
#==========================================================#
def checkBounds( it, jt, kt, Nx, Ny, Nz ):
  # Check bounds
  jt = np.minimum( jt, Ny-1 ); jt = np.maximum( jt , 0 )
  it = np.minimum( it, Nx-1 ); it = np.maximum( it , 0 )
  kt = np.minimum( kt, Nz-1 ); kt = np.maximum( kt , 0 )

  return it, jt, kt
#==========================================================#
def mapPoints( Sx, it, jt, kt ):
  # Map values ... remember that j runs along Northing (i.e. negative y-direction)
  print('Map values onto 3D mesh ... ')
  for i,j,k in zip(it,jt,kt):
    Sx[k,j,i] = 1

  return Sx
#==========================================================#
parser = argparse.ArgumentParser(prog='mapPointsTo3dNpz.py')
parser.add_argument("-f", "--filename",type=str,required=True,\
  help="Name of the input .csv file containing points.")
parser.add_argument("-fo", "--fileout",type=str, default="out.npz", \
  help="Name of the output npz file.")
parser.add_argument("-fn", "--filenetcdf",type=str, default=None,\
  help="Name of the (optional) output NetCDF file.")
parser.add_argument("-s", "--scale",type=float, nargs=3, default=[1.,1.,1.],\
  help="Scale factors [sx sy sz] for the point coordinate values. Default=[1.,1.,1.] ")
parser.add_argument("-N", "--Nxyz",type=int, nargs=3, default=[None,None,None],\
  help="Dimensions of the npz block: Nx Ny Nz ")
parser.add_argument("-Dx", "--Dxyz",type=float, nargs=3, default=[None,None,None],\
  help="Resolution of the npz block: dx dy dz ")
parser.add_argument("-ea", "--eulerAngles",type=float, nargs=3, default=[None,None,None],\
  help="Euler angles [ea_z ea_y ea_x] for rotating the points around the center coord (deg).")
parser.add_argument("-Ic", "--centerPixel",type=int, nargs=3, default=[None,None,None],\
  help="Center pixel for the point placement: iE jN kZ ")
parser.add_argument("-IcF", "--centerPixelFile",type=str, default=None,\
  help="Name of file containing all center pixels [iE jN kZ] for the point placement. Overrides -Ic")
parser.add_argument("-d", "--dGapFill",type=float, default=0.,\
  help="Distance used for filling gaps. Default=None.")
parser.add_argument("-na", "--nans", action="store_true", default=False,\
  help="Initialize 3d numpy array with nans. By default, 3d array will be initialzed with zeros.")
parser.add_argument("-k0", "--kZeroHeight", action="store_true", default=False,\
  help="Position object vertically s.t. min(z)=0. Option -Ic kZ is applied afterwards.")

args = parser.parse_args()
#==========================================================#
# rename 
filename      = args.filename
fileout       = args.fileout
filenetcdf    = args.filenetcdf
sx, sy, sz    = np.array( args.scale )
kZero         = args.kZeroHeight
nansInit      = args.nans
Nx,Ny,Nz      = args.Nxyz
dx,dy,dz      = args.Dxyz
dm            = args.dGapFill
iCx, iCy, iCz = args.centerPixel # northing, easting, elevation
fileICntCoord = args.centerPixelFile

# - - - - - - - - - - #
if( dm == 0. ):
  dgf  = None   # delta gap fill 
else:
  dgf = np.zeros(4)

# - - - - - - - - - - #
if( iCx is not None ):
  if( (iCx>=Nx) or (iCy>=Ny) or (iCz>=Nz) ):
    sys.exit('Center pixel offset goes out of bounds. Exiting ...')
  else:
    iCx = [iCx]; iCy = [iCy]; iCz = [iCz] # convert to list

# - - - - - - - - - - #

if( ~np.any( np.array(args.eulerAngles) == None )):
  ea_z, ea_y, ea_x = np.array( args.eulerAngles ) * (np.pi/180.)
  ea_z = [ea_z]; ea_y = [ea_y]; ea_x = [ea_x]
else:
  ea_z = [0.]; ea_y = [0.]; ea_x = [0.]


# - - - - - - - - - - #
if( fileICntCoord is not None ):
  iCx, iCy, iCz, ea_z, ea_y, ea_x \
    = np.loadtxt( fileICntCoord, usecols=(0,1,2,3,4,5), unpack=True, dtype=float )
  
  ea_z *= (np.pi/180.); ea_y *= (np.pi/180.); ea_x *= (np.pi/180.)
  iCx   = iCx.astype(int); iCy = iCy.astype(int); iCz = iCz.astype(int)

nCntCoord = len(iCx)



#==========================================================#

print('Reading coordinates from {} ... '.format(filename))
x,y,z = np.loadtxt(filename, delimiter=',', unpack=True, skiprows=1)
print('... done!')

print('Determining center coordinate ')
xc, yc, zc = centerCoords(x,y,z, True)

# Apply scaling 
x *= sx; y *= sy; z *= sz
xc, yc, zc = centerCoords(x,y,z, True)
r1 = np.transpose( np.c_[(x-xc), (y-yc), (z-zc)] )

# 3d block indices for the point placements
# N/Y direction is special as Northing advances in negative y-direction.
c1 = np.ones(2); c1[1] = -1.

# Create the 3d block 
S = np.zeros((Nz, Ny, Nx), int )
if( nansInit ):
  S[:,:,:] = np.nan


for ic in range(nCntCoord):
  
  if( ea_z[ic] != 0. or ea_y[ic] != 0. or ea_x[ic] != 0. ):
    print('Performing 3D rotation ... ')
    rr = rotation_by_euler_angles( r1, [ea_z[ic], ea_y[ic], ea_x[ic]] )
    xp = rr[0,:]; yp = rr[1,:]; zp = rr[2,:]; rr = None
    print('... done!')
  else:
    xp = x; yp = y; zp = z
  
  
  print(' Mapping points to center coordinate (I,J,K) = ({},{},{}) ...'.format(iCx[ic], iCy[ic], iCz[ic]))

  if( dgf is None ):
    ia, ja, ka = coordsToIJK(xp, yp, zp, iCx[ic], iCy[ic], iCz[ic], dx, dy, dz, kZero, dgf)
    ia, ja, ka = checkBounds( ia, ja, ka, Nx, Ny, Nz )
    S = mapPoints( S, ia, ja, ka )
  else:
    for l in range(2):
      for m in range(4):
        if( m > 0 ): dgf[m] = dm * c1[l]
        print('Computing indices for delta={}'.format(dgf))
        ia, ja, ka = coordsToIJK(xp, yp, zp, iCx[ic], iCy[ic], iCz[ic], dx, dy, dz, kZero, dgf)
        ia, ja, ka = checkBounds( ia, ja, ka, Nx, Ny, Nz )
        S = mapPoints( S, ia, ja, ka )
        dgf[:] = 0.
  print('... done!')


print('Perform binary closing  ... ') 
for k in range(Nz):
  Sx = binary_closing( S[k,:,:] )
  S[k,:,:] = Sx[:,:]
print('... done! ')

# Make sure the flanks do not leak
S[:,0,:] = S[:,1,:]; S[:,-1,:] = S[:,-2,:]
S[:,:,0] = S[:,:,1]; S[:,:,-1] = S[:,:,-2]

xs = np.linspace(0., Nx*dx, Nx)
ys = np.linspace(0., Ny*dy, Ny)
zs = np.linspace(0., Nz*dz, Nz)

#print('xs: {}, ys: {}, zs: {}'.format(np.max(xs), np.max(ys), np.max(zs)))



if( filenetcdf is not None ):
  dso = netcdfOutputDataset(filenetcdf)
  parameter = True; variable  = False # For NetCDF output
  xv = createNetcdfVariable( dso, xs, 'x', len(xs), 'm', 'f4', ('x',), parameter )
  yv = createNetcdfVariable( dso, ys, 'y', len(ys), 'm', 'f4', ('y',), parameter )
  zv = createNetcdfVariable( dso, zs, 'z', len(zs), 'm', 'f4', ('z',), parameter )
  # NOTE: j direction needs to be mirrored.
  Sv = createNetcdfVariable(dso, S[:,::-1,:], 'S', 0, '-', 'i2', ('z','y','x'), False, False)  
  netcdfWriteAndClose(dso)


# Save the 3D npz array in S[j,i,k] order to maintain compatibility with R[j,i] rasters. 
# take first axis (k) and move it last to obtain (j,i,k) ordering: S[k,j,i] -> S[j,i,k]
S = np.rollaxis(S, 0, 3) 

print('Saving the npz dataset as a dictionary ...')
Sdict = dict()
Sdict['S'] = S
Sdict['Sdims'] = np.array( S.shape )
Sdict['GlobOrig']   = np.array( [   0. , 0., 0.] )    # Top left origo [N,E,Z]
Sdict['GlobOrigBL'] = np.array( [ Ny*dy, 0., 0.] )   # Bottom left origo [N,E,Z]
Sdict['dPx'] = np.array( [dy, dx, dz] )

fileout = fileout.split('.npz')[0]+'.npz'
np.savez_compressed(fileout, **Sdict)
print(' {} saved successfully!'.format(fileout))
