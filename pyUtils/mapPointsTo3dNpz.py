#!/usr/bin/env python3
import sys
import argparse
import numpy as np
from numTools import rotation_by_euler_angles
from netcdfTools import *
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
parser.add_argument("-Ic", "--centerPixel",type=int, nargs=3, required=True,\
  help="Center pixel for the point placement: iE jN kZ ")
parser.add_argument("-na", "--nans", action="store_true", default=False,\
  help="Initialize 3d numpy array with nans. By default, 3d array will be initialzed with zeros.")
parser.add_argument("-k0", "--kZeroHeight", action="store_true", default=False,\
  help="Position object vertically s.t. min(z)=0. Option -Ic kZ is applied afterwards.")

args = parser.parse_args()
#==========================================================#
# rename 
filename     = args.filename
fileout      = args.fileout
filenetcdf   = args.filenetcdf
sx, sy, sz   = np.array( args.scale )
kZero        = args.kZeroHeight
nansInit     = args.nans
Nx,Ny,Nz     = args.Nxyz
dx,dy,dz     = args.Dxyz

eaOn = False
if( ~np.any( np.array(args.eulerAngles) == None )):
     eaOn = True
     ea_z, ea_y, ea_x = np.array( args.eulerAngles ) * (np.pi/180.)

iCx, iCy, iCz    = args.centerPixel # northing, easting, elevation
if( (iCx>=Nx) or (iCy>=Ny) or (iCz>=Nz) ):
  sys.exit('Center pixel offset goes out of bounds. Exiting ...')

#==========================================================#


x,y,z = np.loadtxt(filename, delimiter=',', unpack=True, skiprows=1)
xc, yc, zc = centerCoords(x,y,z, True)

# Rotate the x and y coords
if( eaOn ):
  r1 = np.transpose( np.c_[(x-xc), (y-yc), (z-zc)] )
  rr = rotation_by_euler_angles( r1, [ea_z, ea_y, ea_x] )
  x = rr[0,:]; y = rr[1,:]; z = rr[2,:]; rr = None; r1 = None

# Apply scaling 
x *= sx; y *= sy; z *= sz
xc, yc, zc = centerCoords(x,y,z, True)

# 3d block indices for the point placements
# N/Y direction is special as Northing advances in negative y-direction.
ja = (Ny-1)-(np.round((y/dy), decimals=0).astype(int) + iCy)
ia = np.round((x/dx), decimals=0).astype(int) + iCx
ka = np.round((z/dz), decimals=0).astype(int)
if( kZero ): iCz += -np.min(ka)
ka += iCz



# Check bounds
ja = np.minimum( ja, Ny-1 ); ja = np.maximum( ja , 0 )
ia = np.minimum( ia, Nx-1 ); ia = np.maximum( ia , 0 )
ka = np.minimum( ka, Nz-1 ); ka = np.maximum( ka , 0 )

# Create the 3d block 
S = np.zeros((Nz, Ny, Nx), int )
if( nansInit ):
  S[:,:,:] = np.nan
xs = np.linspace(0., Nx*dx, Nx)
ys = np.linspace(0., Ny*dy, Ny)
zs = np.linspace(0., Nz*dz, Nz)

#print('xs: {}, ys: {}, zs: {}'.format(np.max(xs), np.max(ys), np.max(zs)))

# Map values ... remember that j runs along Northing (i.e. negative y-direction)
for i,j,k in zip(ia,ja,ka):
  S[k,j,i] = 1


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

Sdict = dict()
Sdict['S'] = S
Sdict['Sdims'] = np.array( S.shape )
Sdict['GlobOrig']   = np.array( [   0. , 0., 0.] )    # Top left origo [N,E,Z]
Sdict['GlobOrigBL'] = np.array( [ Ny*dy, 0., 0.] )   # Bottom left origo [N,E,Z]
Sdict['dPx'] = np.array( [dy, dx, dz] )

fileout = fileout.split('.npz')[0]+'.npz'
np.savez_compressed(fileout, **Sdict)
print(' {} saved successfully!'.format(fileout))
