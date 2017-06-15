#!/usr/bin/env python
from netcdfTools import *
import sys
import argparse
import numpy as np
from mapTools import readNumpyZTile

'''
Description:
Reads topography data from Numpy Z file and exports it as NetCDF4.
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='topographyToNetCdf.py')
parser.add_argument("-f", "--filename",type=str, help="Name of the input topography raster data file.")
parser.add_argument("-fo", "--fileout",type=str, help="Name of the output NetCDF file.", default='output.ncdf')
parser.add_argument("-z", "--zlength",type=int, help="Number of grid points in z direction.")
parser.add_argument("-r", "--zresolution",type=int, help="Resolution of z axis. Defaults to resolution of N axis.")
args = parser.parse_args()
#==========================================================#

# Read input raster data file.
Rdict = readNumpyZTile(args.filename)
Rtopo = Rdict['R']
Rdims = np.shape(Rtopo)
Rdims = np.append(Rdims, args.zlength)
Rdpx = Rdict['dPx']

# Set z axis resolution
if (args.zresolution):
    Rdpx = np.append(Rdpx, args.zresolution)
else:
    Rdpx = np.append(Rdpx, Rdpx[0])

print(' Input raster data:')
print(' Size: [N,E] = [{}, {}]'.format(*Rdims))
print(' Resolution: [dPy,dPx] = [{}, {}] \n'.format(*Rdpx))
'''
Establish two boolean variables which indicate whether the created variable is an
independent or dependent variable in function createNetcdfVariable().
'''
parameter = True;  variable  = False

# Create netCDF output dataset
dso = netcdfOutputDataset( args.fileout )

'''
Create coordinate parameter arrays. These are 1D arrays containing
information on the position of the data point in metres.
'''
x = np.empty(Rdims[1]);
for i in xrange(Rdims[1]):
    x[i] = i*Rdpx[1] # dpx is in [N,E], see getGeoTransform() in gdalTools.py
xv = createNetcdfVariable(dso, x, 'x', len(x), 'm', 'f4', ('x',), parameter)
x = None

y = np.empty(Rdims[0])
for i in xrange(Rdims[0]):
    y[i] = i*Rdpx[0]
yv = createNetcdfVariable(dso, y, 'y', len(y), 'm', 'f4', ('y',), parameter)
y = None

z = np.empty(Rdims[2])
for i in xrange(Rdims[2]):
    z[i] = i*Rdpx[2]
zv = createNetcdfVariable(dso, z, 'z', len(z), 'm', 'f4', ('z',), parameter)
z = None

'''
Fill in a 3D array of topography data.
topo(z,y,x) containing 0 for air and 1 for land.
Loop through horizontal grid and use slices to fill the z grid.
'''
topodims = np.array([Rdims[2],Rdims[0],Rdims[1]])
topo = np.zeros(topodims, dtype=int)

print('\n Output data array:')
print(' Dimensions [z,y,x]: [{}, {}, {}]'.format(*topodims))
print(' Total number of data points: {}'.format(np.prod(topodims)))
print(' Filling the output array...')
for x in xrange(Rdims[1]):
    for y in xrange(Rdims[0]):
        maxind = int(round(Rtopo[y][x]))
        topo[0:maxind, y, x] = 1
print(' ...done. \n')

topovar = createNetcdfVariable( dso, topo, 'topography', 0, '', 'i4',('z','y','x',) , variable )
netcdfWriteAndClose(dso)
