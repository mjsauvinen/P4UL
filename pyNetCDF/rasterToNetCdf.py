#!/usr/bin/env python
from netcdfTools import *
import sys
import argparse
import numpy as np
from mapTools import readNumpyZTile

'''
Description:
Reads data from Numpy npz file and exports it as an array in NetCDF.
Note that dealing with larger rasters require a lot of memory.
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='rasterToNetCdf.py')
parser.add_argument("-f", "--filename", type=str,
                    help="Name of the input topography raster data file.")
parser.add_argument("-fo", "--fileout", type=str,
                    help="Name of the output NetCDF file.", default='output.ncdf')
parser.add_argument("-N", "--NdZ", type=int,
                    help="Number of grid points in z direction. Leave empty to calculate automatically.")
parser.add_argument("-dz", "--dZ", type=float,
                    help="Resolution of z axis. Defaults to resolution of N axis.")
parser.add_argument("-vn", "--varname", type=str,
                    help="Name of the variable in NetCDF. Default 'topography'.", default='topography')
args = parser.parse_args()
#==========================================================#

# Read input raster data file.
Rdict = readNumpyZTile(args.filename)
Rtopo = Rdict['R']
Rdims = np.shape(Rtopo)
Rdpx = Rdict['dPx']

# Set z axis resolution
if (args.dZ):
    Rdpx = np.append(Rdpx, args.dZ)
else:
    Rdpx = np.append(Rdpx, Rdpx[0])

# Set vertical grid dimensions
if (args.NdZ):
    Rdims = np.append(Rdims, args.NdZ)
else:
    Rdims = np.append(Rdims, int(round(np.amax(Rtopo) / Rdpx[2])))

print(' Input raster data:')
print(' Size: [N,E] = [{}, {}]'.format(*Rdims))
print(' Resolution: [dPy,dPx] = [{}, {}] \n'.format(*Rdpx))
'''
Establish two boolean variables which indicate whether the created variable is an
independent or dependent variable in function createNetcdfVariable().
'''
parameter = True
variable = False

'''
Available external data types for NetCDF variables. Used data type has
a significant effect on file size and memory usage.
'''
int16 = 'i2'  # 16-bit signed integer
int32 = 'i4'  # 32-bit signed integer
int64 = 'i8'  # 64-bit signed integer
float32 = 'f4'  # 32-bit floating point
float64 = 'f8'  # 64-bit floating point
byte = 'b'  # One byte (8-bit)

'''
Create the dataset and coordinate parameter arrays. These are 1D
arrays containing information on the position of the data point in metres.
'''

dso = netcdfOutputDataset(args.fileout)
xv = createCoordinateAxis(dso, Rdims, Rdpx, 1, 'x', float32, 'm', parameter)
yv = createCoordinateAxis(dso, Rdims, Rdpx, 0, 'y', float32, 'm', parameter)
zv = createCoordinateAxis(dso, Rdims, Rdpx, 2, 'z', float32, 'm', parameter)

'''
Fill in a 3D array of topography data.
topo(z,y,x) containing 0 for air and 1 for land.
Loop through horizontal grid and use slices to fill the z grid.
'''
topo = fillTopographyArray(Rtopo, Rdims, Rdpx, bool)

topovar = createNetcdfVariable(
    dso, topo, args.varname, 0, '', int32, ('z', 'y', 'x',), variable)
netcdfWriteAndClose(dso)
