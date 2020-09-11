#!/usr/bin/env python3

from netcdfTools import *
import argparse
import numpy as np
import sys

#==========================================================#
#parser = argparse.ArgumentParser(prog='rasterToNetCdf.py')
#parser.add_argument("-f", "--filename", type=str, help="Name of the input topography raster data file.")
#parser.add_argument("-fo", "--fileout", type=str, help="Name of the output NetCDF file.", default='output.ncdf')
#args = parser.parse_args()
#writeLog( parser, args )
#==========================================================#
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
#===============================================================#

#A=np.load(args.filename)

#print(A.files)

# Dims will be in the order y,x,z
A=np.random.randint(10,size=(2,3,5))
Rdims=A.shape

print(A)

dso = netcdfOutputDataset('testi.nc')

Rdpx=np.ones(3)

xv = createCoordinateAxis(dso, Rdims, Rdpx, 1, 'x', float32, 'm', parameter,False)#, args.compress)
yv = createCoordinateAxis(dso, Rdims, Rdpx, 0, 'y', float32, 'm', parameter,False)#, args.compress)
zv = createCoordinateAxis(dso, Rdims, Rdpx, 2, 'z', float32, 'm', parameter,False)#, args.compress)

#topo = fillTopographyArray(A, Rdims, Rdpx, int)
topovar = createNetcdfVariable(dso, A, 'testihomma', 0, 'm', int32, ('z', 'y', 'x',), variable, False)
topovar.lod = 2

netcdfWriteAndClose(dso)

