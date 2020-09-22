#!/usr/bin/env python3

from netcdfTools import *
import argparse
import numpy as np
import sys
from utilities import writeLog
#=======================================================================
'''A script to transform 3D numpy array to netCDF. 

The 3D array is expected as the S file of the npz archive. The other
required file is dPx.

This could be merged at some point with rasterToNetCDF.py. 

Author: Jukka-Pekka Keskinen, FMI, 2020

'''
#==========================================================#
parser = argparse.ArgumentParser(prog='npz3DToNetCDF.py',description="Transform a 3D npz file to netCDF. The input file is "
"required to have an array with the name S and coordinate system multipliers with the name dPx.")
parser.add_argument("-f", "--filename", type=str, help="Name of the input 3D npz file.")
parser.add_argument("-fo", "--fileout", type=str, help="Name of the output NetCDF file.", default='output.ncdf')
parser.add_argument("-vn", "--varname", type=str, help="Name of the variable in NetCDF. Default 'buildings_0'.", 
                    default='buildings_0')
parser.add_argument("-c", "--compress", help="Compress netCDF variables with zlib.", action="store_true", default=False)
parser.add_argument("-i", "--integer", help="Set datatype in NetCDF to integer. If not set, floats will be used.", 
                    action="store_true", default=False)
args = parser.parse_args()
writeLog( parser, args )
#==========================================================#
'''Some settings for the netCDF outputs.

Establish two boolean variables which indicate whether the created variable is an
independent or dependent variable in function createNetcdfVariable().

'''
parameter = True
variable = False

'''
Available external data types for NetCDF variables. Used data type has
a significant effect on file size and memory usage.
'''
if args.integer:
    netcdftype = 'i2'  # 16-bit signed integer
    #int32 = 'i4'  # 32-bit signed integer
    #int64 = 'i8'  # 64-bit signed integer
else:
    netcdftype = 'f4'  # 32-bit floating point
    #float64 = 'f8'  # 64-bit floating point
    #byte = 'b'  # One byte (8-bit)
#===============================================================#

# Load file
dat = np.load(args.filename)
A= dict(dat)
dat = None

# Get 3D scalar values (S)
S=A['S']
Rdims=S.shape  # Note, these [j,i,k] dimensions will be used below
Rdpx=A['dPx']

# Change order for NetCDF S[j,i,k] -> S[k,j,i]
S = np.rollaxis( S, 2, 0 )

# Get coordinate system dimensions and multipliers
dso = netcdfOutputDataset(args.fileout)
xv = createCoordinateAxis(dso, Rdims, Rdpx, 1, 'x', 'f4', 'm', parameter, args.compress)
yv = createCoordinateAxis(dso, Rdims, Rdpx, 0, 'y', 'f4', 'm', parameter, args.compress)
zv = createCoordinateAxis(dso, Rdims, Rdpx, 2, 'z', 'f4', 'm', parameter, args.compress)

tv = createNetcdfVariable(dso, S, args.varname, 0, 'm', netcdftype, ('z', 'y', 'x',), variable, False)

#tv.lod = 2

netcdfWriteAndClose(dso)

