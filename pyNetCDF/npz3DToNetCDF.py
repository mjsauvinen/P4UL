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

Author: Jukka-Pekka Keskinen and Mikko Auvinen, FMI, 2020

'''
#==========================================================#
dMsg = '''Transform a 3D npz file to netCDF. 
The input file is required to have an array with the name S 
and coordinate system multipliers with the name dPx.
'''

parser = argparse.ArgumentParser(prog='npz3DToNetCDF.py', description=dMsg)
parser.add_argument("-f", "--filename", type=str, \
  help="Name of the input 3D npz file.")
parser.add_argument("-fo", "--fileout", type=str, default='output.nc',\
  help="Name of the output NetCDF file.")
parser.add_argument("-vn", "--varname", type=str, default='buildings_0',\
  help="Name of the variable in NetCDF. Default 'buildings_0'.")
parser.add_argument("-c", "--compress", action="store_true", default=False,\
  help="Compress netCDF variables with zlib.")
parser.add_argument("-i", "--integer", action="store_true", default=False,\
  help="Set datatype in NetCDF to integer. If not set, floats will be used.")
args = parser.parse_args()
writeLog( parser, args )
#==========================================================#
# Rename ... nothing more
filename = args.filename
fileout  = args.fileout.split('.nc')[0]+'.nc' # Ensure .nc is the suffix
integer  = args.integer

# = = = = = = = = = = = = = = = = = = = = = #
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
if(integer):
    netcdftype = 'i2'  # 16-bit signed integer
    #int32 = 'i4'  # 32-bit signed integer
    #int64 = 'i8'  # 64-bit signed integer
else:
    netcdftype = 'f4'  # 32-bit floating point
    #float64 = 'f8'  # 64-bit floating point
    #byte = 'b'  # One byte (8-bit)
#===============================================================#

# Load file
dat = np.load(filename)
A= dict(dat)
dat = None

# Get 3D scalar values (S)
S=A['S']
Rdims=S.shape  # Note, these [j,i,k] dimensions will be used below
Rdpx=A['dPx']

# Change order for NetCDF S[j,i,k] -> S[k,j,i]
S = np.rollaxis( S, 2, 0 )

# Get coordinate system dimensions and multipliers
dso = netcdfOutputDataset(fileout )
xv = createCoordinateAxis(dso, Rdims, Rdpx, 1, 'x', 'f4', 'm', parameter, args.compress)
yv = createCoordinateAxis(dso, Rdims, Rdpx, 0, 'y', 'f4', 'm', parameter, args.compress)
zv = createCoordinateAxis(dso, Rdims, Rdpx, 2, 'z', 'f4', 'm', parameter, args.compress)
# NOTE: j direction is mirrored as it, by default, advances in +N (i.e. -y)
tv = createNetcdfVariable(dso, S[:,::-1,:], args.varname, 0, 'm', netcdftype, ('z', 'y', 'x',), variable, False)

#tv.lod = 2

netcdfWriteAndClose(dso)
