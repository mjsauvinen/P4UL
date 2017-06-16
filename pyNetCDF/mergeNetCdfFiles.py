#!/usr/bin/env python
from netcdfTools import *
import sys
import argparse
import numpy as np
from mapTools import readNumpyZTile

'''
Description:
Merges given NetCDF files together.
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='mergeNetCdfFiles.py')
parser.add_argument("-f", "--files", metavar='FILES', type=str, nargs='+', default=None,
                    help="List of netCDF files to be merged to output file.")
parser.add_argument("-fo", "--fileout", type=str,
                    help="Name of the output netCDF file. The script append to the file if it already exists.")
parser.add_argument("-c", "--compress", help="Compress netCDF variables with zlib.",
                    action="store_true", default=False)
args = parser.parse_args()
#==========================================================#

'''
Establish two boolean variables which indicate whether the created variable is an
independent or dependent variable in function createNetcdfVariable().
'''
parameter = True
variable = False

# Append or create new
dso = netcdfOutputDataset(args.fileout, 'w')

dsList = []
groupsList = []
for filename in args.files:
  # Create a data group for individual data sets
  print(' Processing file {}...'.format(filename))
  group = dso.createGroup("/" + filename.strip('.ncdf'))
  groupsList.append(group)
  ds, varList, paramList = netcdfDataset(filename)
  
  for v_name in varList:
    # For some reason parameters are also listed in variables (?)
    if v_name in paramList:
      vartype = parameter
    else:
      vartype = variable
    v_arr = ds.variables[v_name]
    vv = createNetcdfVariable(group, v_arr[:], v_name, len(
        v_arr), v_arr.units, v_arr.dtype, v_arr.dimensions, vartype, zlib=args.compress)
    v_arr = None
  ds = None
print(groupsList)
netcdfWriteAndClose(dso)
