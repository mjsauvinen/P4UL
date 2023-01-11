#!/usr/bin/env python3
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
parser.add_argument("-f", "--files", metavar='FILES', type=str, nargs='+', default=None, \
  help="List of netCDF files to be merged to output file.")
parser.add_argument("-fo", "--fileout", type=str,\
  help="Name of the output netCDF file. The script append to the file if it already exists.")
parser.add_argument("-c", "--compress", action="store_true", default=False,\
  help="Compress netCDF variables with zlib.")
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

# Read parameters from the first file
ds, varList, paramList = netcdfDataset(args.files[0], False)
print(' Reading parameters from {}...'.format(args.files[0]))
paramLengths = {}
matchGrid = ['x', 'y']
for pstr in paramList:
  pv = ds.variables[pstr]
  paramLengths[pstr] = len(pv)
  pv = createNetcdfVariable(\
    dso, pv[:], pstr, len(pv), pv.units, pv.dtype, pv.dimensions, parameter, zlib=args.compress )
  pv = None
ds.close()
print(' ...done.')

savedVars = []
for filename in args.files:
  # Create a data group for individual data sets
  print(' Processing file {}...'.format(filename))
  ds, varList, paramList = netcdfDataset(filename, False)

  # Check if x and y dimensions match.
  for pstr in matchGrid:
    pv = ds.variables[pstr]
    if (len(pv) != paramLengths[pstr]):
      sys.exit(' Error: Size mismatch in \'{}\' parameter'.format(pstr))

  # varList contains also the parameters
  varList = [x for x in varList if x not in paramList]
  for v_name in varList:
    if v_name in savedVars:
      sys.exit(' Error: Variable \'{}\' already saved.'.format(v_name))
    savedVars.append(v_name)
    vartype = variable
    v_arr = ds.variables[v_name]
    vv = createNetcdfVariable(dso, v_arr[:], v_name, len(v_arr), v_arr.units, v_arr.dtype, \
      v_arr.dimensions, vartype, zlib=args.compress)
    if (v_name=='buildings_0'):
      vv.lod=2
      print("vv.lod")
    v_arr = None
  ds.close()
netcdfWriteAndClose(dso)
