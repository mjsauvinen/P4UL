#!/usr/bin/env python3

import argparse
from netCDF4 import Dataset
import sys
import numpy as np

#==========================================================#
parser = argparse.ArgumentParser(prog='forceConstGridNetCdf.py')
parser.add_argument("-f", "--filename",type=str,\
  help="Name of the NETCDF file.")
parser.add_argument("-c", "--coords",type=str, nargs='+', required=True,\
  help="Coordinates to be force with a uniform grid spacing. Example: x y z")
parser.add_argument("-c0", "--c0values",type=float, nargs='+', required=True,\
  help="First coordinate values. Example: 0 0 -1")
args = parser.parse_args() 
#==========================================================#

filename = args.filename
cNames   = args.coords
c0Vals   = args.c0values

if( len( cNames ) != len( c0Vals ) ):
  sys.exit(" The number of entries for --coords and --x0values need to match. Exiting ...")

with Dataset(filename, "r+") as ds:
  for xn, x0 in zip( cNames, c0Vals ):
    x = ds.variables[xn][:]
    dx = np.round(  np.ma.median( np.diff(x) ) , decimals=2 )
    x_new = x0 + dx * np.arange(len(x))
    ds.variables[xn][:] = x_new

print(' Done ')
