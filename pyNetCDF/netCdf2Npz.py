#!/usr/bin/env python3
from netcdfTools import *
import sys
import argparse
import numpy as np
from utilities import partialMatchFromList
''' 
Description:


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='netCDF2Npz.py')
parser.add_argument("-f", "--filename",type=str, help="Name of the input NETCDF file.")
parser.add_argument("-fo", "--fileout",type=str, help="Name of the output .npz file.")
parser.add_argument("-m", "--mean", help="Extract mean values.", action="store_true", default=False) 
parser.add_argument("-c", "--coarse", type=int, help="Coarsening level. Int > 1.",\
  default=1) 
args = parser.parse_args() 
#==========================================================#
# Initial renaming operations and variable declarations

filename = args.filename
fileout  = args.fileout
meanOn   = args.mean
cl       = abs(int(args.coarse))

'''
Establish two boolean variables which indicate whether the created variable is an
independent or dependent variable in function createNetcdfVariable().
'''
parameter = True;  variable  = False

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
''' 
Create a NETCDF input dataset (ds), and its associated lists of dependent (varList)
and independent (dimList) variables. 
'''
ds, varList, paramList = netcdfDataset(filename)

'''
Read cell center coordinates and time.
Create the output independent variables right away and empty memory.

'''
uStr = partialMatchFromList( 'u', varList )
if( meanOn ):
  time = np.array([0,]); time_dims = np.shape(time)
else:
  time, time_dims = read1DVariableFromDataset('time', uStr, ds, 0, 0, 1 )

x, x_dims = read1DVariableFromDataset( 'x', uStr, ds, 0, 0, cl ) 
y, y_dims = read1DVariableFromDataset( 'y', uStr, ds, 0, 0, cl ) 
z, z_dims = read1DVariableFromDataset( 'z', uStr, ds, 0, 0, cl )

# - - - - velocity components - - - - - - - - - - - -
sa = ''
if( meanOn ): sa = 'm'
  
u, u_dims = read3DVariableFromDataset( 'u', ds, 0, 0, 0, cl, meanOn )
v, v_dims = read3DVariableFromDataset( 'v', ds, 0, 0, 0, cl, meanOn )
w, w_dims = read3DVariableFromDataset( 'w', ds, 0, 0, 0, cl, meanOn )

infoStr = '''
  time_dims = {}
  x_dims = {}
  y_dims = {}
  z_dims = {}
  u_dims = {}
  v_dims = {}
  w_dims = {}

'''.format( time_dims, x_dims, y_dims, z_dims, u_dims, v_dims, w_dims )

print(infoStr)

fstr = fileout.strip('.npz')+'.npz'
print(' Writing file {} ... '.format(fstr))
np.savez_compressed( fstr , time=time, u=u, v=v, w=w, x=x, y=y, z=z )
time = u = v = w = x = y = z = None
print(' ... done!')

