#!/usr/bin/env python
import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
from netcdfTools import read3dDataFromNetCDF
''' 
Description: 

Author:

'''
#==========================================================#
#==========================================================#
sepStr = ' # = # = # = # = # = # = # = # = '
parser = argparse.ArgumentParser()
parser.add_argument("-f","--filename", type=str, help="Name of input NETCDF data file.")
parser.add_argument("-v", "--varname",  type=str, default='u',\
  help="Name of the variable in NETCDF file. Default='u' ")
parser.add_argument("-p", "--printOn", action="store_true", default=False,\
  help="Print the numpy array data.") 
parser.add_argument("-pp", "--printOnly", action="store_true", default=False,\
  help="Only print the numpy array data. Don't save.")
parser.add_argument("-c", "--coarse", type=int, default=1,\
  help="Coarsening level. Int > 1.")
args = parser.parse_args()    
#==========================================================# 
# Rename ...
filename  = args.filename
varname   = args.varname
cl        = abs(args.coarse)
#==========================================================# 


# Read in data dictionary
dataDict = read3dDataFromNetCDF( filename , varname, cl )
vel = dataDict['v']
x   = dataDict['x']
y   = dataDict['y']
z   = dataDict['z']
time = dataDict['time']

infoStr = '''
  Coord. range:
  min(x)={0} ... max(x)={1}, nx = {2}
  min(y)={3} ... max(y)={4}, ny = {5}
  min(z)={6} ... max(z)={7}, nz = {8}
  min(time) = {9} ... max(time) = {10}, dt = {11}\n'''.format(\
    np.min(x), np.max(x), len(x),\
    np.min(y), np.max(y), len(y),\
    np.min(z), np.max(z), len(z), np.min(time), np.max(time), time[1]-time[0] )

print(infoStr)


print(' v_dims = {} '.format( np.shape(  vel[:,0,0,1]   ) ))


P = np.abs( np.fft.fft(vel) )**2

print(' Done ! ')

