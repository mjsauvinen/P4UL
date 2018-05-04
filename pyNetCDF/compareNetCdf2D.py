#!/usr/bin/env python

import netCDF4 as nc
import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse
import scipy.ndimage as sn # contains the filters
from plotTools import addImagePlot
from netcdfTools import read3dDataFromNetCDF
from utilities import selectFromList

#==========================================================#

parser = argparse.ArgumentParser(prog='compareNetCdf2D.py')
parser.add_argument("-f1", "--filename1",type=str, help="Name of the first (1) input NETCDF file.")
parser.add_argument("-f2", "--filename2",type=str, help="Name of the second (2) input NETCDF file.")
parser.add_argument("-v", "--varname",  type=str, default='u',\
  help="Name of the variable in NETCDF file. Default='u' ")
parser.add_argument("-r", "--relative", help="Compare relative differences: (v1-v2)/|v1|.",\
  action="store_true", default=False)
parser.add_argument("-p", "--printOn", help="Print the numpy array data.",\
  action="store_true", default=False)
parser.add_argument("-pp", "--printOnly", action="store_true", default=False,\
  help="Only print the numpy array data. Don't save.")
parser.add_argument("--lims", help="User specified limits.", action="store_true", default=False)
parser.add_argument("--grid", help="Turn on grid.", action="store_true", default=False)
args = parser.parse_args()    

#==========================================================#
# Rename ... that's all.
f1       = args.filename1      # './DATA_2D_XY_AV_NETCDF_N02-1.nc'
f2       = args.filename2      # './DATA_2D_XY_AV_NETCDF_N02-2.nc'
varname  = args.varname
relative = args.relative
printOn  = args.printOn
printOnly= args.printOnly
limsOn   = args.lims
gridOn   = args.grid

#----------------------------------------------------------#

d1Dict = read3dDataFromNetCDF( f1 , varname , 1 )
v1 = d1Dict['v']; x1 = d1Dict['x']; y1 = d1Dict['y']; z1 = d1Dict['z']
d1Dict = None 
dims1  = np.array( v1.shape )


d2Dict = read3dDataFromNetCDF( f2 , varname , 1 )
v2 = d2Dict['v']; x2 = d2Dict['x']; y2 = d2Dict['y']; z2 = d2Dict['z']
d2Dict = None 
dims2  = np.array( v2.shape )

if( all( dims1 == dims2 ) ):
  print(' Dimensions of the two datasets match!: dims = {}'.format(dims1))
else:
  print(' Caution! Dataset dimensions do not match. dims_1 = {} vs. dims_1 = {}'.format(dims1, dims2))


idk = selectFromList( z1 )

for k1 in idk:
  k2 = np.where(z2==z1[k1])[0] # This outputs a list 
  k2 = k2[0]                  # Take always the first term
  
  #Correct the scales to avoid systematic differences
  vm1 = np.mean( v1[0,k1,:,0] )  # Mean < >_y value at inlet 
  vm2 = np.mean( v2[0,k2,:,0] )
  f2  = vm1/vm2 

  if( relative ):
    dv = (v1[0,k1,:,:] - f2 * v2[0,k2,:,:])/np.abs( v1[0,k1,:,:] )
  else:
    dv = (v1[0,k1,:,:] - f2 * v2[0,k2,:,:])


  if( printOn or printOnly ):
    xydims = dims1[2:]
    figDims = 13.*(xydims[::-1].astype(float)/np.max(xydims))
    fig = plt.figure(num=1, figsize=figDims)
    fig = addImagePlot( fig, dv[:,:], 'dv', gridOn, limsOn )
    plt.show()
