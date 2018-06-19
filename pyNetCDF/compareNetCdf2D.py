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
parser.add_argument("-f1", "--filename1",type=str, help="Name of the first (ref) input NETCDF file.")
parser.add_argument("-f2", "--filename2",type=str, help="Name of the second input NETCDF file.")
parser.add_argument("-v", "--varname",  type=str, default='u',\
  help="Name of the variable in NETCDF file. Default='u' ")
parser.add_argument("-v0", "--vref", type=float, nargs=2, default=[0.,0.],\
  help="Reference values 'v0' in v+ = (v - v0)/v* for -f1 and -f2. Default = [0,0]")
parser.add_argument("-vs", "--vstar", type=float, nargs=2, default=[1.,1.],\
  help="Characteristic value 'v*' in v+ = (v - v0)/v* for -f1 and -f2. Default = [1,1]")
parser.add_argument("-m", "--mode", type=str, default='d', choices=['d', 'r', 's'],\
  help="Diff mode: 'd': delta, 'r': relative, 's': scaled.")
parser.add_argument("-w", "--writeRMS", help="Write the root-mean-square of the differences to a file.",\
  action="store_true", default=False)
parser.add_argument("-p", "--printOn", help="Print the numpy array data.",\
  action="store_true", default=False)
parser.add_argument("-s", "--save", action="store_true", default=False,\
  help="Save figures. Default=False")
parser.add_argument("--lims", help="User specified limits.", action="store_true", default=False)
parser.add_argument("--grid", help="Turn on grid.", action="store_true", default=False)
args = parser.parse_args()    

#==========================================================#
# Rename ... that's all.
f1       = args.filename1      # './DATA_2D_XY_AV_NETCDF_N02-1.nc'
f2       = args.filename2      # './DATA_2D_XY_AV_NETCDF_N02-2.nc'
varname  = args.varname
v0       = np.array(args.vref )
vs       = np.array(args.vstar)
mode     = args.mode
writeRMS = args.writeRMS
printOn  = args.printOn
saveOn   = args.save
limsOn   = args.lims
gridOn   = args.grid

#----------------------------------------------------------#

# Shorter name
vn = varname.split('_')[0]

d1Dict = read3dDataFromNetCDF( f1 , varname , 1 )
v1 = d1Dict['v']; x1 = d1Dict['x']; y1 = d1Dict['y']; z1 = d1Dict['z']
v1 -= v0[0]; v1 /= vs[0]
d1Dict = None 
dims1  = np.array( v1.shape )


d2Dict = read3dDataFromNetCDF( f2 , varname , 1 )
v2 = d2Dict['v']; x2 = d2Dict['x']; y2 = d2Dict['y']; z2 = d2Dict['z']
v2 -= v0[1]; v2 /= vs[1]
d2Dict = None 
dims2  = np.array( v2.shape )

if( all( dims1 == dims2 ) ):
  print(' Dimensions of the two datasets match!: dims = {}'.format(dims1))
else:
  print(' Caution! Dataset dimensions do not match. dims_1 = {} vs. dims_1 = {}'.format(dims1, dims2))


idk = selectFromList( z1 )


if( writeRMS ):
  fout = file('RMS_d{}.dat'.format(vn), 'wb')
  fout.write('# file1 = {}, file2 = {}\n'.format(f1, f2))
  fout.write('# z_coord \t RMS(d{})\n'.format(vn))
  #fout.write('{:.2f}\t{:.2e}'.format( z1[k1], dv ))

for k1 in idk:
  
  k2 = np.where(z2==z1[k1])[0] # This outputs a list 
  if( len(k2) == 0 ):
    print(' Coordinate {} not in file {}. Skipping.'.format(z1[k1],filename2))
    continue
  else:
    k2 = k2[0]    # Take always the first term
  
  
  #Correct the scales to avoid systematic differences
  #vm1 = np.mean( v1[0,k1,:,0] )  # Mean < >_y value at inlet 
  #vm2 = np.mean( v2[0,k2,:,0] )
  #f2  = vm1/vm2
  v1x  = v1[0,k1,:,:]
  idx  = ( np.abs(v1x) > 1E-3 )
  vm1  = np.mean( v1x[idx] )
  f2 = 1.

  if( mode == 'r' ):
    dv = (f2 * v2[0,k2,:,:] - v1[0,k1,:,:])/np.abs( v1[0,k1,:,:] + 1E-5 )
  elif( mode == 's' ):
    dv = (f2 * v2[0,k2,:,:] - v1[0,k1,:,:])/( vm1 + 1E-5 )
  else:
    dv = (f2 * v2[0,k2,:,:] - v1[0,k1,:,:])

  RMSDiff = np.sqrt(np.sum(dv**2)/float(np.prod(dv.shape)))
  print(' RMS (d{}) = {}'.format( vn , RMSDiff ))
  if( writeRMS ):
    fout.write('{:.2f}\t{:.2e}\n'.format( z1[k1], RMSDiff ))

  

  if( printOn ):
    xydims = dims1[2:]
    figDims = 13.*(xydims[::-1].astype(float)/np.max(xydims))
    fig = plt.figure(num=1, figsize=figDims)
    labelStr = '({0}_1 - {0}_2)(z={1} m)'.format(vn, z1[k1])
    fig = addImagePlot( fig, dv[::-1,:], labelStr, gridOn, limsOn )
    
    
    if( saveOn ):
      figname = 'RMSDiff_{}_z{}.jpg'.format(vn, int(z1[k1]))
      print(' Saving = {}'.format(figname))
      fig.savefig( figname, format='jpg', dpi=150)
    plt.show()

if( writeRMS ): fout.close()