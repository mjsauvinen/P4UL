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
def readVar( fn, vstr, cl=1 ):
  xDict = read3dDataFromNetCDF( fn , vstr , cl )
  v = xDict['v']; x = xDict['x']; y = xDict['y']; z = xDict['z']
  xDict = None
  return v, x, y, z
#==========================================================#
def U_hd( fn, cl=1, direction=False ):
  ut, xu, yu, zu = readVar( fn, 'u_xy', cl )
  vt, xv, yv, zv = readVar( fn, 'v_xy', cl )
  x = xv[:-1]; y = yu[:-1]; z = 0.5*(zu+zv)
  uc = 0.5*( ut[:,:,:-1,1:] + ut[:,:,:-1,:-1] )
  vc = 0.5*( vt[:,:,1:,:-1] + ut[:,:,:-1,:-1] )
  if( direction ):
    v = np.arctan( vc/(uc+1.E-4) ) * (180./np.pi)
  else:
    v = np.sqrt( uc**2 + vc**2 )
    
  return v, x, y, z

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
parser.add_argument("-nx", "--nexcl", type=int, default=0,\
  help="Exclude the first, [nx] number of nodes from analysis in x-direction.")
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
nx       = args.nexcl
writeRMS = args.writeRMS
printOn  = args.printOn
saveOn   = args.save
limsOn   = args.lims
gridOn   = args.grid

#----------------------------------------------------------#

# Shorter name
vn = varname.split('_')[0]
dirOn   = 'Ud' in varname
horizOn = 'Uh' in varname


if( (not horizOn) and (not dirOn) ):
  #print('{}'.format(varname))
  v1, x1, y1, z1 = readVar( f1, varname, 1 )
  v2, x2, y2, z2 = readVar( f2, varname, 1 )
  
else:
  v1, x1, y1, z1 = U_hd( f1, 1, dirOn )  
  v2, x2, y2, z2 = U_hd( f2, 1, dirOn )

dims1  = np.array( v1.shape )
dims2  = np.array( v2.shape )

if( not dirOn ):
  v1 -= v0[0]; v1 /= vs[0]
  v2 -= v0[1]; v2 /= vs[1]


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
  
  
  nx2 = 50
  v1x  = v1[0,k1,:,nx:-nx2]
  idx  = ( np.abs(v1x) > 1E-3 )
  vm1  = np.mean( v1x[idx] )
  f2 = 1.

  if( mode == 'r' ):
    dv = (f2 * v2[0,k2,:,nx:-nx2] - v1[0,k1,:,nx:-nx2])/np.abs( v1[0,k1,:,nx:-nx2] + 1E-5 )
  elif( mode == 's' ):
    dv = (f2 * v2[0,k2,:,nx:-nx2] - v1[0,k1,:,nx:-nx2])/( vm1 + 1E-5 )
  else:
    dv = (f2 * v2[0,k2,:,nx:-nx2] - v1[0,k1,:,nx:-nx2])


  N = float( np.prod( dv.shape ) )
  SkewDiff = (1./N)*np.sum(dv**3) * ( 1./(N-1.)*np.sum(dv**2) )**(-1.5) 
  RMSDiff  = np.sqrt(np.sum(dv**2)/N)
  print(' RMS (d{}) = {}, Sk(d{}) = {} '.format( vn , RMSDiff, vn, SkewDiff ))
  
  if( writeRMS ):
    fout.write('{:.2f}\t{:.2e}\n'.format( z1[k1], RMSDiff ))

  

  if( printOn ):
    dimsf  = np.array( np.shape( dv ) )
    xydims = dimsf
    figDims = 13.*(xydims[::-1].astype(float)/np.max(xydims))
    fig = plt.figure(num=1, figsize=figDims)
    labelStr = '({0}_2 - {0}_1)(z={1} m)'.format(vn, z1[k1])
    fig = addImagePlot( fig, dv[::-1,:], labelStr, gridOn, limsOn )
    
    
    if( saveOn ):
      figname = 'RMSDiff_{}_z{}.jpg'.format(vn, int(z1[k1]))
      print(' Saving = {}'.format(figname))
      fig.savefig( figname, format='jpg', dpi=150)
    plt.show()

if( writeRMS ): fout.close()