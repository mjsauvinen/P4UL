#!/usr/bin/env python3
import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
from analysisTools import sensibleIds, groundOffset, quadrantAnalysis
from netcdfTools import read3dDataFromNetCDF, netcdfOutputDataset, \
  createNetcdfVariable, netcdfWriteAndClose
from utilities import filesFromList, inputIfNone
from txtTools import openIOFile
''' 
Description: Reynolds stress calculator. 
In case of PALM-generated results (featuring staggered grid), the velocity data must first be
interpolated onto cell-centers (i.e. scalar grid) with groupVectorDataNetCdf.py script.

Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''
#==========================================================#
sepStr = ' # = # = # = # = # = # = # = # = '
parser = argparse.ArgumentParser()
parser.add_argument("strKey", type=str,nargs='?', default=None,\
  help="Search string for collecting input NETCDF files.")
parser.add_argument("-v", "--varnames",  type=str, nargs=2, default=['u','w'],\
  help="Name of the variables in NETCDF file. Default=[u, w]")
parser.add_argument("-i1", "--ijk1",type=int, nargs=3,\
  help="Starting indices (ix, iy, iz) of the considered data. Required.")
parser.add_argument("-i2", "--ijk2",type=int, nargs=3,\
  help="Final indices (ix, iy, iz) of the considered data. Required.")
parser.add_argument("-vs", "--vstar",type=float, nargs=2, default=[1.,1.],\
  help="Characteristic value v* (vs) used in (v+ =(v-v0)/v*). Default=[1,1].")
parser.add_argument("-v0", "--vref",type=float, nargs=2, default=[0.,0.],\
  help="Reference value v0 (vref) used in (v+ =(v-v0)/v*). Default=[0,0].")
parser.add_argument("-xs", "--xscale",type=float, default=1.,\
  help="Coordinate scaling value (xs) used in (x+ =x/xs). Default=1.")
parser.add_argument("-of", "--outputToFile", type=str, default=None, \
  help="Name of the file to output analysis results. Default=None")
parser.add_argument("-p", "--printOn", action="store_true", default=False,\
  help="Print the numpy array data.")
args = parser.parse_args()    
#==========================================================# 
# Rename ...
strKey    = args.strKey
varnames  = args.varnames
v0        = np.array( args.vref  )  # Convert to numpy array
vs        = np.array( args.vstar )
xs        = args.xscale
ijk1      = args.ijk1
ijk2      = args.ijk2
printOn   = args.printOn
#==========================================================# 
'''
Establish two boolean variables which indicate whether the created variable is an
independent or dependent variable in function createNetcdfVariable().
'''
parameter = True;  variable  = False


strKey = inputIfNone( strKey , " Enter search string: " )
fileNos, fileList = filesFromList( strKey+"*")

for fn in fileNos:
  # - - - - - - - - - - - - - - - - - - - - - - - - - -  #
  # First fluctuation component
  cl = 1
  ncDict = read3dDataFromNetCDF( fileList[fn] , varnames[0], cl )
  v1 = ncDict['v']   # 'v' is a generic name for a variable in ncDict

  # Second fluctuation component
  ncDict = read3dDataFromNetCDF( fileList[fn] , varnames[1], cl )
  v2 = ncDict['v']


  # Dims
  nt, nz, ny, nx = np.shape( v1 )

  # - - - - - - - - - - - - - - - - - - - - - - - - - -  #
  # Spatial coords and time
  x = ncDict['x']; y = ncDict['y']; z = ncDict['z']
  time = ncDict['time']
  ncDict = None

  # Plot coord. information. This aids the user in the beginning.
  infoStr = '''
    Coord. range:
    min(x)={0} ... max(x)={1}, nx = {2}
    min(y)={3} ... max(y)={4}, ny = {5}
    min(z)={6} ... max(z)={7}, nz = {8}
  '''.format(\
      np.min(x), np.max(x), len(x),\
      np.min(y), np.max(y), len(y),\
      np.min(z), np.max(z), len(z) )
  #print(infoStr)

  # - - - - - - - - - - - - - - - - - - - - - - - - - -  #
  # Non-dimensionalize the time series
  v1 -= v0[0]; v2 -= v0[1]
  v1 /= vs[0]; v2 /= vs[1]

  # - - - - - - - - - - - - - - - - - - - - - - - - - -  #
  # Extract the fluctuations
  v1m = np.mean(v1, axis=(0))
  v2m = np.mean(v2, axis=(0))

  # Extract fluctuating part and normalize by variance
  # Reuse the v1 and  v2 variables to store values
  v1 -= v1m; v2 -= v2m

  # - - - - - - - - - - - - - - - - - - - - - - - - - -  #
  # Now check whether the given indices make sense 
  ijk1 = sensibleIds( np.array( ijk1 ), x, y, z )
  ijk2 = sensibleIds( np.array( ijk2 ), x, y, z )
  print(' Check (1): i, j, k = {}'.format(ijk1))
  print(' Check (2): i, j, k = {}'.format(ijk2))

  nvz = (ijk2[2]-ijk1[2])+1; idz = range(ijk1[2],ijk2[2]+1)
  nvy = (ijk2[1]-ijk1[1])+1; idy = range(ijk1[1],ijk2[1]+1)
  nvx = (ijk2[0]-ijk1[0])+1; idx = range(ijk1[0],ijk2[0]+1)
  Cv = np.zeros( ( nt, nvz, nvy, nvx ) )
  d  = np.zeros( (     nvz, nvy, nvx ) )
  # - - - - - - - - - - - - - - - - - - - - - - - - - -  #
  # Compute covariance
  for i in range(nvx):
    for j in range(nvy):
      for k in range(nvz):
        Cv[:,k,j,i] = v1[ :,idz[k], idy[j], idx[i] ] * v2[ :,idz[k], idy[j], idx[i] ]
        d[k,j,i] = np.sqrt( (z[idz[k]]-z[0])**2 + (y[idy[j]]-y[0])**2 + (x[idx[i]]-x[0])**2 )

  # - - - - - - - - - - - - - - - - - - - - - - - - - -  #
  # Reynolds stress
  Rs = np.mean( Cv, axis=(0) )

  hStr = " Reynolds averaged {}'{}' between ijk {} and {} ".format(varnames[0],varnames[1], ijk1,ijk2)
  fileout = '{}{}_UEX'.format(varnames[0],varnames[1]) + fileList[fn].split('UEX')[-1]
  fileout = fileout.strip('.nc') + '.dat'
  np.savetxt(fileout, np.c_[ (1./xs)*d.ravel(), Rs.ravel() ], fmt='%3.6e', header=hStr)
