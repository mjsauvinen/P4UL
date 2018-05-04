#!/usr/bin/env python
import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
from plotTools import addToPlot
from spectraTools import spectraAnalysis
from netcdfTools import read3dDataFromNetCDF
from analysisTools import sensibleIds, groundOffset
from utilities import filesFromList
''' 
Description: A script to perform quadrant analysis on velocity data stored in a NETCDF file.
The analysis is performed for all points along a z-direction.
In case of PALM-generated results (featuring staggered grid), the velocity data must first be
interpolated onto cell-centers (i.e. scalar grid) with groupVectorDataNetCdf.py script.

Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''
#==========================================================#
#==========================================================#
sepStr = ' # = # = # = # = # = # = # = # = '
parser = argparse.ArgumentParser()
parser.add_argument("fileKey", default=None,\
  help="Search string for collecting files.")
parser.add_argument("-v", "--varname",  type=str, default='u',\
  help="Name of the variable in NETCDF file. Default='u' ")
parser.add_argument("-m", "--mode", type=str, default='mean', choices=['mean', 'std', 'var'],\
  help="Mode: mean, std, or var.")
parser.add_argument("-n", "--normalize", action="store_true", default=False,\
  help="Normalize.")
parser.add_argument("-p", "--printOn", action="store_true", default=False,\
  help="Print the numpy array data.")
parser.add_argument("-pp", "--printOnly", action="store_true", default=False,\
  help="Only print the numpy array data. Don't save.")
parser.add_argument("-wa", "--writeAscii", action="store_true", default=False,\
  help="Save profile data to an ascii file.")
parser.add_argument("-c", "--coarse", type=int, default=1,\
  help="Coarsening level. Int > 1.")
args = parser.parse_args()    
#==========================================================# 
# Rename ...
fileKey    = args.fileKey
normalize  = args.normalize
mode       = args.mode
cl         = abs(args.coarse)
varname    = args.varname
writeAscii = args.writeAscii

#==========================================================# 

# Obtain a list of files to include.
fileNos, fileList = filesFromList( fileKey+'*' )

fig = plt.figure(num=1, figsize=(12,10))

for fn in fileNos:
  if('mag' in varname):
    
    dataDict = read3dDataFromNetCDF( fileList[fn] , 'u', cl )
    u = dataDict['v']
    dataDict = read3dDataFromNetCDF( fileList[fn] , 'v', cl )
    v = dataDict['v']
    
    x = dataDict['x']; y = dataDict['y']; z = dataDict['z']
    
    # vr := Umag
    vr = np.sqrt( u**2 + v**2 )
    
  else:
    dataDict = read3dDataFromNetCDF( fileList[fn] , varname, cl )
    vr = dataDict['v']
    x  = dataDict['x']; y = dataDict['y']; z = dataDict['z']
    time = dataDict['time']
    
  dataDict = None
  
  
  # Process data vr --> vp 
  if( mode == 'mean'):
    vp = np.mean( vr, axis=(0,2,3) ); zp = z
    plotStr  = ["mean({}) vs z ".format(varname), varname ,"z"]
  elif( mode == 'std'):
    vp = np.std( vr, axis=(0,2,3) ); zp = z
    N = len( vr[:,0,0,0] )
    vmerr = vp/np.sqrt(N)
    plotStr  = ["std. error of mean({}) vs z ".format(varname), varname ,"z"]
    fig = addToPlot(fig, vmerr, zp,'{}({}), {}'.format('std error of mean',varname,fileList[fn]), plotStr, False )
    
    plotStr  = ["std({}) vs z ".format(varname), varname ,"z"]
  elif( mode == 'var' ):
    vp = np.var( vr, axis=(0,2,3) ); zp = z
    plotStr  = ["var({}) vs z ".format(varname), varname ,"z"]
  

  if( writeAscii ):
    print(' (2) Writing data to ascii file: {}.dat'.format(varname))
    print(' x.shape = {} vs y.shape = {}'.format(np.shape(zp), np.shape(vp)))
    hStr = ' {} '.format(varname)
    fstr = fileList[fn].split('_')[-1]
    fstr = fstr.split('.')[0]
    np.savetxt(varname+'_'+mode+'_'+fstr+'.dat', np.c_[zp, vp], header=hStr)


  fig = addToPlot(fig, vp, zp,'{}({}), {}'.format(mode,varname,fileList[fn]), plotStr, False )
  

plt.legend(loc=0)
plt.show()
