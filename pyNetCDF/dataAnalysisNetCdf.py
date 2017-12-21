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
parser.add_argument("-xn", "--xname",type=str, default='xu',\
  help="Specify the x coordinate. e.g. xu or x. Default='xu' ")
parser.add_argument("-yn", "--yname",type=str, default='y',\
  help="Specify the y coordinate. e.g. yv or y. Default='y' ")
parser.add_argument("-zn", "--zname",type=str, default='zu_3d',\
  help="Specify the z coordinate. e.g. zu_3d or zw_3d. Default='zu_3d' ")
parser.add_argument("-p", "--printOn", action="store_true", default=False,\
  help="Print the numpy array data.")
parser.add_argument("-pp", "--printOnly", action="store_true", default=False,\
  help="Only print the numpy array data. Don't save.")
parser.add_argument("-c", "--coarse", type=int, default=1,\
  help="Coarsening level. Int > 1.")
args = parser.parse_args()    
#==========================================================# 
# Rename ...
fileKey   = args.fileKey
normalize = args.normalize
mode      = args.mode
cl        = abs(args.coarse)
xname     = args.xname
yname     = args.yname
zname     = args.zname
varname   = args.varname


#==========================================================# 
# Create a dict that is passed into the function read3dDataFromNetCDF
nameDict = dict()
nameDict['xname']   = xname
nameDict['yname']   = yname
nameDict['zname']   = zname
nameDict['varname'] = varname


# Obtain a list of files to include.
fileNos, fileList = filesFromList( '*'+fileKey+'*' )

first = True
fig = plt.figure(num=1, figsize=(12,10))

for fn in fileNos:
  if('mag' in varname):
    nameDict['varname'] = 'u'
    nameDict['xname']   = 'x'; nameDict['yname'] = 'y'; nameDict['zname'] = 'z'
    dataDict = read3dDataFromNetCDF( fileList[fn] , nameDict, cl )
    u = dataDict['v']
    x = dataDict['x']; y = dataDict['y']; z = dataDict['z']
    nameDict['varname'] = 'v'
    dataDict = read3dDataFromNetCDF( fileList[fn] , nameDict, cl )
    v = dataDict['v']
    
    Umag = np.sqrt( u**2 + v**2 )
    vr = Umag
    
  else:
    dataDict = read3dDataFromNetCDF( fileList[fn] , nameDict, cl )
    vr = dataDict['v']
    x = dataDict['x']; y = dataDict['y']; z = dataDict['z']
    time = dataDict['time']
    
  dataDict = None
  
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
  
  fig = addToPlot(fig, vp, zp,'{}({}), {}'.format(mode,varname,fileList[fn]), plotStr, False )
  

plt.legend(loc=0)
plt.show()
