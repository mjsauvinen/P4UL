#!/usr/bin/env python3
import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
from plotTools import addToPlot
from netcdfTools import netcdfDataset, readVariableFromDataset
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
parser.add_argument("-vn", "--vnames",  type=str, nargs=2, default=['u','v'],\
  help="Name of the horizontal velocity variables in the NETCDF file. Default=['u','v'] ")
parser.add_argument("-ptn", "--ptname",  type=str, default='pt',\
  help="Name of the potential temperature variable in the NETCDF file. Default='pt' ")
parser.add_argument("-c", "--coarse", type=int, default=1,\
  help="Coarsening level. Int > 1.")
args = parser.parse_args()    
#==========================================================# 
# Rename ...
fileKey   = args.fileKey
cl        = abs(args.coarse)
vnames    = args.vnames
ptname    = args.ptname

#==========================================================# 

# Obtain a list of files to include.
fileNos, fileList = filesFromList( fileKey+'*' )

fig = plt.figure(num=1, figsize=(12,10))
fig2= plt.figure(num=2, figsize=(12,10))
fig3= plt.figure(num=3, figsize=(12,10))

for fn in fileNos:
  ds, varList, paramList = netcdfDataset(fileList[fn], verbose=False)
  u1, veldict = readVariableFromDataset( vnames[0], ds, cl=1 )
  u2, veldict = readVariableFromDataset( vnames[1], ds, cl=1 )
  
  umag = np.sqrt( u1**2 + u2**2 ); u1 = None; u2 = None 

  pt, ptdict = readVariableFromDataset( ptname, ds, cl=1 )
  ptdict = None

  if( len(veldict.keys()) == 4 ):
    umag = np.mean( umag, axis=(0,2,3) )
    pt   = np.mean( pt  , axis=(0,2,3) )
  if( len(veldict.keys()) == 2 ):
    umag = np.mean( umag, axis=(0) )
    pt   = np.mean( pt  , axis=(0) )

  for dstr in veldict.keys():
    if( 'z' in dstr ): z = veldict[dstr]
  
  veldict = None
  
  dudz  = (umag[1:]-umag[:-1])/(z[1:]-z[:-1])
  dptdz = (pt[1:]-pt[:-1])/(z[1:]-z[:-1])
  zm    = 0.5*(z[1:]+z[:-1])
  
  Ri   = (9.81/np.mean(pt))*(dptdz)/(dudz**2+1e-9)  * ( dudz > 1e-3 ).astype(float)
  
  plotStr = ["Local Ri vs z ", "Ri" ,"z"]
  fig = addToPlot(fig, Ri[2:-4], zm[2:-4],'{}'.format(fileList[fn]), plotStr, False )
  
  plotStr = ["dpt/dz vs z ", "dpt/dz" ,"z"]
  fig2 = addToPlot(fig2, dptdz[2:-4], zm[2:-4],'{}'.format(fileList[fn]), plotStr, False )
  
  plotStr = ["du/dz vs z ", "du/dz" ,"z"]
  fig3 = addToPlot(fig3, dudz[2:-4], zm[2:-4],'{}'.format(fileList[fn]), plotStr, False )

plt.legend(loc=0)
plt.show()
