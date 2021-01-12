#!/usr/bin/env python3
from netcdfTools import *
import sys
import argparse
import numpy as np
from utilities import filesFromList, inputIfNone
''' 
Description:


Author: Mikko Auvinen
        mikko.auvinen@fmi.fi
        Finnish Meteorological Institute
'''
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

def readCoords_nc(ds, tskip, sStr):
  time, time_dims = read1DVariableFromDataset('time', sStr, ds, tskip, 0, 1 )
  x, x_dims = read1DVariableFromDataset( 'x', sStr,  ds, 0, 0, cl=1 )
  y, y_dims = read1DVariableFromDataset( 'y', sStr, ds, 0, 0, cl=1 )
  z, z_dims = read1DVariableFromDataset( 'z', sStr, ds, 0, 0, cl=1 )
  
  print(' dims(time,z,y,x) = ({},{},{},{})'.format(time_dims,z_dims,y_dims,x_dims))  

  return time,x,y,z

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

def write1dProfile( st, zt, fn, maskOn=False):
    saz = np.mean( st[-1,:,:,:], axis=(1,2))  # compute all, and correct afterwards

    #print(' saz.shape = {}, saz = {} '.format(np.shape(saz), saz))
    hStr = ' z-coord (m), < s >_xy '
    np.savetxt( fn , np.c_[ zt.ravel() , saz.ravel() ], fmt='%3.6e', header=hStr)
    saz = None
    

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

def writeTimeSeries( st, time, fn, maskOn=False ):
  ntimes = len(time)
  savg = np.zeros( ntimes )
  ids = (st[-1,:,:,:] > 0.) # Take the permanent zeros from the last time step.
  for i in range(ntimes):
    savg[i] = np.mean( st[i,ids] )
    
  hStr = ' time (s), < s >_xyz '
  np.savetxt( fn , np.c_[ time.ravel() , savg.ravel() ], fmt='%3.6e', header=hStr)
  savg = None


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#==========================================================#
parser = argparse.ArgumentParser(prog='accumulateScalarNetCDF.py')
parser.add_argument("strKey", nargs='?', default=None,\
  help="Search string for collecting files.")
parser.add_argument("-s", "--scalar",type=str, required=True,\
  help="Name of the NETCDF scalar in the file. e.g. e or p, pt")
parser.add_argument("-c", "--coeff",type=float, default=1.0,\
  help="Multiplication coefficient for the chosen scalar. Default=1.0")
parser.add_argument("-nt", "--ntimeskip", type=int, help="Skip <nt> number of time steps.",\
  default=0)
parser.add_argument("-ts", "--timeSeries", action="store_true", default=False,\
  help="Compute time series of spatial mean values.")
parser.add_argument("-zp", "--zprofile", action="store_true", default=False,\
  help="Compute vertical profile of final accumulation.")
args = parser.parse_args() 
#==========================================================#
# Initial renaming operations and variable declarations
strKey      = args.strKey
nt          = args.ntimeskip
scalarName  = args.scalar
zprofile    = args.zprofile
timeSeries  = args.timeSeries
coeff       = args.coeff
cl          = 1
useMask     = False
#useMask    = args.useMask

# - - - - Scalar components - - - - - - - - - -
strKey = inputIfNone( strKey , " Enter search string: " )
fileNos, fileList = filesFromList( strKey+"*" )

for fn in fileNos:
  ds, varList, paramList = netcdfDataset(fileList[fn])
  time,x,y,z = readCoords_nc(ds, nt, scalarName )
  ntimes = len(time)
  
  s, s_dims = read3DVariableFromDataset( scalarName, ds,  nt, 0, 0, cl ) # All values
  netcdfWriteAndClose( ds )
  
  idnan = np.isnan(s)
  if( np.count_nonzero(idnan) > 0 ): s[idnan] = 0.; idnan = None
    
  if( np.ma.isMaskedArray(s) ):
    idm = np.ma.getmask(s); print(' Nm = {}'.format(np.count_nonzero(idm)))
    s[idm] = 0.; idm = None
    s = np.array(s)
    
  if( timeSeries ):
    tsFileName = 'ts_mean'+fileList[fn].split('/')[-1]
    tsFileName = tsFileName.strip('.nc')+'.dat'
    writeTimeSeries(s, time, tsFileName, useMask )
