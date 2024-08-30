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

def write1dProfile( st, zt, time, Np, fn, excl0=False):
  dt = (len(time)//Np)-1
  tlist = time[::dt]
    
  ntimes  = len(tlist)
    
  for n in range(1,ntimes):
    saz = np.zeros( len(zt) )
    for k in range(len(zt)):
      if( excl0 ):
        ids = (st[-1,k,:,:]>0.)
        saz[k] = np.mean( st[n*dt,k,ids] )
      else:
        saz[k] = np.mean( st[n*dt,k,:,:] )
      
    #print(' saz.shape = {}, saz = {} '.format(np.shape(saz), saz))
    hStr = ' z-coord (m), < s >_xy '
    filename = fn.strip('.nc')+'_'+str(int(tlist[n]))+'.dat'
    print(' Writing {}'.format(filename))
    np.savetxt( filename , np.c_[ zt.ravel() , saz.ravel() ], fmt='%3.6e', header=hStr)
      

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

def writeTimeSeries( st, time, fn, excl0=False ):
  ntimes = len(time)
  savg = np.zeros( ntimes )
  if( excl0 ): 
    ids = (st[-1,:,:,:] > 0.) # Take the permanent zeros from the last time step.
    for i in range(ntimes):
      savg[i] = np.mean( st[i,ids] )
  else:
    for i in range(ntimes):
      savg[i] = np.mean( st[i,:,:,:] )
    
  hStr = ' time (s), < s >_xyz '
  np.savetxt( fn , np.c_[ time.ravel() , savg.ravel() ], fmt='%3.6e', header=hStr)
  savg = None


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#==========================================================#
parser = argparse.ArgumentParser(prog='scalarDistributionProfiles.py')
parser.add_argument("strKey", nargs='?', default=None,\
  help="Search string for collecting files.")
parser.add_argument("-s", "--scalar",type=str, required=True,\
  help="Name of the NETCDF scalar in the file. e.g. e or p, pt")
parser.add_argument("-c", "--coeff",type=float, default=1.0,\
  help="Multiplication coefficient for the chosen scalar. Default=1.0")
parser.add_argument("-nt", "--ntimeskip", type=int, default=0, \
  help="Skip <nt> number of time steps.")
parser.add_argument("-ts", "--timeSeries", action="store_true", default=False,\
  help="Compute time series of spatial mean values.")
parser.add_argument("-zp", "--zprofile", action="store_true", default=False,\
  help="Compute vertical profile of final accumulation.")
parser.add_argument("-xz", "--exclZeros", action="store_true", default=False,\
  help="Exclude permanent zeros from computations. Default=False")
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
xz          = args.exclZeros
useMask     = False
#useMask    = args.useMask

if( xz ): print(' NOTE: Permanent zeros will be excluded.')

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
    tsFileName = 's_ts_mean_'+fileList[fn].split('/')[-1]
    tsFileName = tsFileName.strip('.nc')+'.dat'
    writeTimeSeries(s, time, tsFileName, xz )

  if( zprofile ):
    zFileName = 's_xy-mean_'+fileList[fn].split('/')[-1]
    write1dProfile( s, z, time, 4, zFileName, xz )
    
print('Done!')
