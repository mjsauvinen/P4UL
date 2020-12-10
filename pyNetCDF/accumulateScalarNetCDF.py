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

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

#==========================================================#
parser = argparse.ArgumentParser(prog='accumulateScalarNetCDF.py')
parser.add_argument("strKey", nargs='?', default=None,\
  help="Search string for collecting files.")
parser.add_argument("-s", "--scalar",type=str, required=True,\
  help="Name of the NETCDF scalar in the file. e.g. e or p, pt")
parser.add_argument("-c", "--coeff",type=float, default=1.0,\
  help="Multiplication coefficient for the chosen scalar. Default=1.0")
parser.add_argument("-pc", "--percentile", type=int, default=None,\
  help="Store given percentile instead. Overrides the standard functionality.")
parser.add_argument("-nt", "--ntimeskip", type=int, help="Skip <nt> number of time steps.",\
  default=0)
parser.add_argument("-zp", "--zprofile", action="store_true", default=False,\
  help="Compute vertical profile of final accumulation.")
parser.add_argument("-M", "--useMask", action="store_true", default=False,\
  help="Use (application specific) mask in computing the vertical profile (w/ --zprofile).")
args = parser.parse_args() 
#==========================================================#
# Initial renaming operations and variable declarations
strKey      = args.strKey
nt          = args.ntimeskip
scalarName  = args.scalar
coeff       = args.coeff
cl          = 1
zprofile    = args.zprofile
useMask     = args.useMask
pc          = args.percentile


parameter = True;  variable  = False
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
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

  
  if( pc is None ): fs = 'ACC_'
  else:             fs = 'P{:02d}_'.format(pc)
  
  dso = netcdfOutputDataset(fs+fileList[fn].split('/')[-1], mode="w")
  tv = createNetcdfVariable( dso, time,'time', ntimes,'s','f4',('time',), parameter )
  xv = createNetcdfVariable( dso, x, 'x', len(x), 'm', 'f4', ('x',), parameter )
  yv = createNetcdfVariable( dso, y, 'y', len(y), 'm', 'f4', ('y',), parameter )
  zv = createNetcdfVariable( dso, z, 'z', len(z), 'm', 'f4', ('z',), parameter )


  sa = np.zeros( np.shape(s) )
  
  if( pc is None ):
    # Work here 
    sx = 'a'
    #rStot = 1./np.sum( s[1,:,:,:] )
    rStot = 1.  #/np.percentile( s[1,:,:,:], 99 )
  
    sa[0,:,:,:] = (time[1]-time[0])*s[0,:,:,:]*rStot
    for i in range(1, ntimes ):
      dt = (time[i]-time[i-1])*rStot
      sa[i,:,:,:] = sa[i-1,:,:,:] + dt*s[i,:,:,:]
    
    
  else: # percentile
    sx = 'p'
    for i in range(1, ntimes ):
      idl = (s[i,:,:,:] > 1.E-5)
      Px  = np.percentile( s[i,idl], pc )
      idl = None
      
      idp = (s[i,:,:,:] >= Px)
      sa[i,idp] = s[i,idp]
    
  s = None # Clear memory
  
  if( zprofile ):
    saz = np.mean( sa[-1,:,:,:], axis=(1,2))  # compute all, and correct afterwards
    
    if( useMask ):
      nz, ny, nx = np.shape( sa[-1,:,:,:] )
      ids = np.ones( (ny,nx) , bool )
      ids[388:410,228:250] = False  # 1) UniqAir
      ids[197:218,248:271] = False  # 2) UniqAir
      for k in range(1,60):
        saz[k] = np.mean( sa[-1,k,ids] )
      ids = None

    print(' saz.shape = {}, saz = {} '.format(np.shape(saz), saz))
    hStr = ' z-coord (m), < sa >_xy '
    np.savetxt('horiz-mean_'+fileList[fn].split('/')[-1], np.c_[ z.ravel() , saz.ravel() ], fmt='%3.6e', header=hStr)
  

  sv = createNetcdfVariable(\
    dso, sa, scalarName+sx, s_dims[0], '', 'f4',('time','z','y','x',), variable)
  
  #hStr = '# time, {0}th percentile({1}), mean({1}), std({1}) '.format(percentile,sname)
  #np.savetxt(sname+'.dat', np.c_[ time, sp, sm, ss ], fmt='%3.6e', header=hStr)

  netcdfWriteAndClose( dso )
  tv = xv = yv = zv = sv = None # Clear everything


print(' Done !')
