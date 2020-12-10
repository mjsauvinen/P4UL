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
    
    if( maskOn ):
      nz, ny, nx = np.shape( st[-1,:,:,:] )
      ids = np.ones( (ny,nx) , bool )
      ids[388:410,228:250] = False  # 1) UniqAir
      ids[197:218,248:271] = False  # 2) UniqAir
      for k in range(1,60):
        saz[k] = np.mean( st[-1,k,ids] )
      ids = None


    #print(' saz.shape = {}, saz = {} '.format(np.shape(saz), saz))
    hStr = ' z-coord (m), < s >_xy '
    np.savetxt( fn , np.c_[ zt.ravel() , saz.ravel() ], fmt='%3.6e', header=hStr)
    saz = None
    

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
parser.add_argument("-tx", "--timeExtension",type=float, nargs=2, default=[None,None],\
  help="Extend time to [timeX] with given increment [dtX] in seconds. Example usage: -tx timeX dtX ")
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
timex, dtx  = args.timeExtension

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
  
  # = = = = = = = = = = = = = = = = = = = = = = #
  sa = np.zeros( s_dims )
  
  if( pc is None ):
    # Work here 
    sffx = 'a'
    #rStot = 1./np.sum( s[1,:,:,:] )
    rStot = 1.  #/np.percentile( s[1,:,:,:], 99 )
  
    sa[0,:,:,:] = (time[1]-time[0])*s[0,:,:,:]*rStot
    for i in range(1, ntimes ):
      dt = (time[i]-time[i-1])*rStot
      sa[i,:,:,:] = sa[i-1,:,:,:] + dt*s[i,:,:,:]
    
    
  else: # percentile
    sffx = 'p'
    for i in range(1, ntimes ):
      idl = (s[i,:,:,:] > 1.E-5)
      Px  = np.percentile( s[i,idl], pc )
      idl = None
      
      idp = (s[i,:,:,:] >= Px)
      sa[i,idp] = s[i,idp]
      
  s1 = np.mean( s[-8:,:,:,:], axis=0 )   
  s  = None # Clear memory
  
  # = = = = = = = = = = = = = = = = = = = = = = #

  if( timex is not None  and dtx is not None ):
    ntimesx = (timex-time[-1])/dtx
    ntimesx = np.round( ntimesx, decimals=0).astype(int)  # Convert to int
    tx = np.linspace(time[-1], timex, ntimesx)
    
    sx_dims = np.array(s_dims); sx_dims[0] = ntimesx
    sx = np.zeros( sx_dims )
    sx[0,:,:,:] = sa[-1,:,:,:]*rStot
    for i in range(1,ntimesx):
      sx[i,:,:,:] = sx[i-1,:,:,:] + s1[:,:,:]*dtx*rStot
    
    s1 = None # Clear memory
    
    dsx = netcdfOutputDataset('ACCX_'+fileList[fn].split('/')[-1], mode="w")
    tvx = createNetcdfVariable( dsx, tx,'time', ntimesx,'s','f4',('time',), parameter )
    xvx = createNetcdfVariable( dsx, x, 'x', len(x), 'm', 'f4', ('x',), parameter )
    yvx = createNetcdfVariable( dsx, y, 'y', len(y), 'm', 'f4', ('y',), parameter )
    zvx = createNetcdfVariable( dsx, z, 'z', len(z), 'm', 'f4', ('z',), parameter )
    svx = createNetcdfVariable(\
      dsx, sx, scalarName+'x', sx_dims[0], '', 'f4',('time','z','y','x',), variable)
    netcdfWriteAndClose( dsx )
    
    if( zprofile ):
      write1dProfile( sx, z,'sx_horiz-mean_'+fileList[fn].split('/')[-1] , useMask ) 
    
    sx = None # Clear memory
    svx = zvx = yvx = xvx = tvx = dsx = None
  
  # = = = = = = = = = = = = = = = = = = = = = = #
  
  if( zprofile ):
    write1dProfile( sa, z,'sa_horiz-mean_'+fileList[fn].split('/')[-1] , useMask )    
  
  # = = = = = = = = = = = = = = = = = = = = = = #
  
  dso = netcdfOutputDataset(fs+fileList[fn].split('/')[-1], mode="w")
  tv = createNetcdfVariable( dso, time,'time', ntimes,'s','f4',('time',), parameter )
  xv = createNetcdfVariable( dso, x, 'x', len(x), 'm', 'f4', ('x',), parameter )
  yv = createNetcdfVariable( dso, y, 'y', len(y), 'm', 'f4', ('y',), parameter )
  zv = createNetcdfVariable( dso, z, 'z', len(z), 'm', 'f4', ('z',), parameter )
  sv = createNetcdfVariable(\
    dso, sa, scalarName+sffx, s_dims[0], '', 'f4',('time','z','y','x',), variable)
  
  sa = None
  #hStr = '# time, {0}th percentile({1}), mean({1}), std({1}) '.format(percentile,sname)
  #np.savetxt(sname+'.dat', np.c_[ time, sp, sm, ss ], fmt='%3.6e', header=hStr)

  netcdfWriteAndClose( dso )
  tv = xv = yv = zv = sv = None # Clear everything


print(' Done !')
