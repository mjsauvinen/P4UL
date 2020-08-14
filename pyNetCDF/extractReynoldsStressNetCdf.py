#!/usr/bin/env python3
from netcdfTools import *
import sys
import argparse
import numpy as np
from utilities import filesFromList, writeLog
''' 
Description:


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='extractReynoldsStressNetCdf.py')
parser.add_argument("fileKey", default=None,\
  help="Search string for collecting files.")
parser.add_argument("-o", "--outstr",type=str, default="RS",\
  help="Prefix for the output NETCDF file. Default=RS.")
parser.add_argument("-vn", "--vnames",type=str, nargs=3, default=['up','vp','wp'],\
  help="Names of the V or V^prime comps in (x,y,z)-order. Default = ['up','vp','wp'].")
parser.add_argument("-np", "--notPrimes", action="store_true", default=False,\
  help="Input data not as V^prime. They should be computed.")
parser.add_argument("-sn", "--sname",type=str, default=None,\
  help="Name of scalar s^prime. Default = None.")
parser.add_argument("-nt", "--ntimeskip", type=int, default=0,\
  help="Skip <nt> number of time steps. Default = 0.")
parser.add_argument("-c", "--coarse", type=int, default=1,\
  help="Coarsening level. Int > 1. Default = 1.")
args = parser.parse_args()
writeLog( parser, args )

#==========================================================#
# Initial renaming operations and variable declarations

fileKey    = args.fileKey 
outstr     = args.outstr
vnames     = args.vnames
sname      = args.sname
notPrimes  = args.notPrimes
nt         = args.ntimeskip
cl         = abs(int(args.coarse))

'''
Establish two boolean variables which indicate whether the created variable is an
independent or dependent variable in function createNetcdfVariable().
'''
parameter = True;  variable = False


# Obtain a list of files to include.
fileNos, fileList = filesFromList( fileKey+'*' )
for fn in fileNos:
  
  fileout = outstr+fileList[fn].split('_')[-1]
  
  # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
  # Read in data.
  dataDict = read3dDataFromNetCDF( fileList[fn] , vnames[0], cl )
  up = dataDict['v']
  
  # Commented out temporarily by Mikko #
  '''
  dataDict = read3dDataFromNetCDF( fileList[fn] , vnames[1], cl )
  vp = dataDict['v']
  '''
  dataDict = read3dDataFromNetCDF( fileList[fn] , vnames[2], cl )
  wp = dataDict['v']

  
  if( notPrimes ):
    # Perform coord. rotation for horizontal components
    um = np.mean( up, axis=(0) ); vm = np.mean( vp , axis=(0) )
    a  = np.arctan( vm/(um+1.e-5) )
    u1  = up * np.cos(a) + vp * np.sin(a)  # Streamwise comp.
    v1  =-up * np.sin(a) + vp * np.cos(a)  # Spanwise comp.
    up = u1; vp = v1
    
    up -= um
    vp -= vm
    wp -= np.mean( wp, axis=(0) )
  
  # Coords and time:
  x = dataDict['x']; y = dataDict['y']; z = dataDict['z']
  time = dataDict['time']; time_dim = len(time)
  dataDict = None

  # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
  # Create a NETCDF output dataset (dso) for writing out the data.
  dso = netcdfOutputDataset( fileout )

  # Create the output independent variables right away and empty memory.
  tv = createNetcdfVariable( dso, time,'time', time_dim,'s','f4',('time',), parameter )
  time = None  

  xv = createNetcdfVariable( dso, x   , 'x'   , len(x)   , 'm', 'f4', ('x',)   , parameter )
  x = None

  yv = createNetcdfVariable( dso, y   , 'y'   , len(y)   , 'm', 'f4', ('y',)   , parameter )
  y = None

  zv = createNetcdfVariable( dso, z   , 'z'   , len(z)   , 'm', 'f4', ('z',)   , parameter )
  z = None

  if( notPrimes ):
    # Streamwise velocity
    u1o = createNetcdfVariable(\
      dso, up, 'u1', time_dim, 'm s^(-1)', 'f4',('time','z','y','x',) , variable )
    u2o = createNetcdfVariable(\
      dso, vp, 'u2', time_dim, 'm s^(-1)', 'f4',('time','z','y','x',) , variable )
    
    u1mo = createNetcdfVariable(\
      dso, um, 'um1', time_dim, 'm s^(-1)', 'f4',('z','y','x',) , variable )
    u2mo = createNetcdfVariable(\
      dso, vm, 'um2', time_dim, 'm s^(-1)', 'f4',('z','y','x',) , variable )


  # First resolved TKE:
  tns = 0.5*( up**2 + wp**2 ) ## 0.5*( up**2 + vp**2 + wp**2 )
  tkeo = createNetcdfVariable(\
    dso, tns, 'tns', time_dim, 'm^2 s^(-2)', 'f4',('time','z','y','x',) , variable )
  tns = None

  tke = 0.5*( np.mean(up**2, axis=0) + np.mean(wp**2, axis=0) )  #0.5*( np.mean(up**2, axis=0) + np.mean(vp**2, axis=0) + np.mean(wp**2, axis=0) )
  rtke = createNetcdfVariable(\
    dso, tke, 'tke', time_dim, 'm^2 s^(-2)', 'f4',('z','y','x',) , variable )
  tke = None 

  # uu 
  uu  = up * up
  uuo = createNetcdfVariable(\
    dso, uu, 'cov_uu', time_dim, 'm^2 s^(-2)', 'f4',('time','z','y','x',) , variable )

  ruu = createNetcdfVariable(\
    dso, np.mean(uu, axis=0), 'r_uu', time_dim, 'm^2 s^(-2)', 'f4',('z','y','x',) , variable )
  uu  = None

  '''
  # uv 
  uv  = up * vp
  uvo = createNetcdfVariable(\
    dso, uv, 'cov_uv', time_dim, 'm^2 s^(-2)', 'f4',('time','z','y','x',) , variable )
  ruv = createNetcdfVariable(\
    dso, np.mean(uv, axis=0), 'r_uv', time_dim, 'm^2 s^(-2)', 'f4',('z','y','x',) , variable )
  uv  = None
  '''

  # uw 
  uw  = up * wp
  uwo = createNetcdfVariable(\
    dso, uw, 'cov_uw', time_dim, 'm^2 s^(-2)', 'f4',('time','z','y','x',) , variable )
  ruw = createNetcdfVariable(\
    dso, np.mean(uw, axis=0), 'r_uw', time_dim, 'm^2 s^(-2)', 'f4',('z','y','x',) , variable )
  uw  = None

  '''
  # vv 
  vv  = vp * vp
  vvo = createNetcdfVariable(\
    dso, vv, 'cov_vv', time_dim, 'm^2 s^(-2)', 'f4',('time','z','y','x',) , variable )
  rvv = createNetcdfVariable(\
    dso, np.mean(vv, axis=0), 'r_vv', time_dim, 'm^2 s^(-2)', 'f4',('z','y','x',) , variable )
  vv  = None

  vw  = vp * wp
  vwo = createNetcdfVariable(\
    dso, vw, 'cov_vw', time_dim, 'm^2 s^(-2)', 'f4',('time','z','y','x',) , variable )
  rvw = createNetcdfVariable(\
    dso, np.mean(vw, axis=0), 'r_vw', time_dim, 'm^2 s^(-2)', 'f4',('z','y','x',) , variable )
  vw  = None
  '''
  
  # ww 
  ww  = wp * wp
  wwo = createNetcdfVariable(\
    dso, ww, 'cov_ww', time_dim, 'm^2 s^(-2)', 'f4',('time','z','y','x',) , variable )
  rww = createNetcdfVariable(\
    dso, np.mean(ww, axis=0), 'r_ww', time_dim, 'm^2 s^(-2)', 'f4',('z','y','x',) , variable )
  ww  = None

  if( sname ):
    dataDict = read3dDataFromNetCDF( fileList[fn] , sname, cl )
    sp = dataDict['v']
    if( notPrimes ):
      sp -= np.mean( sp, axis=(0) )
    
    dataDict = None
    sstr = sname[:-1]   # Remove the 'p' (or similar) from the end.
    
    us = up * sp
    uso = createNetcdfVariable(\
      dso, us, 'cov_u'+sstr, time_dim, 'm s^(-1) []', 'f4',('time','z','y','x',), variable)
    rus = createNetcdfVariable(\
      dso, np.mean(us, axis=0), 'r_u'+sstr, time_dim, 'm s^(-1) []', 'f4',('z','y','x',), variable)
    us = None 
    
    vs = vp * sp
    vso = createNetcdfVariable(\
      dso, vs, 'cov_v'+sstr, time_dim, 'm s^(-1) []', 'f4',('time','z','y','x',), variable)
    rvs = createNetcdfVariable(\
      dso, np.mean(vs, axis=0), 'r_v'+sstr, time_dim, 'm s^(-1) []', 'f4',('z','y','x',), variable)
    vs = None
    
    ws = wp * sp
    wso = createNetcdfVariable(\
      dso, ws, 'cov_w'+sstr , time_dim, 'm s^(-1) []', 'f4',('time','z','y','x',), variable)
    rws = createNetcdfVariable(\
      dso, np.mean(ws, axis=0), 'r_w'+sstr, time_dim, 'm s^(-1) []', 'f4',('z','y','x',), variable)
    ws = None 


  # - - - - Done , finalize the output - - - - - - - - - -
  netcdfWriteAndClose( dso )
