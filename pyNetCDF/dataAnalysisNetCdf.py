#!/usr/bin/env python
import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
from plotTools import addToPlot
from spectraTools import spectraAnalysis
from netcdfTools import read3dDataFromNetCDF
from analysisTools import sensibleIds, groundOffset, calc_ts_entropy_profile
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

def resample(X, n=None):
  nt, nk, nj, ni = X.shape
  Nn  = nt
  if(n is None): n = Nn
  else:          n = min( n, Nn )
  
  Xr = np.zeros( (n,nk,nj,ni) )
  #print(' Resampled array size = {} '.format( Xr.shape ))
  
  for i in xrange(ni):
    for j in xrange(nj):
      for k in xrange(nk):
        if(n is None): 
          n = Nn
        ix = np.floor(np.random.rand(n)*Nn).astype(int)
        Xr[:,k,j,i] = X[ix,k,j,i]
  
  return Xr

#==========================================================#


def calc_skew( V, axs=(0) ):
  import scipy.stats as st # contains st.entropy
  V -= np.mean( V, axis=(0) )
  vo = st.skew( V, axis=axs )
  
  return vo

#==========================================================#
sepStr = ' # = # = # = # = # = # = # = # = '
parser = argparse.ArgumentParser()
parser.add_argument("fileKey", default=None,\
  help="Search string for collecting files.")
parser.add_argument("-v", "--varname",  type=str, default='u',\
  help="Name of the variable in NETCDF file. Default='u' ")
parser.add_argument("-m", "--mode", type=str, default='mean', \
  choices=['mean', 'std', 'var','skew','entropy'],\
  help="Mode: mean, std, var, or entropy.")
parser.add_argument("-rs", "--resample", action="store_true", default=False,\
  help="Include resampling of the time series.")
parser.add_argument("-Ns", "--Nsample", type=int, default=None,\
  help="Number of terms used in resampling. Default=None")
parser.add_argument("-n", "--normalize", action="store_true", default=False,\
  help="Normalize.")
parser.add_argument("-me", "--meanError", action="store_true", default=False,\
  help="Plot std of mean error (with std option).")
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
fileKey     = args.fileKey
normalize   = args.normalize
meanErrorOn = args.meanError
mode        = args.mode
resampleOn  = args.resample
Ns          = args.Nsample
cl          = abs(args.coarse)
varname     = args.varname
writeAscii  = args.writeAscii


#==========================================================#


# Obtain a list of files to include.
fileNos, fileList = filesFromList( fileKey+'*' )

fig = plt.figure(num=1, figsize=(12,10))

for fn in fileNos:
  VNU = varname.upper()
  if('MAG' in VNU or 'U1' in VNU or 'U2' in VNU or 'DIR' in VNU):
    dataDict = read3dDataFromNetCDF( fileList[fn] , 'u', cl )
    u = dataDict['v']
    dataDict = read3dDataFromNetCDF( fileList[fn] , 'v', cl )
    v = dataDict['v']
    
    if('MAG' in VNU ): # vr := Umag
      vr = np.sqrt( u**2 + v**2 ); u = None; v = None
    else:
      um = np.mean( u, axis=(0) ); vm = np.mean( v, axis=(0) )
      a  = np.arctan( vm/(um+1.e-5) ); um = None; vm = None 
      if( 'U1' in VNU ):
        vr = u * np.cos(a) + v * np.sin(a)
      elif('U2' in VNU): # U2
        vr =-u * np.sin(a) + v * np.cos(a)
      else:# direction
        vr = np.arctan( v/(u+1.e-5) ) * (180./np.pi)
    
  elif('TKE' in VNU):
    try:
      dataDict = read3dDataFromNetCDF( fileList[fn] , 'e', cl )
      e_sgs = dataDict['v']
    except:
      print(' No e_sgs -> Result is RESOLVED TKE! ')
      e_sgs = None
    dataDict = read3dDataFromNetCDF( fileList[fn] , 'u', cl )
    u = dataDict['v']; up=u-np.mean(u, axis=0); u = None
    dataDict = read3dDataFromNetCDF( fileList[fn] , 'v', cl )
    v = dataDict['v']; vp=v-np.mean(v, axis=0); v = None
    dataDict = read3dDataFromNetCDF( fileList[fn] , 'w', cl )
    w = dataDict['v']; wp=w-np.mean(w, axis=0); w = None 
    
    e_res = 0.5*(np.mean(up**2,axis=0)+np.mean(vp**2,axis=0)+np.mean(wp**2,axis=0))
    up = None; vp = None; wp = None 
    
    # vr := TKE
    vr = e_res
    if( e_sgs is not None ): vr += e_sgs
    
  else:
    dataDict = read3dDataFromNetCDF( fileList[fn] , varname, cl )
    vr = dataDict['v']
  
  x  = dataDict['x']; y = dataDict['y']; z = dataDict['z']
  time = dataDict['time']
  dataDict = None
  axs = (0,2,3)
  axs = (0)

  if( resampleOn ):
    vr2 = resample( vr, Ns )
  else:
    vr2 = None

  
  # Process data vr --> vp 
  if( mode == 'mean'):
    vp  = np.mean( vr, axis=axs ); zp = z
    if( vr2 is not None ): vp2 = np.mean( vr2, axis=axs )
    plotStr  = ["mean({}) vs z ".format(varname), varname ,"z"]

  elif( mode == 'std'):
    vp  = np.std( vr, axis=axs ); zp = z
    if( vr2 is not None ): vp2 = np.std( vr2, axis=axs ) 
    
    if(meanErrorOn):
      N = len( vr[:,0,0,0] )
      vmerr = vp/np.sqrt(N)
      if( len(vmerr.shape)  == 3 ): vmerr  = vmerr[:,0,0]
      plotStr  = ["std. error of mean({}) vs z ".format(varname), varname ,"z"]
      fig = addToPlot(fig, vmerr, zp,'{}({}), {}'\
        .format('std error of mean',varname,fileList[fn].split('_')[-1]), plotStr, False )
    
      '''
      N2 = len( vr2[:,0,0,0] )
      vmerr2 = vp2/np.sqrt(N2)
      if( len(vmerr2.shape) == 3 ): vmerr2 = vmerr2[:,0,0]
      fig = addToPlot(fig, vmerr2, zp,'{}({}), {}'\
      .format('std error of mean',varname,fileList[fn]), plotStr, False )
      '''
    plotStr  = ["std({}) vs z ".format(varname), varname ,"z"]
    
  elif( mode == 'var' ):
    vp  = np.var( vr, axis=axs ); zp = z
    if( vr2 is not None ): vp2 = np.var( vr2, axis=axs )
    plotStr  = ["var({}) vs z ".format(varname), varname ,"z"]
  
  elif( mode == 'entropy' ):
    vp = calc_ts_entropy_profile( vr, z ); zp = z
    if( vr2 is not None ): vp2 = calc_ts_entropy_profile( vr2, z )
    plotStr  = ["entropy({}) vs z ".format(varname), varname ,"z"]
  
  elif( mode == 'skew' ):
    vp = calc_skew( vr, axs ); zp = z
    if( vr2 is not None ): vp2 = calc_skew( vr2, axs )
    plotStr  = ["skew({}) vs z ".format(varname), varname ,"z"]

# ================================================================= #

  if( len(vp.shape) == 3 ):  
    try: vp  = vp[:,1,1]
    except: vp = vp[:,0,0]

  if( writeAscii ):
    print(' (2) Writing data to ascii file: {}_{}.dat'.format(varname,mode))
    print(' x.shape = {} vs y.shape = {}'.format(np.shape(zp), np.shape(vp)))
    hStr = ' {} '.format(varname)
    fstr = fileList[fn].split('_')[-1]
    fstr = fstr.split('.')[0]
    np.savetxt(varname+'_'+mode+'_'+fstr+'.dat', np.c_[zp, vp], header=hStr)

  
  if( vr2 is not None ):
    if( len(vp2.shape) == 3 ): vp2 = vp2[:,1,1]
  
  fig = addToPlot(fig, vp,  zp,' {}({}), {}, N = {}'\
    .format(mode,varname,fileList[fn].split('_')[-1], len(time)), plotStr, False )
  
  if( vr2 is not None ):
    fig = addToPlot(fig, vp2, zp,' {}({}), {}, Resampled with N = {}'\
      .format(mode,varname,fileList[fn].split('_')[-1], Ns), plotStr, False )
    fig = addToPlot(fig, np.abs(vp-vp2), zp,' {}({}), {}, Resampling error = |(v_o-v_rs)/v_o|'\
      .format(mode,varname,fileList[fn].split('_')[-1]), plotStr, False )

plt.legend(loc=0)
plt.show()
