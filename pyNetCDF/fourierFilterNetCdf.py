#!/usr/bin/env python3
from netcdfTools import *
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import rfft, irfft, fftfreq
from utilities import filesFromList, writeLog
''' 
Description:


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='fourierFilterNetCdf.py')
parser.add_argument("fileKey", default=None,\
  help="Search string for collecting files.")
parser.add_argument("-o", "--outstr",type=str, default="FFlt_",\
  help="Prefix for the output NETCDF file. Default=FFlt_.")
parser.add_argument("-vn", "--varnames",type=str, nargs='+', default=['u'],\
  help="Names of the V or V^prime comps in (x,y,z)-order. Default = ['u'].")
parser.add_argument("-vc", "--varcopy",type=str, nargs='+', default=None,\
  help="Names of the variables which copied to the output file without filtering.")
parser.add_argument("-lf", "--lowfreq", type=float, default=0.01,\
  help="Low frequency cutoff. FFT coefs will be zeroed for frequecies below this value.")
parser.add_argument("-c", "--coarse", type=int, default=1,\
  help="Coarsening level. Int > 1. Default = 1.")
args = parser.parse_args()
writeLog( parser, args )
#==========================================================#
# Initial renaming operations and variable declarations

fileKey    = args.fileKey 
outstr     = args.outstr
varnames   = args.varnames
varcopy    = args.varcopy
lowfreq    = args.lowfreq
cl         = abs(int(args.coarse))

'''
Establish two boolean variables which indicate whether the created variable is an
independent or dependent variable in function createNetcdfVariable().
'''

voDict    = dict()
parameter = True;  variable = False


# Obtain a list of files to include.
fileNos, fileList = filesFromList( fileKey+'*' )
for fn in fileNos:
  
  fileout = outstr+fileList[fn].split('_')[-1]
  parameter = True
  
  # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
  # Create a NETCDF output dataset (dso) for writing out the data.
  dso = netcdfOutputDataset( fileout )

  for vn in varnames:
    # Read in data.
    dataDict = read3dDataFromNetCDF( fileList[fn] , [vn], cl )
    vp = dataDict[vn]
    
    
    if( parameter ):
      # Coords and time:
      x  = dataDict.pop('x'); y = dataDict.pop('y'); z = dataDict.pop('z')
      time = dataDict.pop('time'); time_dim = len(time)
      
    dataDict = None


    if( parameter ):
      # Create the output independent variables right away and empty memory.
      tv = createNetcdfVariable( dso, time,'time', time_dim,'seconds','f4',('time',), parameter )

      xv = createNetcdfVariable( dso, x   , 'x'   , len(x)   , 'meters', 'f4', ('x',)   , parameter )
      x = None

      yv = createNetcdfVariable( dso, y   , 'y'   , len(y)   , 'meters', 'f4', ('y',)   , parameter )
      y = None

      zv = createNetcdfVariable( dso, z   , 'z'   , len(z)   , 'meters', 'f4', ('z',)   , parameter )
      z = None
      
      parameter = False
      
    # If our original signal time was in seconds, this is now in Hz    
    vfreq = fftfreq(len(vp[:,10,0,0]), d=time[1]-time[0])
    
    Fvp = rfft(vp, axis=(0))
    Fvp[(np.abs(vfreq)<lowfreq),:,:,:] = 0   # Filter step.
    vpf = irfft(Fvp, axis=(0)) + np.mean( vp, axis=(0) )
    Fvp = None 
    
    '''
    plt.figure(1)
    vm = np.mean( vp[:,10,0,0] )
    plt.plot(time,vp[:,10,0,0],'b', time, vm+vpf[:,10,0,0],'r')
    
    plt.figure(2)
    plt.semilogy(vfreq, Fvp[:,0,0,0]); plt.show()
    '''
    
    # Filtered value:
    voDict[vn] = createNetcdfVariable(\
      dso, vpf, vn, time_dim, '[-]', 'f4',('time','z','y','x',) , variable )


    # - - - - Done , finalize the output - - - - - - - - - -
  
  for vc in varcopy:
    dataDict = read3dDataFromNetCDF( fileList[fn] , [vc], cl )
    vpc = dataDict.pop(vc)
    if( len(np.shape( vpc )) == 4 ):
      voDict[vc] = createNetcdfVariable(dso, vpc, vc, time_dim,'[-]','f4',('time','z','y','x',), variable)
    elif( len(np.shape(vpc)) == 3 ):
      voDict[vc] = createNetcdfVariable(dso, vpc, vc, time_dim, '[-]', 'f4',('z','y','x',) , variable)
    else:
      print(' Unable to write {} into the output file. Skipping this step. '.format(vc))
      pass
    
  netcdfWriteAndClose( dso )
  
print(' Done! ')
