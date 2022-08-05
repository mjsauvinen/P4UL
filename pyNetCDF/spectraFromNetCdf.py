#!/usr/bin/env python3
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
Description: A script to perform spectral analysis on velocity data stored in a NETCDF file. 
This script adopts most of its content from old Matlab scripts that circled around the 
atmospheric science department.

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
#parser.add_argument("-T", "--DeltaT", help="Length of period to be analyzed in seconds [s].",\
#  type=float)
parser.add_argument("-nb", "--nbins", type=int, default=76,\
  help="Number of frequency bins. Default = 76")
parser.add_argument("-m", "--mode", type=str, default='S', choices=['S', 'E', 'P'],\
  help="Mode: 'S': power spectral density, 'E': energy spectrum, 'P': power spectrum.")
parser.add_argument("-n", "--normalize", action="store_true", default=False,\
  help="Compute f*S/sigma^2.")
parser.add_argument("-i1", "--ijk1",type=int, nargs=3,\
  help="Starting indices (ix, iy, iz) of the considered data. Required.")
parser.add_argument("-i2", "--ijk2",type=int, nargs=3,\
  help="Final indices (ix, iy, iz) of the considered data. Required.")
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
ijk1      = args.ijk1
ijk2      = args.ijk2
normalize = args.normalize
Nbins     = args.nbins
mode      = args.mode
cl        = abs(args.coarse)
varname   = args.varname

#==========================================================# 

# Obtain a list of files to include.
fileNos, fileList = filesFromList( fileKey+'*' )

first = True
fig   = None

for fn in fileNos:
  dataDict = read3dDataFromNetCDF( fileList[fn] , [varname] , cl )
  v = dataDict.pop(varname)
  x = dataDict.pop('x')
  y = dataDict.pop('y')
  z = dataDict.pop('z')
  time = dataDict.pop('time')
  
  if( first ):
    infoStr = '''
    Coord. range:
    min(x)={0} ... max(x)={1}, nx = {2}
    min(y)={3} ... max(y)={4}, ny = {5}
    min(z)={6} ... max(z)={7}, nz = {8}
    '''.format(\
      np.min(x), np.max(x), len(x),\
      np.min(y), np.max(y), len(y),\
      np.min(z), np.max(z), len(z) )
    print(infoStr)

    ijk1 = sensibleIds( np.array( ijk1 ), x, y, z )
    ijk2 = sensibleIds( np.array( ijk2 ), x, y, z )

    iList = np.arange(ijk1[0],ijk2[0]+1)
    jList = np.arange(ijk1[1],ijk2[1]+1)
    kList = np.arange(ijk1[2],ijk2[2]+1)
    try: 
      Np = int( input(" Number of plots per interval (empty -> all), Np = ") )
      stride = max( ((kList[-1]-kList[0])/Np)+1 , 2 )
    except:
      stride = 1
    first = False

  koff = groundOffset( v )
  if( koff > 0 ):
    print(' {}: koffset = {}'.format(fileList[fn], koff))
  
  
  for i in iList:
    for j in jList:
      for k in kList[::stride]:
        vt = v[:,k+koff,j,i]
        vName = varname+'(z={} m), {}'.format(z[k], fileList[fn].split('_NETCDF_')[-1])
        print(' Processing {} ...'.format(vName))
        fig = spectraAnalysis(fig, vt, time, vName, Nbins, mode, normalize)
        
  #print(' Failed to execute spectraAnalysis for {} ... '.format(fileList[fn]))
  #pass 

plt.legend(loc=0)
plt.show()
