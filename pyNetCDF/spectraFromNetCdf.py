#!/usr/bin/env python
import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
from plotTools import addToPlot
from spectraTools import spectraAnalysis
from netcdfTools import read3dDataFromNetCDF

#==========================================================#
def sensibleIds( ixyz, x, y, z ):
  ixyz[0] = np.minimum( ixyz[0] , len(x)-1 ); ixyz[0] = np.maximum( ixyz[0], 0 )
  ixyz[1] = np.minimum( ixyz[1] , len(y)-1 ); ixyz[1] = np.maximum( ixyz[1], 0 )
  ixyz[2] = np.minimum( ixyz[2] , len(z)-1 ); ixyz[2] = np.maximum( ixyz[2], 0 )
  
  return ixyz
#==========================================================#

'''
Kaimal & Finnigan:
The requirements for averaging time T with T >> Tau_{\alpha} can then be expressed
in terms of \sigma^{2}_{\bar{\alpha}}, the variance of the measured time mean \bar{\alpha}
about the expected ensemple mean, and \sigma^{2}_{\alpha}, the ensemble variance of \alpha.
'''
#==========================================================#
sepStr = ' # = # = # = # = # = # = # = # = '
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filename", type=str,\
  help="Name of the input NETCDF file.")
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
parser.add_argument("-xn", "--xname",type=str, default='x',\
  help="Specify the x coordinate. e.g. xu or x. Default='x' ")
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
filename  = args.filename
normalize = args.normalize
Nbins     = args.nbins
mode     = args.mode
cl        = abs(args.coarse)

#==========================================================# 
# Create a dict that is passed into the function read3dDataFromNetCDF
nameDict = dict()
nameDict['xname'] = args.xname
nameDict['yname'] = args.yname
nameDict['zname'] = args.zname
nameDict['varname'] = args.varname


dataDict = read3dDataFromNetCDF( filename , nameDict, cl )
v = dataDict['v']
x = dataDict['x']; y = dataDict['y']; z = dataDict['z']
time = dataDict['time']


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

ixyz1 = input(" (1) Enter starting indices: ix, iy, iz = ")
ixyz2 = input(" (2) Enter final indices:    ix, iy, iz = ")

if( len(ixyz1) != 3  or len(ixyz2) != 3 ):
  sys.exit(' Error! You must provide 3 values for ix, iy and iz. Exiting ...')
ixyz1 = sensibleIds( np.array( ixyz1 ), x, y, z )
ixyz2 = sensibleIds( np.array( ixyz2 ), x, y, z )

ixL = np.arange(ixyz1[0],ixyz2[0]+1)
iyL = np.arange(ixyz1[1],ixyz2[1]+1)
izL = np.arange(ixyz1[2],ixyz2[2]+1)
stride = max( (izL[-1]-izL[0])/5 , 2 )

fig = None

for i in ixL:
  for j in iyL:
    for k in izL[::stride]:
      vt = v[:,int(k),int(j),int(i)]
      vName = nameDict['varname']+'(z={} m)'.format(z[k])
      fig = spectraAnalysis(fig, vt, time, vName, Nbins, mode, normalize)

plt.legend(loc=0)
plt.show()
