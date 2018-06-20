#!/usr/bin/env python
import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
from analysisTools import sensibleIds, groundOffset, continuousWaveletAnalysis
from netcdfTools import read3dDataFromNetCDF, netcdfOutputDataset, \
  createNetcdfVariable, netcdfWriteAndClose
from utilities import filesFromList
from txtTools import openIOFile
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
sepStr = ' # = # = # = # = # = # = # = # = '
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filename", type=str,\
  help="Name of the input NETCDF file.")
parser.add_argument("-fo", "--fileout", type=str, default="cwtout.nc", \
  help="Name of the output NETCDF file. Default=cwtout.nc")
parser.add_argument("-wt", "--wavelet", type=str, default='morl',\
  help="Name of the wavelet. See documentation of PyWavelet package. Default='morl' ")
parser.add_argument("-nf", "--nfreqs", type=int, default=256,\
  help="Number of frequency levels in continuous wavelet transform. Default=256")
parser.add_argument("-v", "--varname",  type=str, default='u',\
  help="Name of the variable in NETCDF file. Default='u' ")
parser.add_argument("-xn", "--xname",type=str,  default='xu',\
  help="Specify the x coordinate. e.g. xu or x. Default='xu' ")
parser.add_argument("-yn", "--yname",type=str, default='y',\
  help="Specify the y coordinate. e.g. yv or y. Default='y' ")
parser.add_argument("-zn", "--zname",type=str, default='zu_3d',\
  help="Specify the z coordinate. e.g. z, zu_3d or zw_3d. Default='zu_3d' ")
parser.add_argument("-nk", "--nkpoints",type=int, default=None,\
  help="Number of data points used within iz1 -> iz2 interval. Default='All' ")
parser.add_argument("-k", "--kIndices",type=int, nargs=2,\
  help="Starting and final index (k_start, k_final) of the considered data. Required.")
parser.add_argument("--linearFreq", action="store_true",\
  help="Use linear frequency distribution.")
parser.add_argument("-p", "--printOn", action="store_true", default=False,\
  help="Print the numpy array data.")
args = parser.parse_args()    
#==========================================================# 
# Rename ...
filename  = args.filename
fileout   = args.fileout
wavelet   = args.wavelet
nfreqs    = args.nfreqs
varname   = args.varname
kIds      = args.kIndices
nkpoints  = args.nkpoints
linearFreq= args.linearFreq
printOn   = args.printOn
#==========================================================# 
'''
Establish two boolean variables which indicate whether the created variable is an
independent or dependent variable in function createNetcdfVariable().
'''
parameter = True;  variable  = False

# Create a dict that is passed into the function read3dDataFromNetCDF
nameDict = dict()
nameDict['xname'] = args.xname
nameDict['yname'] = args.yname
nameDict['zname'] = args.zname

# First fluctuation component
nameDict['varname'] = varname[0]
cl = 1
ncDict = read3dDataFromNetCDF( filename , nameDict, cl )
v = ncDict['v']   # 'v' is a generic name for a variable in ncDict

# Spatial coords and time
x = ncDict['x']; y = ncDict['y']; z = ncDict['z']
time = ncDict['time']

# Now check whether the given indices make sense 
# Here we set i = 0 and j = 0.
ijk1 = sensibleIds( np.array([0,0,kIds[0]]), x, y, z )
ijk2 = sensibleIds( np.array([0,0,kIds[1]]), x, y, z )
print(' Check (1): i, j, k = {}'.format(ijk1))
print(' Check (2): i, j, k = {}'.format(ijk2))

# = = = = = = = = = = = = =
# Mean and variance
vmean = np.mean(v, axis=(0))

# Extract fluctuating part and normalize by variance
# Reuse the v variable
v -= vmean

# Assemble the quadrant analysis dict 
wlDict = dict()
wlDict['wavelet'] = wavelet
wlDict['nfreqs']  = nfreqs
wlDict['dt']      = np.round(np.mean((time[1:]-time[0:-1])) , decimals=2)
wlDict['linearFreq'] = linearFreq 

# = = = = = = = = = = = = = = =
# Allocate the 3D array for data storage
xa = np.arange(len(time))
ya = 1./np.arange(1, nfreqs)
ya = 1./ya
Xa,Ya = np.meshgrid(xa,ya)
ka    = np.arange(ijk1[2],ijk2[2]+1)
dims  = list(np.shape(ka))  
dims.extend(list(np.shape( Xa )))
Cfs   = np.zeros( dims, float )
print(' Cfs.shape={}'.format(np.shape(Cfs)))

# = = = = = = = = = = = = = = =
scales = ya  # Rename, that's all
kList  = ka
k      = 0
for kt in kList:
  vt = v[:,kt, ijk1[1], ijk1[0]]
  
  C, freqs = continuousWaveletAnalysis( vt, wlDict )
  
  print(' Cfs[kt,:,:].shape = {}, Q.shape = {}'.format(np.shape(Cfs[k,:,:]), np.shape(C)))
  Cfs[k,:,:] = C.copy()
  k += 1


# = = output file = = = = =
# Create a NETCDF output dataset (dso) for writing out the data.
dso = netcdfOutputDataset( fileout )
xv = createNetcdfVariable( dso, xa  , 'x'   , len(xa)   , 'm', 'f4', ('x',)   , parameter )
yv = createNetcdfVariable( dso, ya  , 'y'   , len(ya)   , 'm', 'f4', ('y',)   , parameter )
zv = createNetcdfVariable( dso, ka  , 'z'   , len(ka)   , 'm', 'f4', ('z',)   , parameter )
Qv = createNetcdfVariable( dso, Cfs, 'Cfs', dims[0], '-', 'f4',('z','y','x',) , variable )

# - - - - Done , finalize the output - - - - - - - - - -
netcdfWriteAndClose( dso )