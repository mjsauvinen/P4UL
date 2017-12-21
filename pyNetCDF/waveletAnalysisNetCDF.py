#!/usr/bin/env python
import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
from plotTools import addContourf
from analysisTools import sensibleIds, groundOffset, discreteWaveletAnalysis
from netcdfTools import read3dDataFromNetCDF, netcdfOutputDataset, \
  createNetcdfVariable, netcdfWriteAndClose
from utilities import filesFromList
from txtTools import openIOFile
''' 
Description: A script to perform wavelet analysis on velocity data stored in a NETCDF file.

Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
        
        Simone Boi
        University of Helsinki
'''
#==========================================================#
sepStr = ' # = # = # = # = # = # = # = # = '
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filename", type=str,\
  help="Name of the input NETCDF file.")
parser.add_argument("-wt", "--wavelet", type=str,\
  help="Name of the wavelet. See documentation of PyWavelet package.")
parser.add_argument("-nl", "--nlevel", type=int, default=4,\
  help="Number of levels in wavelet analysis. Default=4")
parser.add_argument("-fo", "--fileout", type=str, default="out.nc", \
  help="Name of the output NETCDF file. Default=out.nc")
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
parser.add_argument("-s", "--save", type=str, default=None, \
  help="Identifier (name) for the saved figures. Default=None")
parser.add_argument("-of", "--outputToFile", type=str, default=None, \
  help="Name of the file to output analysis results. Default=None")
parser.add_argument("-p", "--printOn", action="store_true", default=False,\
  help="Print the numpy array data.")
args = parser.parse_args()    
#==========================================================# 
# Rename ...
filename  = args.filename
fileout   = args.fileout
wavelet   = args.wavelet
nlevel    = args.nlevel
varname   = args.varname
kIds      = args.kIndices
nkpoints  = args.nkpoints
ofile     = args.outputToFile
saveFig   = args.save
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

# Plot coord. information. This aids the user in the beginning.
infoStr = '''
  Coord. range:
  min(x)={0} ... max(x)={1}, nx = {2}
  min(y)={3} ... max(y)={4}, ny = {5}
  min(z)={6} ... max(z)={7}, nz = {8}
'''.format(\
    np.min(x), np.max(x), len(x),\
    np.min(y), np.max(y), len(y),\
    np.min(z), np.max(z), len(z) )
#print(infoStr)

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
wlDict['nlevel']  = nlevel
# = = = = = = = = = = = = = = =
kList = np.arange(kIds[0],kIds[1]+1)
for kt in kList:
  vt = v[:,kt, ijk1[1], ijk1[0]]
  
  values, labels = discreteWaveletAnalysis( vt, wlDict )


  fig = plt.figure()
  fig.subplots_adjust(hspace=0.2, bottom=.03, left=.07, right=.97, top=.92)
  plt.subplot(2, 1, 1)
  plt.title(" signal")
  plt.plot(time, vt, 'b')
  plt.xlim(time[0], time[-1])

  ax = plt.subplot(2, 1, 2)
  plt.title("Wavelet packet coefficients at level {}".format(nlevel))
  plt.imshow(values, interpolation='nearest', aspect="auto",\
    origin="lower", extent=[0, 1, 0, len(values)])
  plt.yticks(np.arange(0.5, len(labels) + 0.5), labels)
  #plt.setp(ax.get_xticklabels(), visible=False)

  #plt.figure(2)
  #plt.specgram(data, NFFT=64, noverlap=32, cmap=cmap)
  #plt.imshow(values, origin='upper', extent=[-1,1,-1,1],
  # interpolation='nearest')
  
  if( saveFig ):
    plt.savefig( saveFig+'_{}m.jpg'.format(int(z[kt])), format='jpg', dpi=300)
    
  fig = None


'''
# = = output file = = = = =
# Create a NETCDF output dataset (dso) for writing out the data.
dso = netcdfOutputDataset( fileout )
xv = createNetcdfVariable( dso, xa  , 'x'   , len(xa)   , 'm', 'f4', ('x',)   , parameter )
yv = createNetcdfVariable( dso, ya  , 'y'   , len(ya)   , 'm', 'f4', ('y',)   , parameter )
zv = createNetcdfVariable( dso, ka  , 'z'   , len(ka)   , 'm', 'f4', ('z',)   , parameter )
Qv = createNetcdfVariable( dso, Qa, 'Q', dims[0], 'm-2', 'f4',('z','y','x',) , variable )

# - - - - Done , finalize the output - - - - - - - - - -
netcdfWriteAndClose( dso )
'''