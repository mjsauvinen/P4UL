#!/usr/bin/env python
import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
from analysisTools import sensibleIds, groundOffset, quadrantAnalysis
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
parser.add_argument("-fo", "--fileout", type=str, default="out.nc", \
  help="Name of the output NETCDF file. Default=out.nc")
parser.add_argument("-v", "--varnames",  type=str, nargs=2, default=['u','w'],\
  help="Name of the variables in NETCDF file. Default=[u, w]")
parser.add_argument("-nk", "--nkpoints",type=int, default=None,\
  help="Number of data points used within iz1 -> iz2 interval. Default='All' ")
parser.add_argument("-nx", "--npixels",type=int, default=40,\
  help="Number of pixels in the quadrant plot (should be even). Default=40 ")
parser.add_argument("-l", "--axisLim",type=float, default=4.,\
  help="Limit (mag of max and min) of the plot axes. Default=4.")
parser.add_argument("-hw", "--holewidth",type=float, default=0.,\
  help="Width of the 'hole' in the quadrant analysis. Default=0.")
parser.add_argument("-i1", "--ijk1",type=int, nargs=3,\
  help="Starting indices (ix, iy, iz) of the considered data. Required.")
parser.add_argument("-i2", "--ijk2",type=int, nargs=3,\
  help="Final indices (ix, iy, iz) of the considered data. Required.")
parser.add_argument("-us", "--ustar",type=float, default=None,\
  help="Normalize by given friction velocity (u* or ustar) value.")
parser.add_argument("-of", "--outputToFile", type=str, default=None, \
  help="Name of the file to output analysis results. Default=None")
parser.add_argument("-w", "--weighted", action="store_true", default=False,\
  help="Plot weighted joint PDF (i.e. covariance integrand).")
parser.add_argument("-p", "--printOn", action="store_true", default=False,\
  help="Print the numpy array data.")
args = parser.parse_args()    
#==========================================================# 
# Rename ...
filename  = args.filename
fileout   = args.fileout
varnames  = args.varnames
ijk1      = args.ijk1
ijk2      = args.ijk2
axisLim   = args.axisLim
nkpoints  = args.nkpoints
ustar     = args.ustar
holewidth = args.holewidth
weighted  = args.weighted
npixels   = args.npixels
ofile     = args.outputToFile
printOn   = args.printOn
#==========================================================# 
'''
Establish two boolean variables which indicate whether the created variable is an
independent or dependent variable in function createNetcdfVariable().
'''
parameter = True;  variable  = False


# First fluctuation component
cl = 1
ncDict = read3dDataFromNetCDF( filename , varnames[0], cl )
v1 = ncDict['v']   # 'v' is a generic name for a variable in ncDict

# Second fluctuation component
ncDict = read3dDataFromNetCDF( filename , varnames[1], cl )
v2 = ncDict['v']

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
ijk1 = sensibleIds( np.array( ijk1 ), x, y, z )
ijk2 = sensibleIds( np.array( ijk2 ), x, y, z )
print(' Check (1): i, j, k = {}'.format(ijk1))
print(' Check (2): i, j, k = {}'.format(ijk2))

# = = = = = = = = = = = = =
# Mean and variance
v1mean = np.mean(v1, axis=(0)); v2mean = np.mean(v2, axis=(0))

# Extract fluctuating part and normalize by variance
# Reuse the v1 and  v2 variables to store values
v1 -= v1mean; v2 -= v2mean  

if( ustar is not None ):
  # Normalize by friction velocity
  v1 /= ustar; v2 /= ustar  
else:  
  # Normalize by variance
  v1var  = np.var(v1, axis=(0)) ; v2var  = np.var(v2, axis=(0))
  v1 /= v1var ; v2 /= v2var

# Assemble the quadrant analysis dict 
qaDict = dict()
qaDict['nkpoints']  = nkpoints
qaDict['npixels']   = npixels
qaDict['axisLim']   = axisLim
qaDict['holewidth'] = holewidth
qaDict['weighted']   = weighted  # covariance integrand

# = = = = = = = = = = = = = = =
xa = np.linspace(-axisLim,axisLim,npixels+1)
ya = xa.copy() 
dx = (2.*axisLim)/(npixels)
Xa,Ya = np.meshgrid(xa,ya)
ka    = np.arange(ijk1[2],ijk2[2]+1)
dims  = list(np.shape(ka))  
dims.extend(list(np.shape( Xa )))
Qa  = np.zeros( dims, float )
print(' Qa.shape={}'.format(np.shape(Qa)))

ix1 = ijk1.copy(); ix2 = ijk2.copy()
for kt in xrange(len(ka)):
  ix1[2] = ka[kt]
  ix2[2] = ka[kt]+1
  qaDict['ijk1'] = ix1; qaDict['ijk2'] = ix2
  
  Q, X, Y, resDict = quadrantAnalysis( v1, v2, qaDict )

  # Extract the results 
  nQ = resDict['nQ']  # Number of quadrant hits (nQ[0] := Ntotal)
  SQ = resDict['SQ']  # Quadrant contributions (e.g. Reynolds stress)

  cn = 100./nQ[0]
  print(' Ejections (%) = {}, Sweeps (%) = {} '.format(cn*nQ[2],cn*nQ[4]))
  print(' Outward Interactions (%)  = {}, Inward Interactions (%) = {} '.format(cn*nQ[1],cn*nQ[3]))

  # Write/append the results to file 
  if( ofile is not None ):
    Sa = np.abs(SQ[0])
    Smag = np.sqrt( np.sum(SQ[1:]**2) )
    zm = z[ix1[2]]
    for i in xrange(1,5):
      fwo = openIOFile('{}_Q{}.dat'.format(ofile,i) , 'a')
      fwo.write("{}\t{}\n".format(zm, SQ[i]))
      fwo.close()
  
  print(' Qa[kt,:,:].shape = {}, Q.shape = {}'.format(np.shape(Qa[kt,:,:]), np.shape(Q)))
  Qa[kt,:,:] = Q.copy()



# = = output file = = = = =
# Create a NETCDF output dataset (dso) for writing out the data.
dso = netcdfOutputDataset( fileout )
xv = createNetcdfVariable( dso, xa  , 'x'   , len(xa)   , 'm', 'f4', ('x',)   , parameter )
yv = createNetcdfVariable( dso, ya  , 'y'   , len(ya)   , 'm', 'f4', ('y',)   , parameter )
zv = createNetcdfVariable( dso, ka  , 'z'   , len(ka)   , 'm', 'f4', ('z',)   , parameter )
Qv = createNetcdfVariable( dso, Qa, 'Q', dims[0], 'm-2', 'f4',('z','y','x',) , variable )

# - - - - Done , finalize the output - - - - - - - - - -
netcdfWriteAndClose( dso )