#!/usr/bin/env python
import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
from plotTools import addContourf
from analysisTools import sensibleIds, groundOffset, quadrantAnalysis
from netcdfTools import read3dDataFromNetCDF
from utilities import filesFromList
from txtTools import openIOFile
''' 
Description: A script to perform quadrant analysis on velocity data stored in a NETCDF file.
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
parser.add_argument("-i1", "--ijk1",type=int, nargs=3, default=[0,0,0],\
  help="Starting indices (ix, iy, iz) of the considered data. Default=[0,0,0].")
parser.add_argument("-i2", "--ijk2",type=int, nargs=3, default=[0,0,1],\
  help="Final indices (ix, iy, iz) of the considered data. Default=[0,0,1].")
parser.add_argument("-us", "--ustar",type=float, default=None,\
  help="Normalize by friction velocity, u_* or ustar.")
parser.add_argument("-s", "--save", type=str, default=None, \
  help="Name of the saved figure. Default=None")
parser.add_argument("-of", "--outputToFile", type=str, default=None, \
  help="Name of the file to output analysis results. Default=None")
parser.add_argument("-w", "--weighted", action="store_true", default=False,\
  help="Plot weighted joint PDF (i.e. covariance integrand).")
parser.add_argument("-p", "--printOn", action="store_true", default=False,\
  help="Print the numpy array data.")
parser.add_argument("-c", "--coarse", type=int, default=1,\
  help="Coarsening level. Int > 1.")
args = parser.parse_args()    
#==========================================================# 
# Rename ...
filename  = args.filename
cl        = abs(args.coarse)
varnames  = args.varnames
saveFig   = args.save
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
first = True
fig   = None


# First fluctuation component
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
print(infoStr)

# Now check whether the given indices make sense 
ijk1 = sensibleIds( np.array( ijk1 ), x, y, z )
ijk2 = sensibleIds( np.array( ijk2 ), x, y, z )
print(' Check (1): i, j, k = {}'.format(ijk1))
print(' Check (2): i, j, k = {}'.format(ijk2))

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
qaDict['ijk1'] = ijk1
qaDict['ijk2'] = ijk2
qaDict['nkpoints']  = nkpoints
qaDict['npixels']   = npixels
qaDict['axisLim']   = axisLim
qaDict['holewidth'] = holewidth
qaDict['weighted']   = weighted  # covariance integrand

Q, X, Y, resDict = quadrantAnalysis( v1, v2, qaDict )

# Extract the results 
nQ = resDict['nQ']  # Number of quadrant hits (nQ[0] := Ntotal)
SQ = resDict['SQ']  # Quadrant contributions (e.g. Reynolds stress)
#klims         = resDict['klims']

# === Plot quadrant analysis output === #
cDict = dict() 
cDict['cmap'] = plt.cm.gist_yarg  # Include the colormap info within a dict
cDict['N'] = 12                   # Number of levels in contour plot
cDict['title'] = "Quadrant Analysis\n{}:  z={}-{} m".format(filename, z[ijk1[2]],z[ijk2[2]])
cDict['label'] = "JPDF"
CO = addContourf( X, Y, Q, cDict )
CO.ax.spines['left'].set_position('zero')
CO.ax.spines['bottom'].set_position('zero')
#CO.ax.set_ylabel(r"$w'/\sigma_w$")
#CO.ax.set_xlabel(r"$u'/sigma_u$")
#plt.clabel(CO, CO.levels[:-2:3], inline=False, fontsize=10)

cn = 100./nQ[0]
print(' Ejections (%) = {}, Sweeps (%) = {} '.format(cn*nQ[2],cn*nQ[4]))
print(' Outward Interactions (%)  = {}, Inward Interactions (%) = {} '.format(cn*nQ[1],cn*nQ[3]))
#plt.legend(loc=0)

if( saveFig ):
  plt.savefig( saveFig, format='jpg', dpi=300)

if( printOn ):
  plt.show()

CO = None

