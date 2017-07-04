#!/usr/bin/env python
import sys
import argparse
import numpy as np
from utilities import writeLog
from mapTools import *
from scipy.special import gamma

'''
Description:
Writes PLANT_CANOPY_DATA_3D for PALM from input raster data.
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='rasterToCanopy3D.py', description='''Writes PLANT_CANOPY_DATA_3D for PALM from input raster data.''')
parser.add_argument("rfile", type=str, nargs='?', default=None, help="Name of the input raster data file.")
parser.add_argument("-fo", "--fileout", type=str, help="Name of the output 3D data file.")
parser.add_argument("-dz", "--dZ", type=float, help="Resolution of z axis. Defaults to resolution of N axis.")
parser.add_argument("-a", "--alpha", type=float, help="Dimensionless coefficient required for constructing the leaf area density (LAD) profile, using beta probability density function (Markkanen et al., 2003, BLM 106, 437-459).")
parser.add_argument("-b", "--beta", type=float, help="Dimensionless coefficient required for constructing the leaf area density (LAD) profile, using beta probability density function.")
parser.add_argument("-l", "--lai", type=float, help="Leaf area index LAI value. LAI is the vertical integral over the LAD profile.")
args = parser.parse_args()
writeLog(parser, args)

#==========================================================#

Rdict = readNumpyZTile(args.rfile)
R = Rdict['R']
nPx = np.shape(R)
dPx = Rdict['dPx']
Rdict = None

# Calculate the shape of the new 3D array, use largest in-canopy value
nPx3D = nPx; RdPx3D = dPx
if (args.dZ):
  dPx3D = np.append(dPx, args.dZ)
else:
  dPx3D = np.append(dPx, dPx[0])

nPx3D = np.append(nPx3D, int(np.floor(np.amax(R)/dPx3D[2])+1))

# Fill the 3D canopy array
canopy = np.zeros([nPx3D[1], nPx3D[0], nPx3D[2]])
print(" 3D grid array dimensions [x,y,z]: {}, {}, {}".format(*np.shape(canopy)))
print(" Calculating vertical profiles of leaf area densities...")

# Calclulate the value of the integral in the denominator only once
d_integral = (gamma(args.alpha)*gamma(args.beta))/gamma(alpha+beta)

# Calculate leaf area density profiles for each horizontal grid tile and fill array vertically
for x in xrange(nPx3D[1]):
  for y in xrange(nPx3D[0]):
    # Reverse the y-axis because of the top-left origo in raster
    canopyHeight = R[-y - 1, x]
    nind = int(np.floor(canopyHeight / float(dPx3D[2])))+1 # Number of layers
    # Calculate LAD profile
    lad = canopyBetaFunction(nind,nPx3D[2],dPx3D[2],args.alpha,args.beta,args.lai,d_integral)
    canopy[x,y,0:nind] = lad[:]
print(" ...done.")

# Write output data file
print(" Writing output output file...")
fx = open(args.fileout, 'w')
# The first row of the file is number of vertical canopy layers
fx.write(str(nPx3D[2])+"\n")
# Loop through the verical canopy columns
for y in xrange(nPx3D[0]):
  for x in xrange(nPx3D[1]):
    # Convert everything into strings and write
    # datatype x y col(:)
    lineStr = str(1)+","+str(x)+","+str(y)+","+",".join(map(str,canopy[x,y,:]))+"\n"
    fx.write(lineStr)
fx.close()
print(" ...{} saved successfully.".format(args.fileout))
