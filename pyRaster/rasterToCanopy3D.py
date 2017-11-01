#!/usr/bin/env python
import sys
import argparse
import numpy as np
from utilities import writeLog
from mapTools import *

'''
Description:
Writes PLANT_CANOPY_DATA_3D for PALM from input raster data.
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='rasterToCanopy3D.py', description='''Writes PLANT_CANOPY_DATA_3D for PALM from input raster data.''')
parser.add_argument("-f","--filename", type=str, help="Name of the input raster data file.")
parser.add_argument("-fo", "--fileout", type=str, help="Name of the output 3D data file.", \
  default='PLANT_CANOPY_DATA_3D')
parser.add_argument("-dz", "--dZ", type=float, help="Resolution of the z axis. Defaults to resolution of N axis.")
parser.add_argument("-a", "--alpha", type=float, help="Dimensionless coefficient required for constructing the leaf area density (LAD) profile, using beta probability density function (Markkanen et al., 2003, BLM 106, 437-459).")
parser.add_argument("-b", "--beta", type=float, help="Dimensionless coefficient required for constructing the leaf area density (LAD) profile, using beta probability density function.")
parser.add_argument("-l", "--lai", type=float, help="Leaf area index LAI value. LAI is the vertical integral over the LAD profile.")
args = parser.parse_args()
writeLog(parser, args)

#==========================================================#

filename = args.filename
alpha    = args.alpha
beta     = args.beta
lai      = args.lai
dZ       = args.dZ
fileout  = args.fileout

Rdict = readNumpyZTile( filename )
R = Rdict['R']
nPx = np.shape(R)
dPx = Rdict['dPx']
Rdict = None

# Calculate the shape of the new 3D array, use largest in-canopy value
nPx3D = nPx; RdPx3D = dPx
if ( dZ ):
  dPx3D = np.append(dPx, args.dZ)
else:
  dPx3D = np.append(dPx, dPx[0])

nPx3D = np.append(nPx3D, int(np.floor(np.amax(R)/dPx3D[2])+1))

# Fill the 3D canopy array
canopy = np.zeros([nPx3D[1], nPx3D[0], nPx3D[2]])
print(" 3D grid array dimensions [x,y,z]: {}, {}, {}".format(*np.shape(canopy)))
print(" Calculating vertical profiles of leaf area densities...")

# Calculate leaf area density profiles for each horizontal grid tile and fill array vertically
for iy in xrange(nPx3D[0]):
  for ix in xrange(nPx3D[1]):
    # Reverse the y-axis because of the top-left origo in raster
    canopyHeight = R[-iy - 1, ix]
    # Check if there is canopy at all in the vertical column
    if (canopyHeight == 0.0):
      continue
    nind = int(np.floor(canopyHeight / float(dPx3D[2])))+1 # Number of layers
    if( nind > 3 ):
      # Calculate LAD profile
      lad = canopyBetaFunction(canopyHeight,dPx3D, alpha, beta, lai)
      # BUG: canopy[ix,iy,0:nind-1] = lad
      iend = min( nind, len(lad) )
      canopy[ix,iy,0:iend] = lad[0:iend] #  Grid point at canopy top level gets value 0
print(" ...done.")

# Write output data file
print(" Writing output file...")
fx = open( fileout, 'w')
# The first row of the file is number of vertical canopy layers
fx.write(str(nPx3D[2])+"\n")
# Loop through the verical canopy columns
for x in xrange(nPx3D[1]):
  for y in xrange(nPx3D[0]):
    if (np.all(canopy[x,y,:]==0)):
      # There is no need to write empty columns
      continue
    # Convert everything into strings and write
    # datatype x y col(:)
    lineStr = str(1)+","+str(x)+","+str(y)+","+",".join(map("{:.3g}".format,canopy[x,y,:]))+"\n"
    fx.write(lineStr)
fx.close()
print(" ...{} saved successfully.".format(fileout))
