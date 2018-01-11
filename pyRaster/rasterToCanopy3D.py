#!/usr/bin/env python
import sys
import argparse
import numpy as np
from utilities import writeLog
from mapTools import *
from netcdfTools import *



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
parser.add_argument("-am", "--asmask", type=str, default=None, help="Output a 3D array mask instead of a text file formatted for PALM. The argument is a file type, use nc for NetCDF4 output and npz for Numpy Z array.")
parser.add_argument("-t", "--threshold", type=float, default=0.0, help="Threshold LAD value to be used when generating a 3D mask. Grid points with a LAD value over the threshold will be set to 1 while the rest is set to 0. Effective only if --asmaks is set.")
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

# Calculate the shape of the new 3D array, use largest in-canopy value
nPx3D = nPx; dPx3D = dPx
if ( dZ ):
  dPx3D = np.append(dPx, args.dZ)
else:
  dPx3D = np.append(dPx, dPx[0])

nPx3D = np.append(nPx3D, int(np.floor(np.amax(R)/dPx3D[2])+1))

# Fill the 3D canopy array
canopy = np.zeros([nPx3D[1], nPx3D[0], nPx3D[2]])
nPx3D=[nPx3D[1], nPx3D[0], nPx3D[2]]
print(" 3D grid array dimensions [x,y,z]: {}, {}, {}".format(*nPx3D))
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
if (args.asmask):
  # Use threshold value to create a mask raster.
  canopymask=np.zeros(np.shape(canopy))
  canopymask[np.where(canopy>args.threshold)]=1
  if (args.asmask=="npz"):
    # Save as Numpy Z file.
    Rdict["R"]=canopymask
    saveTileAsNumpyZ( fileout, Rdict )
  elif (args.asmask=="nc"):
    # Save as netCDF4
    # Maybe someday PALM can read canopy data from NetCDF4
    dso = netcdfOutputDataset(args.fileout)
    # Create dimensions
    xv = createCoordinateAxis(dso, nPx3D, dPx3D, 1, 'x', 'f4', 'm', parameter=True)
    yv = createCoordinateAxis(dso, nPx3D, dPx3D, 0, 'y', 'f4', 'm', parameter=True)
    zv = createCoordinateAxis(dso, nPx3D, dPx3D, 2, 'z', 'f4', 'm', parameter=True)
    # Due to a bug Paraview cannot read x,y,z correctly so rolling to z,y,x
    canopymask=np.rollaxis(canopymask,2)
    canopymask=np.swapaxes(canopymask,1,2)
    masknc = createNetcdfVariable(dso, canopymask, "canopy_0", 0, 'm', 'i4', ('z', 'y', 'x'), parameter=False)
    netcdfWriteAndClose(dso)

else:
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
