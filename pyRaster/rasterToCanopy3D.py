#!/usr/bin/env python3
import sys
import argparse
import numpy as np
from utilities import writeLog
from mapTools import *
from netcdfTools import *



'''
Description:
Generates a 3D LAD array from tree height data based on given tree profile.
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='rasterToCanopy3D.py',
                                 description='''Writes PLANT_CANOPY_DATA_3D for PALM from input raster data.''')
parser.add_argument("-f","--filename", type=str, 
                    help="Name of the input raster data file.")
parser.add_argument("-fo", "--fileout", type=str, 
                    help="Name of the output 3D data file. ")
parser.add_argument("-dz", "--dz", type=float, default=None, 
  help="Resolution of the z axis. Defaults to resolution of N axis.")
parser.add_argument("-m", "--method", type=str, default='const', choices=['prof', 'const','file'],
                    help="Method for LAD distribution. Opt 'prof': Alfa-beta profile in "
                    "z-direction. Opt 'const': Constant value per cell. Opt. 'file "
                    "<LAD_profile_file>' read vertical LAD distribution from file "
                    "<LAD_profile_file>.")
parser.add_argument('LAD_file', nargs='?', help="File with the LAD profile, required if mode is 'file'")
parser.add_argument("-a", "--alpha", type=float, default=None, 
                    help="Method 'prof': Dimensionless coefficient required for constructing "
                    "the leaf area density (LAD) profile, using beta probability density "
                    "function (Markkanen et al., 2003, BLM 106, 437-459).")
parser.add_argument("-b", "--beta", type=float, default=None, 
                    help="Method 'prof': Dimensionless coefficient required for constructing "
                    "the leaf area density 2 (LAD) profile, using beta probability density function.")
parser.add_argument("-l", "--lai", type=float, default=6.,
                    help="Reference leaf area index (LAI) value. Method 'prof': Reference "
                    "LAI is the vertical integral over the reference tree's LAD profile. Method "
                    "'const': Reference LAI, which will be used to evaluate constant <LAD>_z = "
                    "LAI_ref/(zref[1]-zref[0]), where zref[0] and zref[1] refer to the values "
                    "given also as input. Default=6.")
parser.add_argument("-zr", "--zref", type=float, nargs=2, metavar=('ZREF[0]','ZREF[1]'),
                    default=[4.,20.], help=" The starting height of the foliage and the maximum "
                    "height of the reference tree whose LAI is given as input. Default=[4,20].")
parser.add_argument("-am", "--asmask", action="store_true", default=False, 
                    help="Output a netCDF4 3D boolean array mask for visualization purposes "
                    "instead of a npz file containing LAD values.")
parser.add_argument("-t", "--threshold", type=float, default=0.0, help="Threshold LAD value to "
                    "be used when generating a 3D mask. Grid points with a LAD value over the "
                    "threshold will be set to 1 while the rest is set to 0. Effective only if "
                    "--asmaks is set.")
args = parser.parse_args()
writeLog(parser, args)

#==========================================================#
# Renaming ... that's all
filename = args.filename
method   = args.method
alpha    = args.alpha
beta     = args.beta
laiRef   = args.lai
dz       = args.dz
fileout  = args.fileout
zref     = args.zref

constantLAD = ( method == 'const' )
profileLAD  = ( method == 'prof'  )

if( profileLAD ):
  if( (alpha is None) or (beta is None) ):
    sys.exit(' Error: alpha and/or beta is None. Exiting ...')

if method == 'file' and not args.LAD_file:
  sys.exit('With file method for LAD distribution, file with LAD profile is required.')

Rdict = readNumpyZTile( filename )
R = Rdict['R']
nPx = np.shape(R)
dPx = np.abs( Rdict['dPx'] )  #  <--- The deltas should be all positive here.

# Calculate the shape of the new 3D array, use largest in-canopy value
nPx3D = nPx; dPx3D = dPx
if ( dz ):
  dPx3D = np.append(dPx, dz)
else:
  dPx3D = np.append(dPx, dPx[0])

nPx3D = np.append(nPx3D, int(np.floor(np.amax(R)/dPx3D[2])+1))

# Fill the 3D canopy array
canopy = np.zeros([nPx3D[1], nPx3D[0], nPx3D[2]])
nPc=[nPx3D[1], nPx3D[0], nPx3D[2]]
dPc=[dPx3D[1], dPx3D[0], dPx3D[2]]
print(" 3D grid array dimensions [x,y,z]: {}, {}, {}".format(*nPc))
print(" Generating vertical distributions of leaf area densities ...")


# Compute <LAD>_z and starting index of the foliage using reference values
lad_const = laiRef/(zref[1]-zref[0])
k1        = int( np.round(zref[0]/float(dPc[2])) )  # starting k index

print(' Rry shape = {} '.format(R.shape))

# Calculate leaf area density profiles for each horizontal grid tile and fill array vertically
for j in range(nPc[1]):
  for i in range(nPc[0]):
    Zi = R[j,i] # Zi := canopy height at [j,i]
    # Check if there is canopy at all in the vertical column
    if (Zi <= zref[0]):
      continue

    # Number of layers
    dZ   = Zi - zref[0]
    nind = int(np.floor( dZ/float(dPc[2]) )) + 1

    if( profileLAD and (nind > 3) ):
      # Calculate LAD profile
      lai = laiRef * (Zi-zref[0])/(zref[1]-zref[0])  # Linear scaling of reference LAI.
      lad = canopyBetaFunction(dZ,dPc, alpha, beta, lai)
      k2C = k1 + min( nind, len(lad) )
      k2L = k2C - k1
      canopy[i,j,k1:k2C] = lad[0:k2L] #  Grid point at canopy top level gets value 0
    else:
      k2 = int(np.ceil(Zi/dPc[2]))+1
      k2 = min( k2, nPc[2] )
      if profileLAD:
        canopy[i,j,k1:k2] = lad_const
      else if method=='file':
        ladf = np.loadtxt(args.LAD_file)
        if ladf.shape[0] != ladf.size:
          sys.exit('Only single column LAD profile files are currently supported. Sorry!')
        # Scale given profile based on zref values and interpolate.
        rz = np.linspace(zref[0],zref[1],ladf.size)
        Pz = np.arange(k1,k2)*dPx3D[2]
        canopy[i,j,k1:k2] = np.interp(Pz, rz, ladf)
      

print(" ... done.\n")

# Write output data file
print(" Writing output file...")
if (args.asmask):
  # Use threshold value to create a mask raster.
  canopymask=np.zeros(np.shape(canopy))
  canopymask[np.where(canopy>args.threshold)]=1

  # Save as netCDF4
  # Maybe someday PALM can read canopy data from NetCDF4
  dso = netcdfOutputDataset(args.fileout)
  # Create dimensions
  xv = createCoordinateAxis(dso, nPc, dPc, 1, 'x', 'f4', 'm', parameter=True)
  yv = createCoordinateAxis(dso, nPc, dPc, 0, 'y', 'f4', 'm', parameter=True)
  zv = createCoordinateAxis(dso, nPc, dPc, 2, 'z', 'f4', 'm', parameter=True)
  # Due to a bug Paraview cannot read x,y,z correctly so rolling to z,y,x
  canopymask=np.rollaxis(canopymask,2)
  canopymask=np.swapaxes(canopymask,1,2)
  masknc = createNetcdfVariable(dso, canopymask, "canopy_0", 0, 'm', 'i4', ('z', 'y', 'x'), parameter=False)
  netcdfWriteAndClose(dso)

else:
  # Save as Numpy Z file.
  canopy = np.rollaxis( canopy, 1, 0 ) # S[i,j,k] -> S[j,i,k] to maintain compatibility with rasters
  Rdict.pop('R', None)
  Rdict.pop('Rdims', None )
  Rdict['S']     = canopy
  Rdict['dPx']   = dPx3D  # This is correct. It's format is [dN,dE,dZ]
  Rdict['Sdims'] = np.array( canopy.shape )
  Rdict['GlobOrig'] = np.append( Rdict['GlobOrig'] , 0. )
  if('GlobOrigBL' in Rdict.keys() ):
    Rdict['GlobOrigBL'] = np.append( Rdict['GlobOrigBL'] , 0. )
  saveTileAsNumpyZ( fileout, Rdict )

