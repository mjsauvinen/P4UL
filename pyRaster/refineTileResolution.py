#!/usr/bin/env python3
import argparse
import numpy as np
from mapTools import *
from utilities import writeLog
from plotTools import addImagePlot
import matplotlib.pyplot as plt
import sys
'''
Description:


Author: Mikko Auvinen & Jukka-Pekka Keskinen
        mikko.auvinen@fmi.fi
        Finnish Meteorological Institute
'''

#==========================================================#
dscStr = '''Refine or coarsen a raster. When coarsening, the new value will be 
the mean of the corresponding cells in the fine raster unless specified otherwise. 
NB! Does not work properly if raster contains missing values or NaNs.
'''

parser = argparse.ArgumentParser(prog='refineTileResolution.py',\
  description=dscStr)
parser.add_argument("-f", "--filename",type=str, help="Name of the .npz data file.")
parser.add_argument("-fo", "--fileout",type=str,\
  help="Name of output Palm/npz topography file.")
parser.add_argument("-N","--refn", type=float,\
  help="Refinement factor N in 2^N. Negative value coarsens.")
parser.add_argument("-m", "--mode", action="store_true", default=False,\
  help="When coarsening, select the most common value (mode) from the finer grid.")
parser.add_argument("-i", "--integer", action="store_true", default=False,\
  help="Output an integer array. Default is a float array.")
parser.add_argument("-p", "--printOn", action="store_true", default=False, \
  help="Print the resulting raster data.")
parser.add_argument("-pp", "--printOnly", action="store_true", default=False, \
  help="Only print the resulting data. Don't save.")
args = parser.parse_args()
writeLog( parser, args, args.printOnly )
#==========================================================#

# Renaming ... nothing else
filename  = args.filename
N         = args.refn
printOn   = args.printOn
printOnly = args.printOnly
fileout   = args.fileout
mode      = args.mode
ints      = args.integer

Rdict = readNumpyZTile( filename )
R1 = Rdict['R']

if( N < 0. ):
  # Take precaution if coarsening and the original raster has an odd numbered dimension.
  nN, nE = R1.shape
  eN = None; eE = None
  if( np.mod( nN, 2 ) != 0 ): eN = -1
  if( np.mod( nE, 2 ) != 0 ): eE = -1
  if( eN is not None or eE is not None ): R1 = R1[:eN,:eE]

R1dims = np.array(np.shape(R1))
print(' R1dims = {}'.format(R1dims))
R1Orig = Rdict['GlobOrig'].astype(float)
dPx1   = Rdict['dPx']
gridRot1 = Rdict['gridRot']

if ints:
  Rtype = int
else:
  Rtype = float

# Resolution ratio (rr).
rr = 2**N

# Create the index arrays. The dims are always according to the larger one.
if( N > 0. ):
  dr1 = rr; fr2 = 1        # dr1 > 1
  R2dims = np.round( dr1 * R1dims ).astype(int)   # Refinement, R2dims > R1dims
  s2 = 1.                   # Scale factor. If we refine, the R1 values always fill a new zero cell, see below.
  n1,e1 = np.ogrid[ 0:R2dims[0] , 0:R2dims[1] ]  # northing, easting 
  n2,e2 = np.ogrid[ 0:R2dims[0] , 0:R2dims[1] ]  # northing, easting
  
else:
  dr1 = 1; fr2 = rr    # fr2 < 1
  R2dims  = np.round(rr * R1dims).astype(int); print(' Coarser dims = {}'.format(R2dims))
  s2 = (2**(2*N))   # Scale factor. If we coarsen, the R1 values are appended. Same value to 2^2n cells.
  n1 = np.arange(R1dims[0]); e1 = np.arange(R1dims[1])
  n2 = np.arange(R1dims[0]); e2 = np.arange(R1dims[1])

# Modify the integer list for refining/coarsening
n1= (n1/dr1).astype(int); e1=(e1/dr1).astype(int)

#n2 = np.round(n2*fr2, decimals=2);  e2 = np.round(e2*fr2, decimals=2)
n2 = np.floor((n2+0.5)*fr2).astype(int);  e2 = np.floor((e2+0.5)*fr2).astype(int)
#print(' n2 = {} '.format(n2))
#print(' n1 = {} '.format(n1))

np.savetxt('n2.dat', n2, fmt='%g')
#np.savetxt('n1.dat', n1, fmt='%g')

if( N > 0 ):
  R2 = np.zeros( R2dims, Rtype ) # Create the output array.
  R2[n2, e2] += R1[n1,e1]
elif  np.isclose(np.around(1/s2),1/s2,0.001):
  n2 = np.minimum( n2 , R2dims[0]-1)
  e2 = np.minimum( e2 , R2dims[1]-1)
  if mode:
    R2=slowCoarsen(R1,R2dims,s2,n1,n2,e1,e2,Rtype)
  else:
    R2=fastCoarsen(R1,R2dims,s2,n1,n2,e1,e2,Rtype)
else:
  sys.exit("ERROR: Attempting to coarsen with an incompatible refinement factor. Exiting.")


#print(' TL:{} TR:{} BL:{} BR:{} '.format( R2[0,0], R2[0,-1], R2[-1,0], R2[-1,-1]))
#print(' TL+1:{} TR+1:{} BL+1:{} BR+1:{} '.format( R2[1,0], R2[1,-1], R2[-1,1], R2[-1,-2]))

# NOTE! The global origin is the coordinate of the top left cell center. 
# Therefore, it must be shifted by in accordance to the top left cc's new location.
R1 = None
Rdict['R'] = R2.astype(Rtype)

# Select the smaller delta
dPx2 = dPx1/rr
dPm  = np.minimum( np.abs(dPx1), np.abs(dPx2) )  # np.abs() just in case negative dPx values sneak in.
dPm[1] = -1.*np.abs( dPm[1] )  # Set dE negative for the offset operation

# Offset the new cell-center origo
# Refine:  +dN, -dE
# Coarsen: -dN, +dE
R2Orig = R1Orig.copy()
R2Orig += np.sign(N)*(2.**np.abs(N) - 1.) * (dPm/2.) # N<0 coarsens, N>0 refines

# Rotate the newly shifted global origin to the global coord. system using the R1Orig as pivot.
R2Orig = rotatePoint(R1Orig, R2Orig, gridRot1)

print(' Old Origin = {} vs.\n New Origin = {}'.format(R1Orig, R2Orig))

Rdict['GlobOrig'] = R2Orig
Rdict['dPx'] = dPx2

if( not args.printOnly ):
  saveTileAsNumpyZ( fileout, Rdict )

if( args.printOn or args.printOnly ):
  fig = plt.figure(num=1, figsize=(9.,9.))
  fig = addImagePlot( fig, R2 , fileout )

  plt.show()
