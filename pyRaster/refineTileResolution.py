#!/usr/bin/env python3
import argparse
import numpy as np
from mapTools import *
from utilities import writeLog
from plotTools import addImagePlot
import matplotlib.pyplot as plt
'''
Description:


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi
        University of Helsinki &
        Finnish Meteorological Institute
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='refineTileResolution.py')
parser.add_argument("-f", "--filename",type=str, help="Name of the .npz data file.")
parser.add_argument("-fo", "--fileout",type=str, help="Name of output Palm/npz topography file.",\
  default="TOPOGRAPHY_MOD")
parser.add_argument("-N","--refn", help="Refinement factor N in 2^N. Negative value coarsens.", type=float)
parser.add_argument("-p", "--printOn", help="Print the resulting raster data.",\
  action="store_true", default=False)
parser.add_argument("-pp", "--printOnly", help="Only print the resulting data. Don't save.",\
  action="store_true", default=False)
args = parser.parse_args()
writeLog( parser, args, args.printOnly )
#==========================================================#

# Renaming ... nothing else
filename  = args.filename
N         = args.refn
printOn   = args.printOn
printOnly = args.printOnly
fileout   = args.fileout


Rdict = readNumpyZTile( filename )
R1 = Rdict['R']
R1dims = np.array(np.shape(R1))
R1Orig = Rdict['GlobOrig']
dPx1   = Rdict['dPx']
gridRot1 = Rdict['gridRot']

# Resolution ratio (rr).
rr = 2**N

# Create the index arrays. The dims are always according to the larger one.
if( N > 0. ):
  dr1 = rr; fr2 = 1        # dr1 > 1
  maxDims = dr1 * R1dims   # Refinement, R2dims > R1dims
  R2dims  = maxDims
  s2 = 1.                   # Scale factor. If we refine, the R1 values always fill a new zero cell, see below.
  n1,e1 = np.ogrid[ 0:maxDims[0] , 0:maxDims[1] ]  # northing, easting 
  n2,e2 = np.ogrid[ 0:maxDims[0] , 0:maxDims[1] ]  # northing, easting
else:
  dr1 = 1; fr2 = rr    # fr2 < 1
  maxDims = R1dims     # Coarsening, R2dims < R1dims
  R2dims  = np.round(rr * R1dims).astype(int); print(' Coarser dims = {}'.format(R2dims))
  s2 = (2**(2*N))   # Scale factor. If we coarsen, the R1 values are appended. Same value to 2^2n cells.
  n1 = np.arange(maxDims[0]); e1 = np.arange(maxDims[1])
  n2 = np.arange(maxDims[0]); e2 = np.arange(maxDims[1])

# Modify the integer list for refining/coarsening
n1/=dr1; e1/=dr1
n1 = n1.astype(int);  e1 = e1.astype(int)

n2 = n2*fr2        ;  e2 = e2*fr2
n2 = n2.astype(int);  e2 = e2.astype(int)

#print ' n2 = {} '.format(n2)
#print ' n1 = {} '.format(n1)



# Create the output array.
R2 = np.zeros( R2dims, float )

if( N > 0 ):
  R2[n2, e2] += R1[n1,e1]
else:
  n2 = np.minimum( n2 , R2dims[0]-1)
  e2 = np.minimum( e2 , R2dims[1]-1)
  for k in range(maxDims[0]):
    for l in range(maxDims[1]):
      R2[ n2[k], e2[l] ] +=  R1[ n1[k] ,e1[l] ]


# NOTE! The global origin is the coordinate of the top left cell center. 
# Therefore, it must be shifted by in accordance to the top left cc's new location.
R1 = None
R2 *= s2
Rdict['R'] = R2

dPx2 = dPx1/rr
dPm  = np.minimum( dPx1, dPx2 )
dPm[1] *= -1.  # When coarsening, the offset should be negative in N and positive in E
R2Orig = R1Orig.copy()
R2Orig += np.sign(N)*(2.*np.abs(N) - 1.) * (dPm/2.)  # N<0 coarsens, N>0 refines

# Rotate the newly shifted global origin to the global coord. system using the R1Orig as pivot.
R2Orig = rotatePoint(R1Orig, R2Orig, gridRot1)

print(' Old Origin = {} vs.\n New Origin = {}'.format(R1Orig, R2Orig))

Rdict['GlobOrig'] = R2Orig
Rdict['dPx'] = dPx2

if( not args.printOnly ):
  fx = open( fileout , 'w' )
  np.savetxt(fx,np.round(R2),fmt='%g')
  fx.close()
  saveTileAsNumpyZ( fileout, Rdict )

if( args.printOn or args.printOnly ):
  fig = plt.figure(num=1, figsize=(9.,9.))
  fig = addImagePlot( fig, R2 , fileout )

  plt.show()
