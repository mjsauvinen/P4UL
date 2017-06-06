#!/usr/bin/env python
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
parser.add_argument("-fo", "--fileOut",type=str, help="Name of output Palm/npz topography file.",\
  default="TOPOGRAPHY_MOD")
parser.add_argument("-n","--refn", help="Refinement factor n in 2^n (int). Negative value coarsens.", type=int)
parser.add_argument("-p", "--printOn", help="Print the resulting raster data.",\
  action="store_true", default=False) 
parser.add_argument("-pp", "--printOnly", help="Only print the resulting data. Don't save.",\
  action="store_true", default=False) 
args = parser.parse_args() 
writeLog( parser, args )
#==========================================================#

# Renaming ... nothing else
n         = args.refn
printOn   = args.printOn
printOnly = args.printOnly
fileOut   = args.fileOut


Rdict1 = readNumpyZTile( args.filename )
R1 = Rdict1['R']
R1dims = np.array(np.shape(R1))
R1Orig = Rdict1['LocalOrig']
dPx1 = Rdict1['dPx']
Rdict1 = None

# Resolution ratio (rr).
rr = 2**n

# Create the index arrays. The dims are always according to the larger one.
if( n > 0 ):
  dr1 = rr; fr2 = 1        # dr1 > 1
  maxDims = dr1 * R1dims   # Refinement, R2dims > R1dims
  R2dims  = maxDims
  s2 = 1.                   # Scale factor. If we refine, the R1 values always fill a new zero cell, see below.
  i1,j1 = np.ogrid[ 0:maxDims[0] , 0:maxDims[1] ]
  i2,j2 = np.ogrid[ 0:maxDims[0] , 0:maxDims[1] ]
else:
  dr1 = 1; fr2 = rr    # fr2 < 1 
  maxDims = R1dims     # Coarsening, R2dims < R1dims
  R2dims  = (rr * R1dims).astype(int)
  s2 = (2**(2*n))   # Scale factor. If we coarsen, the R1 values are appended. Same value to 2^2n cells.
  i1 = np.arange(maxDims[0]); j1 = np.arange(maxDims[1])
  i2 = np.arange(maxDims[0]); j2 = np.arange(maxDims[1])

# Modify the integer list for refining/coarsening
i1/=dr1; j1/=dr1
i1 = i1.astype(int);  j1 = j1.astype(int)

i2*=fr2            ;  j2*=fr2
i2 = i2.astype(int);  j2 = j2.astype(int)

#print ' i2 = {} '.format(i2)
#print ' i1 = {} '.format(i1)

# Create the output array. 
R2 = np.zeros( R2dims, float )

if( n > 0 ):
  R2[i2, j2] += R1[i1,j1]
else:
  for k in xrange(maxDims[0]):
    for l in xrange(maxDims[1]):
      R2[ i2[k], j2[l] ] +=  R1[ i1[k] ,j1[l] ]

R1 = None
R2 *= s2
Rdict2 = {'R' : R2, 'LocalOrig' : R1Orig, 'dPx' : dPx2}

if( not args.printOnly ):
  fx = open( fileOut , 'w' )
  np.savetxt(fx,np.round(R2),fmt='%g')
  fx.close()
  dPx2 = dPx1/rr
  saveTileAsNumpyZ( fileOut, Rdict2 )
  
if( args.printOn or args.printOnly ):
  fig = plt.figure(num=1, figsize=(9.,9.))
  fig = addImagePlot( fig, R2 , fileOut )
  
  plt.show()
