#!/usr/bin/env python
import sys
import argparse
import numpy as np
from mapTools import *
from utilities import filesFromList, writeLog
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
parser = argparse.ArgumentParser(prog='rasterTileSuperposition.py')
parser.add_argument("-f1", "--file1",type=str, help="Name of raster data file (1).")
parser.add_argument("-f2", "--file2",type=str, help="Name of raster data file (2).",\
  default='Rout')
parser.add_argument("-s1", "--scale1",type=float,\
  help="Scaling factor for the data in file 1. Default=1.", default=1.)   
parser.add_argument("-s2", "--scale2",type=float,\
  help="Scaling factor for the data in file 2. Default=1.", default=1.)  

helpFlt = ''' Filter type and its associated number for file 1.
 The raster data will be migrated into the max(f1,f2) resolution 
 before filtering is applied. 
 Available filters: median, percentile, rank, gaussian, local. 
 Entering \"user, num\" allows the user to specify <num> different filters consecutively.
 Example entry: median 5'''
parser.add_argument("-ft1","--filter1",type=str,nargs=2,default=[None,None],\
  help=helpFlt)
parser.add_argument("-ft2","--filter2",type=str,nargs=2,default=[None,None],\
  help="See help for -ft1 above.")
parser.add_argument("-fo", "--fileOut",type=str, help="Name of output .npz file.")
parser.add_argument("-p", "--printOn", help="Print the resulting raster data.",\
  action="store_true", default=False) 
parser.add_argument("-pp", "--printOnly", help="Only print the resulting data. Don't save.",\
  action="store_true", default=False) 
args = parser.parse_args() 
writeLog( parser, args, args.printOnly )
#==========================================================#

# Rename some args for brevity.
s1 = args.scale1;       s2 = args.scale2
flt1 = args.filter1;    flt2 = args.filter2
printOn = args.printOn; printOnly = args.printOnly

figN = 1      # Figure number ... to be appended on demand.

# Read in the data.
dataOnly = False
R1, R1dims, R1Orig, dPx1 = readNumpyZTile(args.file1, dataOnly)
R2, R2dims, R2Orig, dPx2 = readNumpyZTile(args.file2, dataOnly)


dPx1 = entry2Int( dPx1 ); dPx2 = entry2Int( dPx2 )
dPc = max( dPx1, dPx2 )   # Coarser resolution.
dPf = min( dPx1, dPx2 )   # Finer resolution.

if( (R1Orig == R2Orig).all() ):
  print ' Excellent! The origos match.'
else:
  print ' The tiles do not have identical origos. Exiting.'
  sys.exit(1)

maxDims = np.array([ max(R1dims[0],R2dims[0]) , max(R1dims[1],R2dims[1]) ]).astype(int)

# Compute resolution ratios (rr). One of the two is going to be 1.
rr1 = int(dPc/dPx2)
rr2 = int(dPc/dPx1)

i1,j1 = np.ogrid[0:maxDims[0], 0:maxDims[1]]
i1/=rr1; j1/=rr1

i2,j2 = np.ogrid[0:maxDims[0], 0:maxDims[1]]
i2/=rr2; j2/=rr2

print ' dims 1: {}, {} '.format(np.max(i1),np.max(j1))
print ' dims 2: {}, {} '.format(np.max(i2),np.max(j2))

# Initialize the new storage arrays which are of fine resolution.
Rt1 = np.zeros( maxDims, float )
Rt2 = np.zeros( maxDims, float )

# Apply filtering if desired and perform the superposition with appropriate scaling.  
Rt1 = filterAndScale(Rt1, R1, flt1, s1 , i1, j1)
Rt2 = filterAndScale(Rt2, R2, flt2, s2 , i2, j2)
Rt = Rt1 + Rt2

# Print the filtered raster maps.
if( printOn or printOnly ):
  if( flt1.count(None) == 0):
    fg1 = plt.figure(num=figN, figsize=(9.,9.)); figN+=1
    fg1 = addImagePlot(fg1, Rt2, args.file1, gridOn=True )
  
  if( flt2.count(None) == 0):
    fg2 = plt.figure(num=figN, figsize=(9.,9.)); figN+=1
    fg2 = addImagePlot(fg2, Rt2, args.file2, gridOn=True )

Rt1 = Rt2 = None

if( not printOnly ):
  saveTileAsNumpyZ( args.fileOut, Rt, maxDims, R1Orig, np.array([dPf,dPf]) )
  
if( printOn or printOnly ):
  fig = plt.figure(num=figN, figsize=(9.,9.)); figN+=1
  fig = addImagePlot( fig, Rt[::2,::2], args.file1+' + '+args.file2, gridOn=True )
  plt.show()

Rt = None
'''
print ' {}  {} '.format( dPc/dPx1 , dPc/dPx2 )
print ' {}  {} '.format( R1Orig ,R2Orig )
print ' {}  {} '.format( R1dims ,R2dims )
'''




