#!/usr/bin/env python3
import sys
import argparse
import numpy as np
from gdalTools import *
from mapTools import saveTileAsNumpyZ
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
parser = argparse.ArgumentParser(prog='geoTif2NumpyZ.py')
parser.add_argument("-f", "--filename", type=str, help="Name of the target .tif file", \
  default=None)
parser.add_argument("-nd", "--ndecimals",type=int,\
   help="Number of decimal places. Default = 1", default=1) 
parser.add_argument("-b", "--bandSelect", help="Raster Band Selection.",\
  action="store_true", default=False)
parser.add_argument("-p", "--printOn", help="Print the extracted tile.",\
  action="store_true", default=False)
parser.add_argument("-pp", "--printOnly", help="Only print the extracted tile. Don't save.",\
  action="store_true", default=False)
args = parser.parse_args()

#==========================================================#
# Rename ... that's all.
filename = args.filename 
bandSelect = args.bandSelect
printOn    = args.printOn
printOnly  = args.printOnly
ndecimals  = args.ndecimals

dataset = openGeoTiff( filename )
nBands = numberOfRasterBands( dataset, True) # PrintOn = True/False
ROrig, dPx = getGeoTransform( dataset ) # Both 2d arrays. XOrig is Top Left!

# Ask for the band ID if the user so wants (given that it makes sense)
ib = selectBand( nBands, bandSelect , 1 ) # last argument is for default
rb = getRasterBand( dataset, ib )
printRasterBandStatistics(rb)

# Read the raster
R = readAsNumpyArray( rb )

print(' Number of non-zero values = {}'.format(np.count_nonzero(R)))

# Construct the Raster dict 
Rdict = dict()
Rdict['dPx'] = dPx
Rdict['rotation'] = 0.
Rdict['R'] = np.round( R , decimals=ndecimals )
Rdict['GlobOrig'] = ROrig

if( not printOnly ):
  saveTileAsNumpyZ( filename.replace('.tif','.npz'), Rdict )
  Rdict = None


if( printOnly or printOn ):
  Rdims = np.array( R.shape )
  figDims = 13.*(Rdims[::-1].astype(float)/np.max(Rdims))
  fig = plt.figure(num=1, figsize=figDims)
  fig = addImagePlot( fig, R, filename )
  plt.show()


