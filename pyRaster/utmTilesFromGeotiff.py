#!/usr/bin/env python
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
parser = argparse.ArgumentParser(prog='utmTilesFromGeoTiff.py')
parser.add_argument("-f", "--filename",type=str, help="Name of the target .tif file", \
  default=None)
parser.add_argument("-b", "--bandSelect", help="Raster Band Selection. Default=1.",\
  action="store_true", default=False) 
parser.add_argument("-p", "--printOn", help="Print the extracted tile.",\
  action="store_true", default=False) 
parser.add_argument("-pp", "--printOnly", help="Only print the extracted tile. Don't save.",\
  action="store_true", default=False) 
parser.add_argument("-t", "--utmTile", help="Utm tile code (2 chars). For example: P5",\
   type=str, default=False)
parser.add_argument("-s", "--scale",type=float,\
   help="Scale factor for the output. Default=1.", default=1.)
args = parser.parse_args() 
writeLog( parser, args, args.printOnly )
#==========================================================#



# If the geotiff file has not been specified, ask specifically.
if( not args.filename ):
  fileNos, fileList = filesFromList( "*.tif" )
  try: args.filename = fileList[0]
  except: sys.exit('Could not obtain a valid GEOTIFF file. Exiting ...')


dataset     = openGeoTiff( args.filename )
XOrig, resolution = getGeoTransform( dataset ) # Both 2d arrays. XOrig is Top Left!

nBands  = numberOfRasterBands( dataset, True) # PrintOn = True/False

# Ask for the band ID if the user so wants (given that it makes sense)
ib = selectBand( nBands, args.bandSelect , 1 ) # last argument is for default

rb = getRasterBand( dataset, ib )
printRasterBandStatistics(rb)

#R, Rdims = readAsNumpyArray( rb )  # This is heavy and undesirable way to go!

# Define the bottom left (BL) origin of this tile, for example, P5.
# This is redundant if the getGeoTransform() works.
XOrig_check = UtmReference( args.utmTile )
print " XOrig: {} vs. {} ".format( XOrig, XOrig_check )

# Determine the resolution of the raster data [dn,de].
#resolution_check = UtmTileDims() / Rdims
#print ' Pixel resolution (m): {} vs. {}'.format(resolution, resolution_check)

R, Rdims, ROrig = extractSubTile( rb, args.utmTile, XOrig, resolution)

R[R>32765] = 0.

if( args.scale != 1.):
  R*=args.scale

if( not args.printOnly ):
  saveTileAsNumpyZ( args.utmTile, R, Rdims, ROrig, resolution )


if( args.printOnly or args.printOn ):
  figDims = 13.*(Rdims[::-1].astype(float)/np.max(Rdims))
  fig = plt.figure(num=1, figsize=figDims)
  fig = addImagePlot( fig, R, args.utmTile )
  plt.show()

