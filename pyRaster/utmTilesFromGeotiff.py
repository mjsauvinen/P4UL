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
parser = argparse.ArgumentParser(prog='utmTilesFromGeoTiff.py')
parser.add_argument("-f", "--filename", type=str, help="Name of the target .tif file", \
  default=None)
parser.add_argument("-b", "--bandSelect", help="Raster Band Selection. Default=1.",\
  action="store_true", default=False)
parser.add_argument("-p", "--printOn", help="Print the extracted tile.",\
  action="store_true", default=False)
parser.add_argument("-pp", "--printOnly", help="Only print the extracted tile. Don't save.",\
  action="store_true", default=False)
parser.add_argument("-t", "--utmTile", help="Utm tile code (2 chars). For example: P5",\
   type=str, default=None)
parser.add_argument("-s", "--scale", type=float,\
   help="Scale factor for the output. Default=1.", default=1.)
args = parser.parse_args()
writeLog( parser, args, args.printOnly )
#==========================================================#
# Rename ... that's all.
filename = args.filename 
utmTile  = args.utmTile
bandSelect = args.bandSelect
scale      = args.scale

# If the geotiff file has not been specified, ask specifically.
if( not filename ):
  fileNos, fileList = filesFromList( "*.tif" )
  try: filename = fileList[0]
  except: sys.exit('Could not obtain a valid GEOTIFF file. Exiting ...')


dataset = openGeoTiff( filename )
XOrig, resolution = getGeoTransform( dataset ) # Both 2d arrays. XOrig is Top Left!

nBands = numberOfRasterBands( dataset, True) # PrintOn = True/False

# Ask for the band ID if the user so wants (given that it makes sense)
ib = selectBand( nBands, bandSelect , 1 ) # last argument is for default

rb = getRasterBand( dataset, ib )
printRasterBandStatistics(rb)

#R, Rdims = readAsNumpyArray( rb )  # This is heavy and undesirable way to go!

if( utmTile is not None ):
  # Define the bottom left (BL) origin of this tile, for example, P5.
  # This is redundant if the getGeoTransform() works.
  XOrig_check = UtmReference( utmTile )
  print " XOrig: {} vs. {} ".format( XOrig, XOrig_check )

# Determine the resolution of the raster data [dn,de].
#resolution_check = UtmTileDims() / Rdims
#print ' Pixel resolution (m): {} vs. {}'.format(resolution, resolution_check)
Rdict = extractSubTile( rb, utmTile, XOrig, resolution)

Rdict['dPx'] = resolution
Rdict['rotation'] = 0.

Rdict['R'][Rdict['R']>32765] = 0.

if( args.scale != 1.):
  Rdict['R']=Rdict['R']*args.scale

if( utmTile is not None):
  fileout = utmTile
else:
  fileout = filename.split('.')[0]


if( not args.printOnly ):
  saveTileAsNumpyZ( fileout , Rdict )


if( args.printOnly or args.printOn ):
  Rdims = np.array( Rdict['R'].shape )
  figDims = 13.*(Rdims[::-1].astype(float)/np.max(Rdims))
  fig = plt.figure(num=1, figsize=figDims)
  fig = addImagePlot( fig, Rdict['R'], fileout )
  plt.show()
