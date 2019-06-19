#!/usr/bin/env python
import sys
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
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

#==========================================================#
parser = argparse.ArgumentParser(prog='newRasterTile.py')
parser.add_argument("-fo", "--fileout",type=str, \
  help="Name of output Palm topography file.")
parser.add_argument("-D","--Dims", type=int, nargs=2,\
  help="Dimensions [N, E] of the raster tile.")
parser.add_argument("-O","--Orig", type=float, nargs=2, default=[0.,0.],\
  help="Global top left origin [N, E] of the raster tile. Default = [0, 0]")
parser.add_argument("-d","--dP", type=float, nargs=2,\
  help="Raster resolution [dn, de] in meters [m].")
parser.add_argument("-r","--rotation", type=float, default=0.,\
  help="Grid rotation in [rad] to be included in the raster file. Default = 0.")
parser.add_argument("-p", "--printOn", help="Print the resulting raster data.",\
  action="store_true", default=False)
parser.add_argument("-pp", "--printOnly", help="Only print the resulting data. Don't save.",\
  action="store_true", default=False)
args = parser.parse_args()
writeLog( parser, args, args.printOnly )
#==========================================================#

fileout    = args.fileout
Rdims      = args.Dims
ROrig      = args.Orig
dP         = args.dP
rotation   = args.rotation
printOn    = args.printOn
printOnly  = args.printOnly


Rdict = dict()
R = np.zeros( Rdims )
Rdict['R']   = R
Rdict['GlobOrig'] = np.array( ROrig )
Rdict['dPx']      = np.array( dP )
Rdict['gridRot']  = rotation

if( not args.printOnly ):
  saveTileAsNumpyZ( fileout, Rdict)


if( args.printOn or args.printOnly ):
  Rdims = np.array( np.shape( R ) )
  figDims = 13.*(Rdims[::-1].astype(float)/np.max(Rdims))
  fig = plt.figure(num=1, figsize=figDims)
  fig = addImagePlot( fig, R, fileout )

  plt.show()


