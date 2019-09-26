#!/usr/bin/env python3
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
parser = argparse.ArgumentParser(prog='rasterToAscii.py')
parser.add_argument("-f", "--filename",type=str, help="Name of the comp domain data file.")
parser.add_argument("-fo", "--fileout",type=str, help="Name of output ASCII file.")
parser.add_argument("-i", "--round2Int", help="Round the output data to nearest integers.",\
  action="store_true", default=False)
parser.add_argument("-p", "--printOn", help="Print also the raster data.",\
  action="store_true", default=False)
parser.add_argument("-pp", "--printOnly", help="Only print raster data. Don't write.",\
  action="store_true", default=False)
args = parser.parse_args()
writeLog( parser, args, args.printOnly )
#==========================================================#

filename   = args.filename
fileout    = args.fileout
round2Int  = args.round2Int
printOn    = args.printOn
printOnly  = args.printOnly

# Read the raster tile to be processed.
Rdict = readNumpyZTile(filename)
R = Rdict['R']
Rdims = np.array(np.shape(R))
ROrig = Rdict['GlobOrig']
print(' Rdims = {} '.format(Rdims))
print(' ROrig = {} '.format(ROrig))

if( not printOnly ):
  fx = open( fileout , 'w' )
  if( round2Int ): np.savetxt(fx,np.round(R),fmt='%g')
  else:            np.savetxt(fx,R,fmt='%g')
  fx.close()

if( args.printOn or args.printOnly ):
  figDims = 13.*(Rdims[::-1].astype(float)/np.max(Rdims))
  #print('Sum = {}'.format(np.sum(R)))
  fig = plt.figure(num=1, figsize=figDims)
  fig = addImagePlot( fig, R, fileout )

  plt.show()


R = Rf = None

