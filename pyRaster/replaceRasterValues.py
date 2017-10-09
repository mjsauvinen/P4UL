#!/usr/bin/env python
import sys
import argparse
import numpy as np
from mapTools import *
from utilities import filesFromList, writeLog
from plotTools import addImagePlot, addContourf, addScatterPlot
import matplotlib.pyplot as plt
'''
Description:


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi
        University of Helsinki &
        Finnish Meteorological Institute
'''


#==========================================================#
parser = argparse.ArgumentParser(prog='replaceRasterValues.py')
parser.add_argument("-f", "--filename",type=str, help="Name of the input raster file.")
parser.add_argument("-fo", "--fileOut",type=str, help="Name of output raster file.")
parser.add_argument("-p1","--pixels1", help="Pixel ids [N,E] for the top left.",\
  type=int,nargs=2,default=[None,None])
parser.add_argument("-p2","--pixels2", help="Pixel ids [N,E] for the bottom right.",\
  type=int,nargs=2,default=[None,None])
parser.add_argument("-v","--value", help="Replacement value. Default=0.",\
  type=float, default=0.)
parser.add_argument("-p", "--printOn", help="Print the resulting raster data.",\
  action="store_true", default=False)
parser.add_argument("-pp", "--printOnly", help="Only print the resulting data. Don't save.",\
  action="store_true", default=False)
args = parser.parse_args()
writeLog( parser, args, args.printOnly )
#==========================================================#

p1      = args.pixels1 
p2      = args.pixels2 
val     = args.value
filename= args.filename
fileOut = args.fileOut

NonesExist = (p1.count(None) != 0) or (p2.count(None) != 0) 
p1 = np.array( p1 ); p2 = np.array( p2 )
WrongOrder = any( p1 > p2 )

if( NonesExist or WrongOrder ):
  print('Error: pixels1 = {} or pixels2 = {} incorrectly specified. Exiting ... '.format(p1,p2))
  sys.exit(1)



# Read the raster tile to be processed.
Rdict = readNumpyZTile( filename )
R = Rdict['R']
Rdims = np.array(np.shape(R))
ROrig = Rdict['GlobOrig']
print(' Rdims = {} '.format(Rdims))
print(' ROrig = {} '.format(ROrig))

print(' Value at top left: {} '.format(R[p1[0],p1[1]]))
print(' Value at bottom right: {} '.format(R[p2[0],p2[1]]))

R[p1[0]:p2[0],p1[1]:p2[1]] = val 
Rdict['R'] = R

if( not args.printOnly ):
  saveTileAsNumpyZ( fileOut, Rdict )

if( args.printOn or args.printOnly ):
  figDims = 13.*(Rdims[::-1].astype(float)/np.max(Rdims))
  fig = plt.figure(num=1, figsize=figDims)
  fig = addImagePlot( fig, R, fileOut )

  plt.show()


R = None

