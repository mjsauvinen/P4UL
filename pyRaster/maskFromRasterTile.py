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
parser = argparse.ArgumentParser(prog='maskFromRasterTile.py')
parser.add_argument("-f", "--filename",type=str, help="Input .npz file name.")
parser.add_argument("-fo", "--fileout",type=str, help="Output .npz file name.")
parser.add_argument("-mv", "--maskvals",type=int, nargs='+', \
   help="Values used to create the mask.")
parser.add_argument("-zv", "--zeroval",type=int, help="Value used in place of zero.")
parser.add_argument("-a1", "--allone", help="All values as one [1].",\
  action="store_true", default=False)
parser.add_argument("-p", "--printOn", help="Print the numpy array data.",\
  action="store_true", default=False) 
parser.add_argument("-pp", "--printOnly", help="Only print the numpy array data. Don't save.",\
  action="store_true", default=False) 
args = parser.parse_args() 
writeLog( parser, args, args.printOnly )
#==========================================================#

# Renaming, nothing more.
filename  = args.filename
fileout   = args.fileout
mvals     = args.maskvals
allOne    = args.allone
zval      = args.zeroval
printOn   = args.printOn
printOnly = args.printOnly

# Read in the raster data.
R, Rdims, ROrig, dPx = readNumpyZTile( filename )

# Create mask raster
Rm = np.zeros( Rdims, 'uint8' )

for vx in mvals:
  fx = vx
  if( fx == 0 ): fx = zval
  if( allOne  ): fx = 1
  
  Rm += (R == vx ).astype(int) * int(fx)
  print(' Rmtype = {} '.format(Rm.dtype))

if( not printOnly ):
  print(' Writing file {} ... '.format(fileout) ) 
  saveTileAsNumpyZ( fileout, Rm, Rdims, ROrig, dPx)
  print(' ... done! ')

if( printOn or printOnly ):
  figDims = 13.*(Rdims[::-1].astype(float)/np.max(Rdims))
  pfig = plt.figure(num=1, figsize=figDims)
  pfig = addImagePlot( pfig, Rm, fileout, gridOn=True )
  plt.show()

