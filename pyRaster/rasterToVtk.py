#!/usr/bin/env python
import argparse
import numpy as np
from mapTools import *
from plotTools import addContourf
from utilities import vtkWriteDataStructured2d, writeLog
import matplotlib.pyplot as plt
''' 
Description:


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='rasterToSTL.py')
parser.add_argument("-f", "--filename",type=str, help="Name of the .npz data file.")
parser.add_argument("-fo", "--fileout",type=str, help="Name of output .vtk file.",\
  default="topography.vtk")
parser.add_argument("-p", "--printOn", help="Print the resulting raster data.",\
  action="store_true", default=False) 
parser.add_argument("-pp", "--printOnly", help="Only print the resulting data. Don't save.",\
  action="store_true", default=False) 
args = parser.parse_args() 
writeLog( parser, args, args.printOnly )
#==========================================================#

# Rename, that's all.
filename  = args.filename
fileout   = args.fileout
printOn   = args.printOn
printOnly = args.printOnly 


# Read in the raster data for mesh.
# NOTE! The Ry and Rx coords are in top left [N,E] format. 
# If y-origin is at zero, the Ry-values run downward into neg. numbers.
Rdict = readNumpyZTileForMesh( filename )
R = Rdict['R']
Ry = Rdict['rowCoords']
Rx = Rdict['colCoords']
Rdims = np.array(np.shape(R))
ROrig = Rdict['GlobOrig']
dPx = entry2Int( Rdict['dPx'] )
Rdict = None

if( Ry[-1] < Ry[0] ): 
  Ry *= -1.  # Make values run into positive direction.
  R = R[::-1,:]

dPx = entry2Int( dPx )

print ' Rdims = {} '.format(Rdims)
print ' ROrig = {} '.format(ROrig)

X, Y = np.meshgrid( Rx, Ry )

if( ROrig[0] == 0 ):
  print ' Resetting y-origin. '
  ROrig[0] = Rdims[0]*dPx


if( not printOnly ):
  vtkWriteDataStructured2d( None , X, Y, R, fileout, 'Z' )

if( printOnly or printOn ):
  C = addContourf( X, Y, R, 'R', fileout+' raster data.' )
  plt.show()


