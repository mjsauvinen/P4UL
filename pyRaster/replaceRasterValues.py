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
parser.add_argument("-fr", "--filereplace",type=str,\
  help="(Optional) Name of raster file from which replacement values are obtained.", default=None)
parser.add_argument("-fo", "--fileout",type=str, help="Name of output raster file.")
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
parser.add_argument("-gt", "--gt", type=float, default=None,\
  help="Replace values greater than the given value.")
parser.add_argument("-lt", "--lt", type=float, default=None,\
  help="Replace values less than the given value.")
args = parser.parse_args()
writeLog( parser, args, args.printOnly )
#==========================================================#

p1      = args.pixels1 
p2      = args.pixels2 
val     = args.value        # Replacement value
gtval   = args.gt           # value greater than which will be replaced by [val]
ltval   = args.lt           # value less than which will be replaced by [val]
filename= args.filename
filereplace = args.filereplace
fileout = args.fileout

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

if( filereplace is not None ):
  Rrdict = readNumpyZTile( filereplace )
  Rr = Rrdict['R']
  Rrdims = np.array( np.shape(Rr) )
  if( any( Rrdims != Rdims ) ):
    sys.exit(' Rasters {} and {} are not the same size. Exiting ...'.format(filename, filereplace))
  print(' Values from {} are used to replace values in {}.'.format(filereplace, filename))
  R[p1[0]:p2[0],p1[1]:p2[1]] = Rr[p1[0]:p2[0],p1[1]:p2[1]]
  Rr = None
  
elif( (gtval is not None) or (ltval is not None) ):
  idx = np.zeros( R.shape, bool )
  if( gtval is not None ):
    idx[p1[0]:p2[0],p1[1]:p2[1]] = ( R[p1[0]:p2[0],p1[1]:p2[1]] > gtval )
    print(' {} values > {} will be replaced.'.format(np.count_nonzero(idx),gtval)) 
    R[idx] = val
    idx[:,:] = False  # Reset
  if( ltval is not None ):
    idx[p1[0]:p2[0],p1[1]:p2[1]] = ( R[p1[0]:p2[0],p1[1]:p2[1]] < ltval )
    print(' {} values < {} will be replaced.'.format(np.count_nonzero(idx),ltval))
    R[idx] = val 
else:
  R[p1[0]:p2[0],p1[1]:p2[1]] = val


Rdict['R'] = R

if( not args.printOnly ):
  saveTileAsNumpyZ( fileout, Rdict )

if( args.printOn or args.printOnly ):
  figDims = 13.*(Rdims[::-1].astype(float)/np.max(Rdims))
  fig = plt.figure(num=1, figsize=figDims)
  fig = addImagePlot( fig, R, fileout )

  plt.show()


R = None

