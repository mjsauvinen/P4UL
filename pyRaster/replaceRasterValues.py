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
def replaceMask( Rx, px1, px2, lineOpt ):
  idm = np.zeros( np.shape(Rx), bool )
  
  if( lineOpt ):
    N = np.max( np.abs(px2-px1) )
    jrows = np.round( np.linspace(px1[0],px2[0],N) ).astype(int)
    icols = np.round( np.linspace(px1[1],px2[1],N) ).astype(int)
    for i in range(N):
      idm[ jrows[i] , icols[i] ] = True
  else:
    idm[ px1[0]:px2[0] , px1[1]:px2[1] ] = True 

  return idm


#==========================================================#
parser = argparse.ArgumentParser(prog='replaceRasterValues.py')
parser.add_argument("-f", "--filename",type=str, help="Name of the input raster file.")
parser.add_argument("-fr", "--filereplace",type=str, default=None,\
  help="(Optional) Name of raster file from which replacement values are obtained.")
parser.add_argument("-fo", "--fileout",type=str, help="Name of output raster file.")
parser.add_argument("-p1","--pixels1", type=int,nargs=2,default=[None,None],\
  help="Pixel ids [N,E] for the top left or start of line if --line opt is used.")
parser.add_argument("-p2","--pixels2",type=int,nargs=2,default=[None,None],\
  help="Pixel ids [N,E] for the bottom right or end of line if --line opt is used.")
parser.add_argument("-v","--value", help="Replacement value. Default=nan",\
  type=float, default=np.nan)
parser.add_argument("-l", "--line", action="store_true", default=False,\
  help="Replace values along line from -p1 to -p2.")
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

p1       = np.array( args.pixels1 )
p2       = np.array( args.pixels2 )
val      = args.value        # Replacement value
gtval    = args.gt           # value greater than which will be replaced by [val]
ltval    = args.lt           # value less than which will be replaced by [val]
filename = args.filename
filereplace = args.filereplace
fileout  = args.fileout
lineMode = args.line 


if( not lineMode ):
  NonesExist = (list(p1).count(None) != 0) or (list(p2).count(None) != 0) 
  WrongOrder = any( p1 > p2 )

  if( NonesExist or WrongOrder ):
    sys.exit('Error: p1 = {} or p2 = {} incorrectly specified. Exiting ...'.format(p1,p2))



# Read the raster tile to be processed.
Rdict = readNumpyZTile( filename )
R = Rdict['R']
Rdims = np.array(np.shape(R))
ROrig = Rdict['GlobOrig']
print(' Rdims = {} '.format(Rdims))
print(' ROrig = {} '.format(ROrig))

print(' Value at top left: {} '.format(R[p1[0],p1[1]]))
print(' Value at bottom right: {} '.format(R[p2[0]-1,p2[1]-1]))

idR = replaceMask( R , p1, p2, lineMode )

if( filereplace is not None ):
  Rrdict = readNumpyZTile( filereplace )
  Rr = Rrdict['R']
  Rrdims = np.array( np.shape(Rr) )
  if( any( Rrdims != Rdims ) ):
    sys.exit(' Rasters {} and {} are not the same size. Exiting ...'.format(filename, filereplace))
  print(' Values from {} are used to replace values in {}.'.format(filereplace, filename))
  R[idR] = Rr[idR]
  Rr = None
  
elif( (gtval is not None) or (ltval is not None) ):
  idx = np.zeros( R.shape, bool )
  if( gtval is not None ):
    idx = (R > gtval) * idR
    print(' {} values > {} will be replaced.'.format(np.count_nonzero(idx),gtval)) 
    R[idx] = val
    idx[:,:] = False  # Reset
  if( ltval is not None ):
    idx = (R < ltval) * idR
    print(' {} values < {} will be replaced.'.format(np.count_nonzero(idx),ltval))
    R[idx] = val 
else:
  R[idR] = val


Rdict['R'] = R

if( not args.printOnly ):
  saveTileAsNumpyZ( fileout, Rdict )

if( args.printOn or args.printOnly ):
  R[idR] = 0.
  figDims = 13.*(Rdims[::-1].astype(float)/np.max(Rdims))
  fig = plt.figure(num=1, figsize=figDims)
  fig = addImagePlot( fig, R, fileout )

  plt.show()


R = None
