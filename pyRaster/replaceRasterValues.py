#!/usr/bin/env python3
import sys
import argparse
import numpy as np
from mapTools import *
from utilities import filesFromList, writeLog
from plotTools import addImagePlot
import matplotlib.pyplot as plt
import scipy.ndimage as sn
'''
Description:


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi
        University of Helsinki &
        Finnish Meteorological Institute
'''
#==========================================================#
def replaceMask( Rx, px1, px2, lineOpt, gtval, ltval, Nbd, Nans=False):
  '''Return the indices that will be replaced

  This function handles all the location specific work. Line mode and rectangle
  mode are complementary. The gt and lt modes can be combined with other
  location selection modes.

  '''
  
  # Treat --fillNans first
  if( Nans ):
    idm = np.isnan(Rx) 
    print('N={} nan values to be filled.'.format( np.count_nonzero(idm) ))
    return idm
  
  # Continue with conventional treatment
  idm = np.zeros( np.shape(Rx), bool )
  
  if( lineOpt ):
    N = np.max( np.abs(px2-px1) )
    jrows = np.round( np.linspace(px1[0],px2[0],N) ).astype(int)
    icols = np.round( np.linspace(px1[1],px2[1],N) ).astype(int)
    for i in range(N):
      idm[ jrows[i] , icols[i] ] = True
  else:
    idm[ px1[0]:px2[0] , px1[1]:px2[1] ] = True

  if( gtval is not None ):
    idm = (Rx > gtval) * idm
    print(' {} values > {} will be replaced with selection mask.'.format(np.count_nonzero(idm),gtval))

  if( ltval is not None ):
    idm = (Rx < ltval) * idm
    print(' {} values < {} will be replaced with selection mask.'.format(np.count_nonzero(idm),ltval))

  if(Nbd > 0):
    print('Applying binary dilation to selection mask.')
    for i in range(Nbd):
      idm = sn.binary_dilation(idm)
      print('iter {}: number of masked points = {} '.format(i, np.count_nonzero(idm) ))

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
parser.add_argument("-v","--value", type=float, default=None, \
  help="Replacement value. If None, other actions may apply. Default=None")
parser.add_argument("-c","--coef", type=float, default=1.0,\
  help="Multiplication coefficient with which masked values are multiplied. Default=1.")
parser.add_argument("-na", "--nans", action="store_true", default=False,\
  help="Use Nans as replacement values. Overrides the --value specified above.")
parser.add_argument("-nf", "--fillNans", help="Fill nans with replace value.",\
  action="store_true", default=False)
parser.add_argument("-l", "--line", action="store_true", default=False,\
  help="Replace values along line from -p1 to -p2.")
parser.add_argument("-p", "--printOn", help="Print the resulting raster data.",\
  action="store_true", default=False)
parser.add_argument("-pp", "--printOnly", action="store_true", default=False,\
  help="Only print the resulting data. Don't save.")
parser.add_argument("-gt", "--gt", type=float, default=None,\
  help="Replace values greater than the given value.")
parser.add_argument("-lt", "--lt", type=float, default=None,\
  help="Replace values less than the given value.")
parser.add_argument("-Nbd","--Nbidil", type=int, default=0,\
  help="Number of binary dilations of nan values.")
parser.add_argument("-a", "--allPoints", action="store_true", default=False,\
  help="Select all points in raster for processing. Overrides p1 and p2.")
args = parser.parse_args()
writeLog( parser, args, args.printOnly )
#==========================================================#

p1          = np.array( args.pixels1 )
p2          = np.array( args.pixels2 )
val         = args.value        # Replacement value
useNans     = args.nans
fillNans    = args.fillNans
cf          = args.coef         # Multiplication coefficient
gtval       = args.gt           # value greater than which will be replaced by [val]
ltval       = args.lt           # value less than which will be replaced by [val]
filename    = args.filename
filereplace = args.filereplace
fileout     = args.fileout
lineMode    = args.line 
Nbd         = args.Nbidil
allPoints   = args.allPoints   

# Read the raster tile to be processed.
Rdict = readNumpyZTile( filename )
R = Rdict['R']
Rdims = np.array(np.shape(R))
ROrig = Rdict['GlobOrig']

if( allPoints ):
  p1=np.zeros(2,int)
  p2=Rdims
  print('Selecting all points for processing')

if( not lineMode ):
  try:
    WrongOrder = any( p1 > p2 )
    if( WrongOrder ):
      sys.exit('Error: p1 = {} or p2 = {} in wrong order. Exiting ...'.format(p1,p2))
  except (TypeError):
    sys.exit('Error: p1 = {} or p2 = {} incorrectly specified. Exiting ...'.format(p1,p2))

print(' Rdims = {} '.format(Rdims))
print(' ROrig = {} '.format(ROrig))

print(' Value at top left: {} '.format(R[p1[0],p1[1]]))
print(' Value at bottom right: {} '.format(R[p2[0]-1,p2[1]-1]))

idR = replaceMask( R , p1, p2, lineMode, gtval, ltval, Nbd, fillNans )

if( useNans ):
  val = np.nan
  R = R.astype(float)

if( filereplace is not None ):
  Rrdict = readNumpyZTile( filereplace )
  Rr = Rrdict['R']
  Rrdims = np.array( np.shape(Rr) )
  if( any( Rrdims != Rdims ) ):
    sys.exit(' Rasters {} and {} are not the same size. Exiting ...'.format(filename, filereplace))
  print(' Values from {} are used to replace values in {}.'.format(filereplace, filename))
  R[idR] = Rr[idR]
  Rr = None
elif( val is not None):
  R[idR] = val
  
if( not useNans ): 
  R[idR] = R[idR]*cf

Rdict['R'] = R

if( not args.printOnly ):
  saveTileAsNumpyZ( fileout, Rdict )

if( args.printOn or args.printOnly ):
  #R[idR] = 0.
  figDims = 13.*(Rdims[::-1].astype(float)/np.max(Rdims))
  fig = plt.figure(num=1, figsize=figDims)
  fig = addImagePlot( fig, R, fileout )

  plt.show()

R = None
