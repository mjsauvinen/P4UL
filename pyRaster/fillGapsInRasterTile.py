#!/usr/bin/env python3
import sys
import argparse
import numpy as np
from utilities import writeLog
from mapTools import readNumpyZTile, interpolateOverNans, saveTileAsNumpyZ
from plotTools import  addImagePlot 
import matplotlib.pyplot as plt
import scipy.ndimage as sn # contains the filters
'''
Description:


Author: Mikko Auvinen
        mikko.auvinen@fmi.fi
        Finnish Meteorological Institute
'''
#==========================================================#
parser = argparse.ArgumentParser(prog='fillGapsInRasterTile.py')
parser.add_argument("-f", "--filename",type=str, help="Name of the input raster data file.")
parser.add_argument("-fo", "--fileout",type=str, help="Name of output raster file.")
parser.add_argument("-Nbd","--Nbidial", type=int, default=0,\
  help="Number of binary dialations of nan values.")
parser.add_argument("-b","--beta", type=float, nargs=2, default=[0.5,0.5],\
  help="Blending factors for the X and Y interpolated rasters. Default=[0.5,0.5].")
parser.add_argument("-p", "--printOn", action="store_true", default=False,\
  help="Print the resulting raster data.")
parser.add_argument("-pp", "--printOnly", help="Only print the resulting data. Don't save.",\
  action="store_true", default=False)
args = parser.parse_args()
writeLog( parser, args, args.printOnly )
#==========================================================#
filename= args.filename
fileout = args.fileout
Nbd     = args.Nbidial
beta    = np.array( args.beta )
printOn    = args.printOn
printOnly  = args.printOnly
#==========================================================#

# Read the raster tile to be processed.
Rdict = readNumpyZTile(filename)
R = Rdict['R']
Rdims = np.array(np.shape(R))
ROrig = Rdict['GlobOrig']
print(' Rdims = {} '.format(Rdims))
print(' ROrig = {} '.format(ROrig))

idm = np.isnan(R)
Rt2 = R.copy()   # Make a copy

if( Nbd > 0 ):
  for i in range(Nbd):
    idm = sn.binary_dilation(idm)
    print('iter {}: number of nonzeros = {} '.format(i, np.count_nonzero(idm) ))

Rt2[idm] = np.nan
Rt1 = interpolateOverNans( Rt2.copy() )
Rt2 = interpolateOverNans( Rt2.T )
Rt1 = Rt1.reshape( R.shape ); Rt2 = Rt2.reshape( R.shape[::-1] ).T

R = beta[0]*Rt1 + beta[1]*Rt2 

if( printOn or printOnly ):
  figDims = 13.*(Rdims[::-1].astype(float)/np.max(Rdims))
  #print('Sum = {}'.format(np.sum(Rf)))
  fig = plt.figure(num=1, figsize=figDims)
  fig = addImagePlot( fig, R, fileout )
  plt.show()


if( not printOnly ):
  Rdict['R'] = np.round( R , decimals=1 )
  saveTileAsNumpyZ( fileout, Rdict )
