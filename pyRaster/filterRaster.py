#!/usr/bin/env python3
import argparse
import numpy as np
import scipy.ndimage as sn # contains the filters
from mapTools import readNumpyZTile, saveTileAsNumpyZ
'''A script for filtering rasters. 

Jukka-Pekka Keskinen
jukka-pekka.keskinen@fmi.fi
Finnish Meteorological Insitute
2020
'''
#=arg parser==================================================================#
parser = argparse.ArgumentParser(prog='filterRaster.py',
                                 description='Filter an npz raster.')
parser.add_argument("-f", "--filename",type=str, 
                    help="Name of the input raster data file.")
parser.add_argument("-fo", "--fileout",type=str, 
                    help="Name of output raster file.")
parser.add_argument("-n","--size", type=int, nargs=2, default=[5,5],
                    help="Size (x,y) of filtering window. Default=5×5.")
parser.add_argument("-s","--sigma", type=float, default=None, help="Standard" 
     " deviation (σ). Invokes the"+"\033[1m"+" Gaussian filter"+"\033[0m"+".")
args = parser.parse_args()
filename= args.filename
fileout = args.fileout
fsize = np.array( args.size)
sigma = args.sigma
#=input=======================================================================#
Rdict = readNumpyZTile(filename)
R = Rdict['R']
Rdims = np.array(np.shape(R))
ROrig = Rdict['GlobOrig']
print(' Rdims = {} '.format(Rdims))
print(' ROrig = {} '.format(ROrig))
#=filtering===================================================================#
if (sigma is not None):
  print('Applying the gaussian filter.')
  R=sn.gaussian_filter(R,sigma)
else:
  print('Applying the median filter.')
  R=sn.median_filter(R,size=fsize)
#=output======================================================================#
Rdict['R'] = np.round( R , decimals=1 )
saveTileAsNumpyZ( fileout, Rdict )
