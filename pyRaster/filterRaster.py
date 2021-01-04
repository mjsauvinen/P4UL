#!/usr/bin/env python3
import argparse
import sys
import numpy as np
import scipy.ndimage as sn # contains the filters
from mapTools import readNumpyZTile, saveTileAsNumpyZ
'''A script for filtering rasters. 

Jukka-Pekka Keskinen
jukka-pekka.keskinen@fmi.fi
Finnish Meteorological Insitute
2020-2021
'''
#=arg parser==================================================================#
parser = argparse.ArgumentParser(prog='filterRaster.py',
    description='Filter a P4UL raster. Only one filter is applied at a time.')
parser.add_argument("-f", "--filename",type=str, 
                    help="Name of the input raster data file.")
parser.add_argument("-fo", "--fileout",type=str, 
                    help="Name of output raster file.")
parser.add_argument("-n","--size", type=int, nargs=2, default=[5,5],
                    help="Size (x,y) of filtering window. Default=5×5.")
parser.add_argument("-s","--sigma", type=float, default=None, help="Standard" 
     " deviation (σ). Invokes the"+"\033[1m"+" Gaussian filter"+"\033[0m"+"."
                    " Does not use filter size even if specified.")
parser.add_argument("-t","--maximum",action="store_true", default=False, 
    help="Invoke the "+"\033[0m"+"maximum filter"+"\033[0m"+".")
parser.add_argument("-m","--median",action="store_true", default=False, 
    help="Invoke the "+"\033[0m"+"median filter"+"\033[0m"+".")
parser.add_argument("-b","--minimum",action="store_true", default=False, 
    help="Invoke the "+"\033[0m"+"minimum filter"+"\033[0m"+".")
parser.add_argument("-r","--rank", type=int, default=None, help="Rank to "
                    " use in the "+"\033[0m"+"rank filter"+"\033[0m"+".")
args = parser.parse_args()
filename= args.filename
fileout = args.fileout
fsize = np.array( args.size)
sigma = args.sigma
maxi  = args.maximum
medi  = args.median
mini  = args.minimum 
rank  = args.rank
#=input=======================================================================#
Rdict = readNumpyZTile(filename)
R = Rdict['R']
Rdims = np.array(np.shape(R))
ROrig = Rdict['GlobOrig']
print(' Rdims = {} '.format(Rdims))
print(' ROrig = {} '.format(ROrig))
#=filtering===================================================================#
if (sigma is not None):
  print(' Applying the gaussian filter with σ='+str(sigma)+'.')
  R=sn.gaussian_filter(R,sigma)
elif maxi:
  print(' Applying the maximum filter.')
  R=sn.maximum_filter(R,size=fsize)
elif medi:
  print(' Applying the median filter.')
  R=sn.median_filter(R,size=fsize)
elif mini:
  print(' Applying the minumum filter.')
  R=sn.minimum_filter(R,size=fsize)
elif (rank is not None):
  print(' Applying the rank filter. Selecting values at rank '+str(rank)+'/'+
        str(np.prod(fsize))+' from the filter window.')
  R=sn.rank_filter(R,rank,size=fsize)
else:
  print(' No filter specified. Nothing to do here.')
  sys.exit()
#=output======================================================================#
Rdict['R'] = np.round( R , decimals=1 )
saveTileAsNumpyZ( fileout, Rdict )
