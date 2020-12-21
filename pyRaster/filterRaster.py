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
args = parser.parse_args()
filename= args.filename
fileout = args.fileout
#=input=======================================================================#
Rdict = readNumpyZTile(filename)
R = Rdict['R']
Rdims = np.array(np.shape(R))
ROrig = Rdict['GlobOrig']
print(' Rdims = {} '.format(Rdims))
print(' ROrig = {} '.format(ROrig))
#=filtering===================================================================#
R=sn.median_filter(R,size=5)
#=output======================================================================#
Rdict['R'] = np.round( R , decimals=1 )
saveTileAsNumpyZ( fileout, Rdict )
