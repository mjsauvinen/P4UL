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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

def maskFromClip( Ri, mclip, mlist ):
  Rm  = np.zeros( Ri.shape )
  idx = (Ri > mclip )
  Rm[idx] = 1
  mlist.append(1)
  
  return Rm, mlist

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

def maskFromData(Ri, mlist, mxN=20):
  for im in xrange(1,mxN+1):
    if( im in Ri ):
      #print(' Mask id = {} found in file.'.format(im))
      mlist.append(im)
      
  return mlist
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

def totalArea( Rdims, dx, fname ):
  Npx  = np.prod( Rdims ) # Number of pixels
  At = Npx*np.abs(np.prod(dx))
  print('\n Total area of {} domain:\n Atot = {} m^2 \n'.format(fname,At))
  return At

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

def frontalAreas( Ri ):
  Ae = 0.
  for i in xrange( Ri.shape[0] ):
    ce = (Ri[i,1:] > Ri[i,:-1]).astype(float)
    he = (Ri[i,1:]-Ri[i,:-1]); he[(he<4.)] = 0. # Height, clip out non-buildings
    Ae += np.sum( ce * he )
  
  An = 0.
  for j in xrange( Ri.shape[1] ):
    cn = (Ri[1:,j] > Ri[:-1,j]).astype(float)
    hn = (Ri[1:,j]-Ri[:-1,j]); hn[(hn<4.)] = 0. 
    An += np.sum( cn* hn )
    
  return Ae, An

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

def maskMeanValues(Rm, Ri, mlist):
  mvals = np.zeros(np.shape(mlist))
  j = 0
  for im in mlist:
    idm = (Rm == im)
    mvals[j] = np.mean( Ri[idm] )
    print(' Mask {} mean value = {} '.format(im, mvals[j]))
    j += 1
  
  return mvals
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

def planAreaFractions( Ri, mlist ):
  Npx = np.prod( np.array(Ri.shape) )
  r = np.zeros( np.shape(mlist) )
  j = 0
  for im in mlist:
    r[j] = np.count_nonzero( Ri == im )/float( Npx )
    print(' Mask {} plan area fraction = {} '.format(im, r[j]))
    j += 1
    
  return r

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

#==========================================================#
parser = argparse.ArgumentParser(prog='maskDataEvaluation.py')
parser.add_argument("-fm", "--filemask",type=str, default=None,\
  help="Input .npz mask file name.")
parser.add_argument("-fd", "--filedata",type=str, default=None,\
  help="(Optional) Topography .npz file name.")
parser.add_argument("-Fb", "--Fafb", action="store_true", default=False, \
  help="Compute frontal area fraction (of buildings) from the topography data.")
parser.add_argument("-mx", "--maxMaskNo",type=int, default=20,\
  help="Maximum mask id value. Default=20")
parser.add_argument("-ma", "--maskAbove", help="Mask all above given value.",\
  type=int, default=None)
parser.add_argument("-p", "--printOn", help="Print the numpy array data.",\
  action="store_true", default=False)
args = parser.parse_args()
#==========================================================#

# Renaming, nothing more.
filemask  = args.filemask
filedata  = args.filedata
maxMaskNo = args.maxMaskNo
maskAbove = args.maskAbove
Fafb      = args.Fafb # frontal area fraction
printOn   = args.printOn

# Set rasters to None 
Rm = None; Rt = None


# Read in the raster data.
if( filemask ):
  Rmdict = readNumpyZTile( filemask )
  Rm = Rmdict['R']
  Rmdims = np.array(np.shape(Rm))
  RmOrig = Rmdict['GlobOrig']
  dPxm = Rmdict['dPx']

# Read the topography file if provided
if( filedata ):
  Rtdict = readNumpyZTile( filedata )
  Rt = Rtdict['R']
  Rtdims = np.array(np.shape(Rt))
  RtOrig = Rtdict['GlobOrig']
  dPxt = Rtdict['dPx']
  if( filemask and not all(Rtdims == Rmdims) ):
    sys.exit(' Error! Dimensions of topography and mask files do not agree. Exiting ...')
    
if( not (filemask or filedata) ):
  sys.exit(' No data files provided, Exiting ...')


if( filemask ):
  Atot = totalArea( Rmdims, dPxm, filemask )
else:
  Atot = totalArea( Rtdims, dPxt, filedata )


# Compute frontal area fractions
if( filedata and Fafb ):
  AE, AN = frontalAreas( Rt )
  print('\n Frontal area fractions of {}:\n Ae/Atot = {}, An/Atot = {}\n'\
    .format(filedata, AE/Atot, AN/Atot)) 


# Create an empty mask id list
mskList = list()
Rxm = None
if( filedata and maskAbove ):
  Rxm, mskList = maskFromClip( Rt, maskAbove, mskList )
elif( maskAbove ):
  Rxm, mskList = maskFromClip( Rm, maskAbove, mskList )
elif( filemask ):
  mskList = maskFromData( Rm, mskList, maxMaskNo )
  Rxm = Rm.copy()
else:
  sys.exit(' Nothing to do. Exiting ...')


# Create an empty mask id list
if( Rxm is not None ):
  ratios = planAreaFractions( Rxm , mskList )
  if( Rt is not None ):
    vm = maskMeanValues( Rxm, Rt, mskList )
  
  if( printOn ):
    Rdims = np.array(Rxm.shape)
    figDims = 13.*(Rdims[::-1].astype(float)/np.max(Rdims))
    pfig = plt.figure(num=1, figsize=figDims)
    pfig = addImagePlot( pfig, Rxm, ' Mask ', gridOn=True )
    plt.show()
