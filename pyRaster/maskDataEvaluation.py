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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

def maskFromClip( Ri, mclip, mlist ):
  Rm  = np.zeros( Ri.shape )
  idx = (Ri > mclip )
  Rm[idx] = 1
  mlist.append(1)

  return Rm, mlist

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

def maskFromData(Ri, mlist, mxN=20):
  for im in range(1,mxN+1):
    if( im in Ri ):
      #print(' Mask id = {} found in file.'.format(im))
      mlist.append(im)

  return mlist

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

def histogram( Rm, Ri, mlist, threshold = 0. ):

  maxv  = int(np.ceil(np.max(Ri)))
  zbins = np.zeros( (len(mlist), maxv) )

  j = 0
  for im in mlist:
    LR, labelCount = labelRaster(Rm, im)
    #Np = float(np.count_nonzero( (LR > 0) ))#; print('Np={}'.format(Np))
    #Np  = float(labelCount)
    Np = 0
    Av = np.zeros( labelCount )
    for l in range(1,labelCount+1):
      idx = (LR==l)
      for v in Ri[idx]:
      #v  = int( np.floor(np.percentile( Ri[idx], 96 )) )
        nw = 1 # = np.count_nonzero(idx) #; print('w={}'.format(w))
        Np += nw
        Av[l-1] = float(v)
        zbins[j,int(v)] += nw
    zbins[j,:] /= Np
    
    idx = ( Av > threshold )
    print(' Mask {}: mean = {}, var = {}, std = {}'\
      .format(im, np.mean(Av[idx]), np.var(Av[idx]), np.std(Av[idx])))
    j += 1

  return zbins

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
parser.add_argument("-ex", "--exBuffer", action="store_true", default=False,\
  help="Consider the effective area excluding frontal buffer.")
parser.add_argument("-p", "--printOn", help="Print the numpy array data.",\
  action="store_true", default=False)
parser.add_argument("--lims", help="User specified colormap and limits for plot.",\
  action="store_true", default=False)
parser.add_argument("--save", help="Save the figure.", action="store_true", default=False)
args = parser.parse_args()
#==========================================================#
# Renaming, nothing more.
filemask  = args.filemask
filedata  = args.filedata
maxMaskNo = args.maxMaskNo
maskAbove = args.maskAbove
Fafb      = args.Fafb # frontal area fraction
exBuffer  = args.exBuffer
printOn   = args.printOn
limsOn    = args.lims
saveFig   = args.save
#==========================================================#

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
  
#----------------------------------------------------------#

ix = 0
if( exBuffer ):
  if( filemask ):
    emx  = int(np.median( Rm[:,0] ))
    print(' Excluding buffer of Mask {} ...'.format(emx))
    Ixm = ( Rm != emx ).astype(int)
    for i in range( len(Rm[0,:]) ):
      nz  = np.count_nonzero( Ixm[:,i] )
      
      if( nz > 3 ): 
        ix = i; break
    
    Rm = Rm[:,ix:].copy(); Rmdims = np.array(np.shape(Rm))
    print(' ix = {}, Rm.shape = {}'.format(ix, Rm.shape))
    
  if( filedata ):
    if( ix == 0 ):
      for i in range( len(Rt[0,:]) ):
        nz = np.count_nonzero( Rt[:,i] )
        if( nz > 0 ): 
          ix = i; break
    Rt = Rt[:,ix:].copy(); Rtdims = np.array(np.shape(Rt))
    print(' ix = {}, Rt.shape = {}'.format(ix, Rt.shape))

#----------------------------------------------------------#

if( filemask ):
  Atot = totalArea( Rmdims, dPxm )
  print('\n Total area of {} domain:\n Atot = {:.4g} m^2 \n'.format(filemask,Atot))
else:
  Atot = totalArea( Rtdims, dPxt )
  print('\n Total area of {} domain:\n Atot = {:.4g} m^2 \n'.format(filedata,Atot))

#----------------------------------------------------------#

# Compute frontal area fractions
if( filedata and Fafb ):
  AE, AN = frontalAreas( Rt )
  print('\n Frontal area fractions of {}:\n Ae/Atot = {:.2f}, An/Atot = {:.2f}\n'\
    .format(filedata, AE/Atot, AN/Atot))

#----------------------------------------------------------#

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

if( filedata ):
  zbins = histogram( Rxm, Rt, mskList )
  np.savetxt('mask_height_histogram.dat', \
    np.c_[np.arange(1,len(zbins[0,:])+1), 100.*np.transpose(zbins) ])
  #print(' sum = {} '.format( np.sum(zbins[0,:])))


# Create an empty mask id list
if( Rxm is not None ):
  ratios = planAreaFractions( Rxm , mskList )
  if( Rt is not None ):
    vmeam, vvar, vstd = maskMeanValues( Rxm, Rt, mskList )

  if( printOn ):
    Rdims = np.array(Rxm.shape)
    figDims = 13.*(Rdims[::-1].astype(float)/np.max(Rdims))
    pfig = plt.figure(num=1, figsize=figDims)
    gridOn = True
    pfig = addImagePlot( pfig, Rxm, ' Mask ', gridOn, limsOn )
    if( saveFig ):
      if( filemask ): figname = filemask.split('.')[0]+'.jpg'
      else:           figname = filedata.split('.')[0]+'.jpg'
      pfig.savefig( figname , format='jpg', dpi=300)
      
    plt.show()
