#!/usr/bin/env python3

import argparse
import numpy as np
import sys
from utilities import writeLog
from mapTools import readNumpyZTile, saveTileAsNumpyZ, initRdict
#================================================================================================#
'''
A script for flattening 3d numpy array into 2d raster.

Author: Mikko Auvinen, mikko.auvinen@fmi.fi, FMI

'''
#================================================================================================#
parser = argparse.ArgumentParser(prog='flatten3dNpz.py',\
  description="Flatten 3d numpy array into 2d.")
parser.add_argument("-f", "--filename", type=str, \
  help="Input file (.npz) containing 3d array as 'S' within a dict.")
parser.add_argument("-fo", "--fileout", type=str, default='output.npz',\
  help="Output file which will contain a dict with a 2d array as 'R'.")
parser.add_argument("-mx", "--maxVal", type=float, default=None,\
  help="Force to maintain given maximum value.")
parser.add_argument("-mn", "--minVal", type=float, default=None,\
  help="Force to maintain given minimum value.")
parser.add_argument("-op", "--operator", type=str, choices=['sum','max','mean','mask'], required=True,\
  help="Operator for the flattening.")
parser.add_argument("-ax", "--axis", type=int, choices=[0,1,2], default=2,\
  help="Axis number to flatten. Default is 2 (i.e. the third axis).")
args = parser.parse_args()
writeLog( parser, args )
#=================================================================================================#

filename = args.filename
fileout  = args.fileout.split('.npz')[0]+'.npz'
operator = args.operator
maxVal   = args.maxVal
minVal   = args.minVal
ax       = args.axis

# = = = = = = = = = = = = = = = = = #

# Load files
dat=readNumpyZTile(filename)
Sdict = dict(dat); dat = None
SdictKeys = Sdict.keys()

if( 'dPx' not in SdictKeys ):
  sys.exit(" Resolution information (key:'dPx') is missing in {}. Exiting.".format(filename)) 

if( 'S' in SdictKeys ):
  S=Sdict['S']
else:
  sys.exit(" 'S' not in dict keys: {}".format(SdictKeys))


if(   operator == 'sum' ):
  oper = np.nansum
elif( operator == 'max' ):
  oper = np.nanmax
elif( operator == 'mean' ):
  oper = np.nanmean
else:
  oper = np.nanmax 

# Perform the operation
R = oper( S, axis=ax )

if( maxVal is not None ):
  R = np.minimum( maxVal , R )

if( minVal is not None ):
  R = np.maximum( minVal, R )

if( operator == 'mask' ):
  # Replace nans with zeros temporarily
  idnn = np.isnan(R) 
  R2 = R.copy()
  R2[idnn] = 0.; idnn = None
  idm = (np.abs(R2) > 0.)
  R2 = None
  
  R = R.astype(int)
  R[idm] = 1


Rdict = dict()
dPx   = Sdict['dPx'][:-1]
if('GlobOrig' in SdictKeys): Rdict['GlobOrig'] = Sdict['GlobOrig'][:-1]

Rdict = initRdict(Rdict, R, dPx)

#Save the npz
saveTileAsNumpyZ(fileout, Rdict)
