#!/usr/bin/env python
import sys
import argparse
import numpy as np
from utilities import writeLog
from mapTools import * 
from plotTools import addImagePlot
import matplotlib.pyplot as plt

''' 
Description:
Labels areas from raster data and generates random values 
(e.g. temperatures) matching given probability density function and mean.
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='distributeValuesForAreas.py')
parser.add_argument("rfile", type=str, nargs='?', default=None,\
  help="Name of the raster data file.")
parser.add_argument("-fo", "--fileout",type=str, help="Name of the output raster data file.")
parser.add_argument("-p", "--printOn", help="Print the resulting raster data.",\
  action="store_true", default=False)
parser.add_argument("-pp", "--printOnly", help="Print resulting data without saving.",\
  action="store_true", default=False)
parser.add_argument("-d", "--distribution",type=str,nargs=2, metavar=('DIST', 'VALUE'), help="Distribution and its associated value used for varying values. Choices: gaussian and uniform. For Gaussian distribution the value required is standard deviation and for the uniform distribution it is the maximum offset. Example: gaussian 4.5 .")
parser.add_argument("-m", "--mean",type=float, help="Mean of the distribution.")
args = parser.parse_args()
writeLog( parser, args )

#==========================================================#

dist            =    args.distribution[0]
scaleValue      =    float(args.distribution[1])
mean            =    args.mean

# Read data into an ndarray
R, Rdims, ROrig, dPx = readNumpyZTile(args.rfile)
Rdims = np.array(np.shape(R))

# Label shapes from 0 to shapeCount-1 with SciPy ndimage package
LR, shapeCount = labelRaster(R)
R = None
R = np.zeros(Rdims)

# Generate new values for areas
if (dist=="gaussian"):
  for i in xrange(shapeCount):
    R[LR==i+1]=np.random.normal(args.mean,scaleValue)
elif (dist=="uniform"):
  for i in xrange(shapeCount):
    R[LR==i+1]=args.mean + (np.random.uniform(-scaleValue,scaleValue))
  
LR=None
# Calculate mean and move nonzero values accordingly
offset = np.nanmean(R[np.nonzero(R)])-mean
R[np.nonzero(R)] = R[np.nonzero(R)]-offset

# plot the resulting raster
if( args.printOn or args.printOnly ):
  R[R==0] = np.nan # Replacing zeros with NaN helps plotting
  figDims=13.*(Rdims[::-1].astype(float)/np.max(Rdims))
  fig = plt.figure(num=1,figsize=figDims)
  fig = addImagePlot(fig,R,args.rfile,False,False)
  plt.show()
  
# save as npz
if( not args.printOnly ):
  saveTileAsNumpyZ( args.fileout, R, Rdims, ROrig, dPx )