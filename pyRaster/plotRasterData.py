#!/usr/bin/env python
import sys
import argparse
import numpy as np
from mapTools import *
from footprintTools import readNumpyZFootprint
from utilities import filesFromList
from plotTools import addImagePlot, userLabels
import matplotlib.pyplot as plt
''' 
Description:


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='plotRasterData.py')
parser.add_argument("rfile", type=str, nargs='?', default=None,\
  help="Name of the comp domain data file.")
parser.add_argument("-s", "--size", type=float, default=13.,\
  help="Size of the figure (length of the longer side). Default=13.")
parser.add_argument("--lims", help="User specified limits.", action="store_true",\
  default=False)
parser.add_argument("--grid", help="Turn on grid.", action="store_true",\
  default=False)
parser.add_argument("--labels", help="User specified labels.", action="store_true",\
  default=False)
parser.add_argument("--footprint", help="Plot footprint data.", action="store_true",\
  default=False)
parser.add_argument("--save", help="Save the figure right away.", action="store_true",\
  default=False)
args = parser.parse_args() 
#writeLog( parser, args )
#==========================================================#
# Renaming ... that's all.
rasterfile  = args.rfile
size        = args.size
limsOn      = args.lims
gridOn      = args.grid
labels      = args.labels
footprintOn = args.footprint
saveOn      = args.save

plt.rc('xtick', labelsize=18); #plt.rc('ytick.major', size=10)
plt.rc('ytick', labelsize=18); #plt.rc('ytick.minor', size=6)
plt.rc('axes', titlesize=20)

if( not footprintOn ):
  Rdict = readNumpyZTile(rasterfile)
  R = Rdict['R']
  Rdims = np.array(np.shape(R))
  ROrig = Rdict['GlobOrig']
  dPx = Rdict['dPx']
  Rdict = None
else:
  R, X, Y, Z, C = readNumpyZFootprint(rasterfile)
  Rdims = np.array(np.shape(R))
  ROrig = np.zeros(2)
  dPx   = np.array([ (X[0,1]-X[0,0]) , (Y[1,0]-Y[0,0]) ])
  X = None; Y = None; Z = None; C = None  # Clear memory
  

info = ''' Info:
 Dimensions [rows, cols] = {}
 Origin (top-left) [N,E] = {}
 Resolution        [N,E] = {}
'''.format(Rdims,ROrig,dPx)

print(info)

figDims = size*(Rdims[::-1].astype(float)/np.max(Rdims))
fig = plt.figure(num=1, figsize=figDims)
fig = addImagePlot( fig, R , rasterfile, gridOn, limsOn)
R = None

if(labels):
  fig = userLabels( fig )

if(saveOn):
  filename = rasterfile.split('/')[-1]  # Remove the path in Linux system
  filename = filename.split('\\')[-1]   # Remove the path in Windows system
  filename = filename.strip('.npz')+'.jpg'
  fig.savefig( filename )
  
plt.show()

