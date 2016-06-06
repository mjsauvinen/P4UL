#!/usr/bin/env python
import sys
import argparse
import numpy as np
from mapTools import *
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
args = parser.parse_args() 
#writeLog( parser, args )
#==========================================================#

# Renaming ... that's all.
rasterfile = args.rfile
size       = args.size
limsOn     = args.lims
gridOn     = args.grid
labels     = args.labels

R, Rdims, ROrig, dPx = readNumpyZTile(rasterfile)

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

plt.show()

