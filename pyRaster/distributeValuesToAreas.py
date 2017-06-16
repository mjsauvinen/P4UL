#!/usr/bin/env python
from utilities import vtkWriteDataStructured2d, vtkWriteHeaderAndGridStructured2d
from utilities import vtkWritePointDataHeader, vtkWritePointDataStructured2D
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
parser = argparse.ArgumentParser(prog='distributeValuesToAreas.py', description='''Labels areas from raster data and generates random values (e.g. temperatures) matching given probability density function and mean.
''')
parser.add_argument("rfile", type=str, nargs='?', default=None,
                    help="Name of the raster data file.")
parser.add_argument("-a", "--add", metavar="DFILE", type=str,
                    help="Add data to an existing raster data file.")
parser.add_argument("-fo", "--fileout", type=str,
                    help="Name of the output raster data file.")
parser.add_argument("-p", "--printOn", help="Print the resulting raster data.",
                    action="store_true", default=False)
parser.add_argument("-pp", "--printOnly", help="Print resulting data without saving.",
                    action="store_true", default=False)
parser.add_argument("-d", "--distribution", type=str, nargs=2, metavar=('TYPE', 'SCALE'),
                    help="Use a statistical distribution function to vary values specific to spearate areas. Types available: gaussian and uniform. For Gaussian distribution the scale value is standard deviation and for the uniform distribution it is the maximum offset. Example: gaussian 4.5 .")
parser.add_argument("-m", "--mean", type=float,
                    help="Mean of the distribution or constant value if not using distributed values.")
parser.add_argument("-n", "--name", default="Temperature", type=str,
                    help="Name of the VTK data array. Leave empty for 'Temperature'.")
parser.add_argument("-v", "--vtk", metavar="VTKFILE", type=str,
                    help="Write the results in VTKFILE with topography.")
parser.add_argument("-ft", "--filetopo", type=str,
                    help="File containing the topography data for VTK results (npz format).", default='')
args = parser.parse_args()
writeLog(parser, args)

#==========================================================#


if(args.vtk and (args.filetopo == '')):
  sys.exit(' Error: VTK results require -ft/--filetopo. Exiting ...')

# Read data into an ndarray
Rdict = readNumpyZTile(args.rfile)
R = Rdict['R']
Rdims = np.array(np.shape(R))
ROrig = Rdict['GlobOrig']
dPx = Rdict['dPx']

# Label shapes from 0 to shapeCount-1 with SciPy ndimage package
if (not(args.distribution == None)):
  LR, shapeCount = labelRaster(R)
else:  # no need for labeling
  LR = R
  shapeCount = 1
R = None

# Initialize a new array or read existing data
if (args.add == None):
  R = np.zeros(Rdims)
else:
  Rdict = readNumpyZTile(args.add)
  R = Rdict['R']
  Rdims2 = np.array(np.shape(R))
  if (all(Rdims != Rdims2)):
    sys.exit(' Error: size mismatch between two data files when appending.')

# Fill the areas with generated values
if (args.distribution == None):  # Fill with a constant value
  R[np.nonzero(LR)] = args.mean
elif (args.distribution[0] == "gaussian"):
  for i in xrange(shapeCount):
    R[LR == i + 1] = np.random.normal(args.mean, args.distribution[1])
elif (args.distribution[0] == "uniform"):
  for i in xrange(shapeCount):
    R[LR == i + 1] = args.mean + \
        (np.random.uniform(-args.distribution[1], args.distribution[1]))
else:
  sys.exit('Error: invalid distribution given.')
LR = None

# Calculate mean and move nonzero values accordingly
offset = np.nanmean(R[np.nonzero(R)]) - args.mean
R[np.nonzero(R)] = R[np.nonzero(R)] - offset

# Read topography data
if (not(args.vtk) == None and not(args.printOnly)):
  topoDict = readNumpyZTile(args.filetopo)
  topo = topoDict['R']
  topoDims = np.array(np.shape(topo))
  topoOrig = topoDict['GlobOrig']
  topoDPX = topoDict['dPx']
  topoDict = None

  if(all(topoDims != Rdims)):
    sys.exit(' Error: mismatch in raster data and topography data shapes, Topo_dims={} vs. Data_dims={}').format(
        topoDims, Rdims)

  # Fill in the coordinate grid
  X = np.zeros(Rdims)
  Y = np.zeros(Rdims)
  for i in xrange(Rdims[0]):
    X[i, :] = i;
  for i in xrange(Rdims[1]):
    Y[:, i] = i

  # Write the data into a VTK file
  # N axis of (N,E) coordinates has to be reversed
  t_vtk = vtkWriteHeaderAndGridStructured2d(
      Y, X, topo[::-1, :], args.vtk, 'VTK map');
  t_vtk = vtkWritePointDataHeader(t_vtk, R[::-1, :], 1)
  t_vtk = vtkWritePointDataStructured2D(t_vtk, R[::-1, :], Y, args.name)

  t_vtk.close()

# Save as npz
if(not args.printOnly):
  Rdict['R'] = R; Rdict['dPx']: dpx; Rdict['GlobOrig']: ROrig; Rdict['ShapeCount']: shapecount
  saveTileAsNumpyZ(args.fileout, Rdict)
  Rdict = None

# Plot the resulting raster
if(args.printOn or args.printOnly):
  R[R == 0] = np.nan  # Replacing zeros with NaN helps plotting
  figDims = 13. * (Rdims[::-1].astype(float) / np.max(Rdims))
  fig = plt.figure(num=1, figsize=figDims)
  fig = addImagePlot(fig, R, args.rfile, False, False)
  plt.show()
