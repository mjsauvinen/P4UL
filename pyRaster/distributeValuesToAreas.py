#!/usr/bin/env python3
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
hdmsg = '''Use a statistical distribution function to vary values 
mfix to spearate areas. Types available: gaussian and uniform. 
For Gaussian distribution the scale value is standard deviation and 
for the uniform distribution it is the maximum offset. 
Example: gaussian 4.5 .
'''
hmmsg = "Mean of the distribution or constant value if not using distributed values."


parser = argparse.ArgumentParser(prog='distributeValuesToAreas.py', description='''Labels areas from raster data 
and generates random values (e.g. temperatures) matching 
given probability density function and mean.''')
parser.add_argument("-fm","--filemask", type=str,\
  help="Name of the raster data file.")
parser.add_argument("-fa", "--fileaugment", type=str, default=None,\
  help="An existing raster data file to augment. Default=None.")
parser.add_argument("-fo", "--fileout", type=str,\
  help="Name of the output raster data file.")
parser.add_argument("-ft", "--filetopo", type=str, default=None,\
  help="File containing the topography data for VTK results (npz format).")
parser.add_argument("-d", "--distribution", type=str, nargs=2, \
  metavar=('TYPE', 'SCALE'), default=[None,None], help=hdmsg)
parser.add_argument("-m", "--mean", type=float,   help=hmmsg)
parser.add_argument("-id", "--maskIds", type=int,nargs='+',\
  help="List of mask numbers to be included in the value distribution.")
parser.add_argument("-n", "--name", default="Temp", type=str,\
  help="Name of the VTK data array. Default='Temp'.")
parser.add_argument("-v", "--vtk", action="store_true", default=False,\
  help="Write the results in VTK file with topography.")
parser.add_argument("-ix", "--imfix", type=int, nargs='+', default=[None],\
  help="Specific mask ids where fixed prescribed values are used. (Optional).")
parser.add_argument("-vx", "--vmfix", type=float, nargs='+', default=[None],\
  help="Prescribed values for the masks given by --imfix. (Optional).")
parser.add_argument("-vz", "--valueForZeros", type=float, default=None, \
  help=" Value to replace remaining zeros. Mean value --mean used as default.")
parser.add_argument("-p", "--printOn", action="store_true", default=False,\
  help="Print the resulting raster data.")
parser.add_argument("-pp", "--printOnly",action="store_true", default=False,\
  help="Print resulting data without saving.")
args = parser.parse_args()
writeLog(parser, args, args.printOnly)

#==========================================================#
filemask     = args.filemask
fileout      = args.fileout
filetopo     = args.filetopo
fileaug      = args.fileaugment
distribution = args.distribution
maskIds      = args.maskIds
vmean        = args.mean
vtkOn        = args.vtk 
printOnly    = args.printOnly
printOn      = args.printOn
varname      = args.name
imfix        = args.imfix
vmfix        = args.vmfix
vfz          = args.valueForZeros
mfixOn       = False
#==========================================================#

if( (None not in imfix) and (None not in vmfix) ):
  if( len(imfix) != len(vmfix) ):
    sys.exit('ERROR: len(imfix) != len(vmfix). Exiting ...')
  mfixOn = True

# - - - - - - - - - - - - - - - - - - - - - - - - - - #
if( vfz is None ):
  vfz = vmean

# - - - - - - - - - - - - - - - - - - - - - - - - - - #
if(vtkOn and (filetopo is None)):
  sys.exit(' Error: VTK results require -ft/--filetopo. Exiting ...')

# - - - - - - - - - - - - - - - - - - - - - - - - - - #

# Read data into an ndarray
Rdict = readNumpyZTile(filemask)
R = Rdict['R']
Rdims = np.array(np.shape(R))
ROrig = Rdict['GlobOrig']
dPx = Rdict['dPx']
Rdict['R'] = None 

# - - - - - - - - - - - - - - - - - - - - - - - - - - #

# Label shapes from 0 to Nshapes-1 with SciPy ndimage package
if (not(distribution == None)):
  LR, Nshapes = labelRaster(R, maskIds)
else:  # no need for labeling
  LR = R
  Nshapes = 1

# - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Initialize a new array or read existing data
if (fileaug is None):
  Rx = np.zeros(Rdims)
  Rxdims = Rdims
else:
  Rxdict = readNumpyZTile(fileaug)
  Rx = Rxdict['R']
  Rxdict = None
  Rxdims = np.array(np.shape(Rx))
  if (all(Rdims != Rxdims)):
    sys.exit(' Error: size mismatch between two data files when appending.')

# - - - - - - - - - - - - - - - - - - - - - - - - - - #

# Fill the areas with generated values
if(None in distribution):  # Fill with a constant value
  Rx[(LR>0)] = vmean
elif (distribution[0] == "gaussian"):
  for i in range(Nshapes):
    rv = np.random.normal(vmean, float(distribution[1]))
    Rx[(LR ==(i + 1))] = np.random.normal(vmean, float(distribution[1]))
elif (distribution[0] == "uniform"):
  for i in range(Nshapes):
    rv = np.random.uniform(-1.*float(distribution[1]), float(distribution[1])) 
    Rx[(LR ==(i + 1))] = vmean + rv
else:
  sys.exit('Error: invalid distribution given.')

LR = None

# - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Gather the indecies for a mfix mask with fixed value.
idmm = None 
if( mfixOn ):
  for i in range( len(imfix) ):
    idxm = (R==imfix[i])
    if( idmm is None):
      idmm = idxm.copy()
    else:
      idmm = np.maximum( idmm, idxm )
    Rx[idxm] = vmfix[i]
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - #

# Clear memory
R = None

# - - - - - - - - - - - - - - - - - - - - - - - - - - #

# Calculate mean and move nonzero values accordingly
if( idmm is None ):
  idmm = (Rx>0.)
  
offset = np.nanmean(Rx[idmm]) - vmean
print(' Offset to ensure prescribed mean values. V_offset = {}'.format(offset))
Rx[idmm] -= offset
print(' V_max, V_min = {}, {}'.format(np.max(Rx[idmm]), np.min(Rx[idmm])))
idmm = None

# - - - - - - - - - - - - - - - - - - - - - - - - - - #

# Replace remaining zero values. These do not influence the area average
# imposed above. 
Rx[~(Rx>0.)] = vfz 


# - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Read topography data
if ( vtkOn and not printOnly):
  topoDict = readNumpyZTile(filetopo)
  topo = topoDict['R']
  topoDims = np.array(np.shape(topo))
  topoOrig = topoDict['GlobOrig']
  topoDPX = topoDict['dPx']
  topoDict = None

  if(all(topoDims != Rxdims)):
    msg = ''' Error: mismatch in raster data and topography data shapes
    Topo_dims={} vs. Data_dims={}'''.format(topoDims, Rxdims)
    sys.exit(msg)

  # Fill in the coordinate grid
  xa = np.arange( Rdims[1] ).astype(float) * dPx[1]
  ya = np.arange( Rdims[0] ).astype(float) * dPx[0]
  X, Y = np.meshgrid(xa,ya)
  
  # Write the data into a VTK file
  # N axis of (N,E) coordinates has to be reversed
  fname = fileout.split('.')[0]+'.vtk'
  t_vtk = vtkWriteHeaderAndGridStructured2d(X, Y, topo[::-1, :], fname, 'VTK map');
  t_vtk = vtkWritePointDataHeader(t_vtk, Rx[::-1, :], 1)
  t_vtk = vtkWritePointDataStructured2D(t_vtk, Rx[::-1, :], X, varname)

  t_vtk.close()

# - - - - - - - - - - - - - - - - - - - - - - - - - - #

# Save as npz
if(not printOnly):
  Rdict['R'] = Rx; Rdict['dPx'] = dPx; Rdict['GlobOrig'] = ROrig
  Rdict['Nshapes'] = Nshapes
  saveTileAsNumpyZ(fileout, Rdict)
  Rdict = None

# - - - - - - - - - - - - - - - - - - - - - - - - - - #

# Plot the resulting raster
if( printOn or printOnly):
  Rx[Rx == 0] = np.nan  # Replacing zeros with NaN helps plotting
  figDims = 13. * (Rxdims[::-1].astype(float) / np.max(Rxdims))
  fig = plt.figure(num=1, figsize=figDims)
  fig = addImagePlot(fig, Rx, fileout, False, False)
  plt.show()
