#!/usr/bin/env python3
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from mapTools import readNumpyZTile, initRdict
from footprintTools import readNumpyZFootprint
from utilities import filesFromList
from plotTools import addImagePlotDict, userLabels, addNests
''' 
Description:


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''

hStr = '''Nests line color. Only effective together with --drawNests. 
Examples: w, r, k, c, g, m. 
Default = k (i.e. black )
'''
#==========================================================#
parser = argparse.ArgumentParser(prog='plotRasterData.py')
parser.add_argument("-f","--filename", type=str, default=None,\
  help="Name of the raster file.")
parser.add_argument("-s", "--size", type=float, default=13.,\
  help="Size of the figure (length of the longer side). Default=13.")
parser.add_argument("-ib", "--ibounds", nargs=2 , type=int, default=[None,None],\
  help="Index bounds in x-direction (easting) for the raster. By default no bounds imposed.")
parser.add_argument("-jb", "--jbounds", nargs=2 , type=int, default=[None,None],\
  help="Index bounds in y-direction (northing) for the raster. By default no bounds imposed.")
parser.add_argument("--abs", action="store_true", default=False,\
  help="Plot absolute values.")
parser.add_argument("--lims", action="store_true", default=False,\
  help="User specified limits.")
parser.add_argument("--grid", help="Turn on grid.", action="store_true", default=False)
parser.add_argument("--cmap", type=str, default=None, \
  help="Matplotlib colormap. Default: Matplotlib default.")
parser.add_argument("--title", type=str, default=None, \
  help="Plot title. By default, the raster filename is used.")
parser.add_argument("--xlabel", type=str, default=None, \
  help="Label for the horizontal axis.")
parser.add_argument("--ylabel", type=str, default=None, \
  help="Label for the vertical axis.")
parser.add_argument("-c", "--coords", choices=["local", "geo", "pixel"], default=None,\
  help="Use raster origin (and resolution) to show either local, geo or pixel coordinates")
parser.add_argument("--labels", action="store_true", default=False,\
  help="User specified labels.")
parser.add_argument("--footprint", action="store_true", default=False,\
  help="Plot footprint data.")
parser.add_argument("-n","--drawNests", action="store_true", default=False,\
  help="Draw nests (rectangles) on top of the raster.")
parser.add_argument("-lc","--nestlinecolor", type=str, default='k', help=hStr)
parser.add_argument("-i","--infoOnly", action="store_true", default=False,\
  help="Print only info to the screen.")
parser.add_argument("--save", metavar="FORMAT" ,type=str, default='', \
  help="Save the figure in specified format. Formats available: jpg, png, pdf, ps, eps and svg")
parser.add_argument("--dpi", metavar="DPI" ,type=int, default=100,\
  help="Desired resolution in DPI for the output image. Default: 100")
args = parser.parse_args() 

#parser.add_argument("--origin", choices=["upper", "lower"], default="upper",\
#  help="Location of the left origin ('upper' or 'lower'). Default: upper")
#writeLog( parser, args )
#==========================================================#

# Renaming ... that's all.
rasterfile  = args.filename
size        = args.size
ib          = args.ibounds
jb          = args.jbounds
absOn       = args.abs
limsOn      = args.lims
gridOn      = args.grid
cmapOn      = args.cmap
title       = args.title
xlabel      = args.xlabel
ylabel      = args.ylabel
coords      = args.coords
infoOnly    = args.infoOnly
labels      = args.labels
drawNests   = args.drawNests
nlcolor     = args.nestlinecolor
footprintOn = args.footprint
save        = args.save
origin      = None

plt.rc('xtick', labelsize=14); #plt.rc('ytick.major', size=10)
plt.rc('ytick', labelsize=14); #plt.rc('ytick.minor', size=6)
plt.rc('axes', titlesize=18)

if( not footprintOn ):
  Rdict = readNumpyZTile(rasterfile)
  Rdict = initRdict( Rdict )
  R = Rdict['R']
  Rdims = np.array(np.shape(R))
  ROrig = Rdict['GlobOrig']
  dPx = Rdict['dPx']
  gridRot = Rdict['gridRot']
  ROrigBL = None
  if( 'GlobOrigBL' in Rdict ):
    ROrigBL = Rdict['GlobOrigBL']
  Rdict = None
else:
  R, X, Y, Z, C = readNumpyZFootprint(rasterfile)
  Rdims = np.array(np.shape(R))
  ROrig = np.zeros(2)
  dPx   = np.array([ (Y[1,0]-Y[0,0]) , (X[0,1]-X[0,0]) ])  # dN, dE
  X = None; Y = None; Z = None; C = None  # Clear memory
  
if( absOn ): R = np.abs(R)

info = ''' Info (Orig):
 Dimensions    [rows, cols] = {0}
 Origin (top-left)    [N,E] = {1}
 Origin (bottom-left) [N,E] = {2}
 Resolution         [dN,dE] = {3}
 Grid rotation (deg)        = {4} deg
 Max(R) / Min(R)            = {5} / {6}
'''.format(Rdims,ROrig,ROrigBL,dPx,gridRot*(180./np.pi),np.nanmax(R),np.nanmin(R))

print(info)

nrows, ncols = R.shape

i0, i1 = 0, ncols
j0, j1 = 0, nrows
if( ib[0] is not None ): i0 = max(0, ib[0])
if( ib[1] is not None ): i1 = min(ncols, ib[1])

if( jb[0] is not None ): j0 = max(0, jb[0])
if( jb[1] is not None ): j1 = min(nrows, jb[1])

if( i0 >= i1 ):
  raise ValueError("Invalid x-index bounds: [{}, {}]".format(i0, i1))

if( j0 >= j1 ):
  raise ValueError("Invalid y-index bounds: [{}, {}]".format(j0, j1))

if( i0 != 0 or i1 != ncols or j0 != 0 or j1 != nrows):
  R = R[j0:j1, i0:i1]
  Rdims = np.array(R.shape)
  print('\n Plot dimensions [rows, cols] = {}'.format(Rdims))

extent = None

if( coords is not None ):
  if( footprintOn) :
    raise ValueError("--coords is currently supported only for raster tile files.")

  # ROrig is [northing, easting] at the top-left raster corner.
  # dPx is [dN, dE].
  if( coords == 'geo' ):
    top    = ROrig[0] - j0 * dPx[0]
    bottom = ROrig[0] - j1 * dPx[0]
    left   = ROrig[1] + i0 * dPx[1]
    right  = ROrig[1] + i1 * dPx[1]
  elif( coords == 'local'):
    top    = j1 * dPx[0]
    bottom = j0 * dPx[0]
    left   = i0 * dPx[1]
    right  = i1 * dPx[1]
  else:
    top    = j1
    bottom = j0
    left   = i0
    right  = i1

  # Matplotlib extent order:
  # [left, right, bottom, top]
  extent = [left, right, bottom, top]

  if( xlabel is None ):
    if(   coords == 'geo'):   xlabel = "Easting"
    elif( coords == 'local'): xlabel = "x-coord. (m)"
    else:                     xlabel = "i coord."
  
  if( ylabel is None ):
    if(   coords == 'geo'):   ylabel = "Northing"
    elif( coords == 'local'): ylabel = "y-coord. (m)"
    else:                     ylabel = "j coord."

if( title is None ): title = rasterfile

plotDict = {
    'R': R,
    'extent': extent,
    'title' : title,
    'xlabel': xlabel,
    'ylabel': ylabel,
    'gridOn': gridOn,
    'limsOn': limsOn,
    'cmap'  : cmapOn,
    'origin': origin
}

if( not infoOnly ):
  
  figDims = size*(Rdims[::-1].astype(float)/np.max(Rdims))
  fig = plt.figure(num=1, figsize=figDims)
  if( drawNests): fig = addNests(fig, nlcolor)
  
  fig = addImagePlotDict(fig, plotDict)
  
  R = None

  if(labels):
    fig = userLabels( fig )

  if(not(save=='')):
    filename = rasterfile.split('/')[-1]  # Remove the path in Linux system
    filename = filename.split('\\')[-1]   # Remove the path in Windows system
    filename = filename.strip('.npz')+'.'+save
    fig.savefig( filename, format=save, dpi=args.dpi)
  
  plt.show()
else:
  R = None

