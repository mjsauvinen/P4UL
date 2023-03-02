#!/usr/bin/env python3
import sys
import argparse
import numpy as np
from mapTools import readNumpyZTile, initRdict
from footprintTools import readNumpyZFootprint
from utilities import filesFromList
from plotTools import addImagePlot, userLabels, addNests
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
parser.add_argument("--labels", action="store_true", default=False,\
  help="User specified labels.")
parser.add_argument("--footprint", action="store_true", default=False,\
  help="Plot footprint data.")
parser.add_argument("-n","--drawNests", action="store_true", default=False,\
  help="Draw nests (rectangles) on top of the raster.")
parser.add_argument("-lc","--nestlinecolor", type=str, choices=['w','r','k','c','g','m'], default='k',\
  help="Nests line color. Only effective together with --drawNests. Default='k' i.e. black ")
parser.add_argument("-i","--infoOnly", action="store_true", default=False,\
  help="Print only info to the screen.")
parser.add_argument("--save", metavar="FORMAT" ,type=str, default='', \
  help="Save the figure in specified format. Formats available: jpg, png, pdf, ps, eps and svg")
parser.add_argument("--dpi", metavar="DPI" ,type=int, default=100,\
  help="Desired resolution in DPI for the output image. Default: 100")
args = parser.parse_args() 
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
infoOnly    = args.infoOnly
labels      = args.labels
drawNests   = args.drawNests
nlcolor     = args.nestlinecolor
footprintOn = args.footprint
save        = args.save

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

imod = False; jmod = False
if( np.count_nonzero( jb ) > 0 ):
  jb[0] = max(         0 , jb[0] )
  jb[1] = min( Rdims[0]-1, jb[1] )
  jmod = True

if( np.count_nonzero( ib ) > 0 ):
  ib[0] = max( 0         , ib[0] )
  ib[1] = min( Rdims[1]-1, ib[1] )
  imod = True

if( imod or jmod ):
  R = R[jb[0]:jb[1],ib[0]:ib[1]]
  Rdims = np.array( R.shape )
  print('\n Plot dimensions [rows, cols] = {}'.format(Rdims))



if( not infoOnly ):
  figDims = size*(Rdims[::-1].astype(float)/np.max(Rdims))
  fig = plt.figure(num=1, figsize=figDims)
  if( drawNests): fig = addNests(fig, nlcolor)
  fig = addImagePlot( fig, R , rasterfile, gridOn, limsOn)
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

