#!/usr/bin/env python3
import sys
import argparse
import numpy as np
from mapTools import *
from utilities import filesFromList, writeLog
from plotTools import addImagePlot, addContourf, addScatterPlot
import matplotlib.pyplot as plt
'''
Description:
Extracts a rotated and/or cropped subdomain from a raster file.

Authors:  Mikko Auvinen
          mikko.auvinen@helsinki.fi
          University of Helsinki &
          Finnish Meteorological Institute

          Sasu Karttunen
          sasu.karttunen@helsinki.fi
          University of Helsinki

          Mona Kurppa
          mona.kurppa@helsinki.fi
          University of Helsinki
'''


#==========================================================#
parser = argparse.ArgumentParser(prog='extractDomainFromTile.py')
parser.add_argument("-f", "--filename",type=str, help="Name of raster data file.")
parser.add_argument("-fo", "--fileout",type=str, help="Name of output raster file.")
parser.add_argument("-iP","--iPivot", help="Local pixel ids [N,E] for the pivot in the raster file.",\
  type=int,nargs=2, required=True)
parser.add_argument("-gc","--useGlobCoords", help="Use global coordinates in iP. ",
                    action="store_true", default=False)
parser.add_argument("-N","--NxG", type=int, nargs=2,\
  help="Number of points [Nx, Ny] in the 2D Palm grid.")
parser.add_argument("-dx","--dxG", type=float, nargs=2, default=[ 2. , 2.],\
  help="Resolution [dx, dy] of the 2D Palm grid.")
parser.add_argument("-r","--rLx", type=float,nargs=2, default=[ 0.5, 0.5],\
  help="Pivot location [rLx, rLy] as ratio of Lx & Ly of grid domain (top left origo).")
parser.add_argument("-wd", "--windDir", type=float,default=270.0,\
  help="Wind direction (deg) --> Rotation angle around the pivot point. North wind = 0deg")
parser.add_argument("-nr", "--noRotation", action="store_true",default=False,\
  help="Do not rotate the grid.")
parser.add_argument("-s", "--scale",type=float,\
  help="Scale factor for the output. Default=1.", default=1.)
parser.add_argument("-p", "--printOn", help="Print the resulting raster data.",\
  action="store_true", default=False)
parser.add_argument("-pp", "--printOnly", help="Only print the resulting data. Don't save.",\
  action="store_true", default=False)
parser.add_argument("-v", "--verbose", help="Print intermediates onto the screen.",\
  action="store_true", default=False)
args = parser.parse_args()
writeLog( parser, args, args.printOnly )
#==========================================================#


# Renaming the argument variables for brevity and clarity:
filename  = args.filename
fileout   = args.fileout
NxG       = args.NxG
iPv       = args.iPivot
dxG       = args.dxG
rLx       = args.rLx
windDir   = args.windDir
noRotation= args.noRotation
verbose   = args.verbose
printOn   = args.printOn
printOnly = args.printOnly

# Read in the underlying topography data and obtain the pivot coordinates.
dataOnly = False
Rdict= readNumpyZTileForMesh( filename )
R = Rdict['R']
rY = Rdict['rowCoords']   # Local (NOT global!) row coords.
cX = Rdict['colCoords']   # Local (NOT global) column coords.


if( verbose ):
  print(' Local: [N] coords = {}...{}, [E] coords = {}...{}'\
    .format(rY[0], rY[-1], cX[0], cX[-1] ))

Rdims = np.array(np.shape(R))
# Retain information about rotation
try:
  gridRot = Rdict['gridRot']
except:
  gridRot = 0.
ROrig = Rdict['GlobOrig']

if (args.useGlobCoords==True):
  iPv[1] = iPv[1] - ROrig[1]
  iPv[0] = ROrig[0]-iPv[0]


dPx = entry2Int( Rdict['dPx'] )
if( verbose ): print(' dPx = {} '.format(dPx))


# Pivot coordinates in the local coord. system
pY = rY[iPv[0]]; pX = cX[iPv[1]]
print(' Origo in the input Topography data: [N,E] = [{}, {}]'.format(ROrig[0],ROrig[1]))
print(' Pivot Coords in the input Topography data: [N,E] = [pY={}, pX={}]'.format(pY,pX))
#NY, EX = np.meshgrid(rY,cX)

'''
Create Palm grid which obeys the X,Y-coordinate layout. This might cause confusion
so let's proceed carefully.

NOTE:
Even though the data is saved as raster array, the data points are now cell centers.
'''
xbegin = 0. # dxG[0]/2. # First cell-centers.
ybegin = 0. # dxG[1]/2.
xend = (NxG[0]-1)*dxG[0] - xbegin   # Last cell-center.
yend = (NxG[1]-1)*dxG[1] - ybegin

XgridCoords = np.linspace(xbegin,xend, NxG[0])
YgridCoords = np.linspace(ybegin,yend, NxG[1])
#Xg, Yg = np.meshgrid( XgridCoords, YgridCoords )

# Location of the pivot (indecies and coords) in the Palm grid. Not going into negative indices.
# The -1 omitted in iPGy due to the computation of Irows (vs. Jcols) to preserve symmetry.
iPGx = np.maximum(int(    rLx[0] *(NxG[0]-1) ), 0)
iPGy = np.maximum(int((1.-rLx[1])*(NxG[1]  ) ), 0)
iPGy = np.minimum( iPGy , (NxG[1]-1) )

pXG = XgridCoords[iPGx]
pYG = YgridCoords[iPGy]

if( verbose ):
  print(' Grid coords X = {} '.format( (XgridCoords[0], XgridCoords[-1]) ))
  print(' Grid coords Y = {} '.format( (YgridCoords[0], YgridCoords[-1]) ))
  print( ' Palm grid pivot indices: iPGx = {}, iPGy = {}'.format( iPGx, iPGy ))
  print( ' Palm grid pivot coords:   pXG = {},  pYG = {}'.format( pXG, pYG ))


'''
From palm coordinates to underlying local topography coordinates.
We use the pivot point which is known for both systems.
'''

dXT = pX - pXG
dYT = pY - pYG
#dXT = (pX - ROrig[1]) - pXG
#dYT = (pY - ROrig[0]) - pYG

XT = XgridCoords + dXT
YT = YgridCoords + dYT

if(verbose):
  print(' Coordinate transform:  dXT = {},  dYT = {}'.format( dXT, dYT ))
  print(' Transformed coords: XT = {}...{}, YT = {}...{}'\
    .format(np.min(XT), np.max(XT),np.min(YT), np.max(YT)))

'''
Rotate the new coordinates (within the local coordinate system) according to the wind direction:
Coordinate transformations for counterclockwise rotation.
'''
# NOTE: At the pivot point XTR = pX
XTM, YTM = np.meshgrid( XT, YT ); XT  = None; YT  = None
if (noRotation):
  theta = 0.
else:
  theta = 270. - windDir

if( theta != 0.):
  XTRM,YTRM = rotateGridAroundPivot(XTM,YTM, pX, pY,theta, deg=True)
else:
  print(' No rotation! ')
  XTRM = XTM.copy(); YTRM = YTM.copy()

'''
Bottom Left  :  XTRM[0,0], YTRM[0,0]
Top Left     :  XTRM[-1,0], YTRM[-1,0]
Bottom Right :  XTRM[0,-1], YTRM[0,-1]
Top Right    :  XTRM[-1,-1], YTRM[-1,-1])
'''

XTM = None; YTM = None

'''
Using the known transformed coordinates, we can extract the pixel values
at those locations and copy them to the Palm grid. The grid arrays origo is
located at the bottom left, which makes things a bit confusing here.
'''

Irow = ((ROrig[0]-YTRM)/dPx ).astype(int)
Jcol = ((XTRM-ROrig[1])/dPx ).astype(int)

if( verbose ):
  print(' Irow = {}'.format( (Irow[0,0], Irow[-1,-1]) ))
  print(' Jcol = {}'.format( (Jcol[0,0], Jcol[-1,-1]) ))


# Make sure the indecies don't run beyond the allowable bounds.
if (np.amin(Irow) < 0 or np.amin(Jcol) < 0):
  # Warn the user about streching edges.
  print("WARNING: Domain out of raster data bounds! Streching edge cells to fill the domain.")

Irow = np.maximum(Irow, 0);          Jcol = np.maximum(Jcol, 0)
Irow = np.minimum(Irow, Rdims[0]-1); Jcol = np.minimum(Jcol, Rdims[1]-1)

#print " np.shape(Irow) = {},  Irow = {} ".format(np.shape(Irow) ,Irow[::4,::4])
#print " Jcol = {} ".format(Jcol[::4,::4] )
Xdims = np.array( np.shape(XTRM) )
PR = np.zeros( Xdims  , float)
PR[::-1,:] = R[Irow,Jcol]    # The row order must be reversed to go back to raster format.
R = None

'''
 Reset the top left origo such that it lies in the global coordinate system.
 This requires that we rotate the new local coordinates XTRM, YTRM to the original
 coord. system using the input raster's top left origin as pivot.
'''
theta2 = gridRot/(np.pi/180.)
XTRM,YTRM = rotateGridAroundPivot(XTRM,YTRM, ROrig[1], ROrig[0],theta2, deg=True)
PROrig = np.array([ YTRM[-1,0], XTRM[-1,0] ])  # Reset top left origo
PROrigBL = np.array([ YTRM[0,0], XTRM[0,0] ])  # Reset top left origo
print(' Top left origo coords. (cell centers!): [N,E] = {}'.format(PROrig))
rotation = (theta+theta2)*(np.pi/180.)
#print((theta+theta2)*(np.pi/180.))

# Retain unused keys from original raster
PRdict = Rdict.copy()
Rdict = None
PRdict['R'] = PR
PRdict['GlobOrig']   = PROrig
PRdict['GlobOrigBL'] = PROrigBL 
PRdict['gridRot'] = rotation
PRdict['dPx'] = np.array([dxG[0],dxG[1]])

if( not args.printOnly ):
  saveTileAsNumpyZ( fileout, PRdict)


# Print the raster map, first, in a coordinate system where x-axis is aligned with the windDir
# and, second, in its original orientation.
if( printOn or printOnly ):
  figDims = 13.*(Xdims[::-1].astype(float)/np.max(Xdims))
  fig = plt.figure(num=1, figsize=figDims)
  fig = addImagePlot( fig, PR, args.fileout )

  CfD = dict()
  CfD['title']=' Z(x,y) '; CfD['label']="PALM DOMAIN ON MAP"; CfD['N']=16
  CO = addContourf( XTRM, YTRM, PR[::-1,:], CfD )
  plt.show()

XTRM = None; YTRM = None
PR   = None; PRDict = None
