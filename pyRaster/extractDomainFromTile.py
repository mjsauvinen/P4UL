#!/usr/bin/env python
import sys
import argparse
import numpy as np
from mapTools import *
from utilities import filesFromList, writeLog
from plotTools import addImagePlot, addContourf, addScatterPlot
import matplotlib.pyplot as plt
''' 
Description:


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''


#==========================================================#
parser = argparse.ArgumentParser(prog='extractDomainFromTile.py')
parser.add_argument("-f", "--filename",type=str, help="Name of raster data file.")
parser.add_argument("-fo", "--fileOut",type=str, help="Name of output Palm mesh file.",\
  default="PalmTopo")
parser.add_argument("-iP","--iPivot", help="Pixel ids [N,E] for the pivot in the raster file.",\
  type=int,nargs=2,required=True)
parser.add_argument("-N","--NxG", help="Number of points [Nx, Ny] in the 2D Palm grid.",\
  type=int,nargs=2, default=[ 2048 , 1024])
parser.add_argument("-dx","--dxG", help="Resolution [dx, dy] of the 2D Palm grid.",\
  type=float,nargs=2, default=[ 2. , 2.])
parser.add_argument("-r","--rLx", type=float,nargs=2, default=[ 0.9, 0.5],\
  help="Pivot location [rLx, rLy] as ratio of Lx & Ly of grid domain.")
parser.add_argument("-wd", "--windDir", type=float,default=0.,\
  help="Wind direction (deg) --> Rotation angle around the pivot point. North wind = 0deg")
parser.add_argument("-s", "--scale",type=float,\
  help="Scale factor for the output. Default=1.", default=1.)
parser.add_argument("-p", "--printOn", help="Print the resulting raster data.",\
  action="store_true", default=False) 
parser.add_argument("-pp", "--printOnly", help="Only print the resulting data. Don't save.",\
  action="store_true", default=False) 
args = parser.parse_args() 
writeLog( parser, args, args.printOnly )
#==========================================================#

# Renaming the argument variables for brevity and clarity:
NxG = args.NxG
iPv = args.iPivot
dxG = args.dxG
rLx = args.rLx
windDir = args.windDir

# Read in the underlying topography data and obtain the pivot coordinates.
dataOnly = False
Rdict= readNumpyZTileForMesh( args.filename )
R = Rdict['R']
nY = Rdict['rowCoords']
eX = Rdict['colCoords']
Rdims = np.array(np.shape(R))
ROrig = Rdict['LocalOrig']
dPx = entry2Int( Rdict['dPx'] )
Rdict = None

# Pivot coordinates
pY = nY[iPv[0]]; pX = eX[iPv[1]] 
print(' Origo in the Topography data: [XO,YO] = [{}, {}]'.format(ROrig[1],ROrig[0]))
print(' Pivot Coords in Topography data: [pX,pY] = [{}, {}]'.format(pX,pY))
#NY, EX = np.meshgrid(nY,eX)

'''
Create Palm grid which obeys the X,Y-coordinate layout. This might cause confusion
so let's proceed carefully.

NOTE: 
Even though the data is saved as raster array, the data points are now cell centers.
'''
xbegin = dxG[0]/2. # First cell-centers.
ybegin = dxG[1]/2. 
xend = NxG[0]*dxG[0] - xbegin   # Last cell-center.
yend = NxG[1]*dxG[1] - ybegin

XgridCoords = np.linspace(xbegin,xend, NxG[0]) 
YgridCoords = np.linspace(ybegin,yend, NxG[1])
#Xg, Yg = np.meshgrid( XgridCoords, YgridCoords )

# Location of the pivot (indecies and coords) in the Palm grid.
iPGx = int(rLx[0]*NxG[0])-1   
iPGy = int(rLx[1]*NxG[1])-1

pXG = XgridCoords[iPGx]
pYG = YgridCoords[iPGy]

#print ' Palm grid: iPGx = {}, iPGy = {}'.format( iPGx, iPGy )
#print ' Palm grid: pXG = {}, pYG = {}'.format( pXG, pYG )

'''
From palm coordinates to underlying topography coordinates.
We use the pivot point which is known for both systems.
'''
dXT = pX - pXG; dYT = pY - pYG
XT = XgridCoords + dXT
YT = YgridCoords + dYT
'''  
Rotate the new coordinates according to the wind direction:
Coordinate transformations for counterclockwise rotation.
'''


# NOTE: At the pivot point XTR = pX 
XTM, YTM = np.meshgrid( XT, YT )
theta = 270. - windDir
if( theta != 0. ):
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


'''
 Reset the top left origo utilizing the NON-rotated coordinates. This 
 allows the relative position of different raster maps (with identical 
 coord. rotation) to be determined easily.
'''
PROrig = np.array([ XTM[-1,0], YTM[-1,0] ])  # Reset top left origo
print(' Top left origo coords. (cell centers!): [X,Y] = {}'.format(PROrig))

XT  = None; YT  = None
XTM = None; YTM = None

''' 
Using the known transformed coordinates, we can extract the pixel values
at those locations and copy them to the Palm grid. The grid arrays origo is 
located at the bottom left, which makes things a bit confusing here.
'''
Irow = (np.abs(YTRM-ROrig[0])/dPx ).astype(int)
Jcol = (np.abs(XTRM-ROrig[1])/dPx ).astype(int) 

# Make sure the indecies don't run beyond the allowable bounds.
Irow = np.maximum(Irow, 0);          Jcol = np.maximum(Jcol, 0)
Irow = np.minimum(Irow, Rdims[0]-1); Jcol = np.minimum(Jcol, Rdims[1]-1)

#print " np.shape(Irow) = {},  Irow = {} ".format(np.shape(Irow) ,Irow[::4,::4])
#print " Jcol = {} ".format(Jcol[::4,::4] )
Xdims = np.array( np.shape(XTRM) )
PR = np.zeros( Xdims  , float)
PR[::-1,:] = R[Irow,Jcol]    # The row order must be reversed. 
R = None
PRdict = {'R' : PR, 'LocalOrig' : PROrig, 'dPx' : np.array([dxG[0],dxG[1]])}

if( not args.printOnly ):
  saveTileAsNumpyZ( args.fileOut, PRdict)


# I'm not fully sure why the row indices have to be fed in reverse order ...
if( args.printOn or args.printOnly ):
  figDims = 13.*(Xdims[::-1].astype(float)/np.max(Xdims))
  fig = plt.figure(num=1, figsize=figDims)
  fig = addImagePlot( fig, PR, args.fileOut )
 
  CO = addContourf( XTRM, YTRM, PR[::-1,:], " Z(X,Y) ", "PALM DOMAIN ON MAP" )
  plt.show()

XTRM = None; YTRM = None
PR   = None; PRDict = None

  