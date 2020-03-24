#!/usr/bin/env python3
import sys
import argparse
import numpy as np
from mapTools import *
from utilities import filesFromList, writeLog
from plotTools import addImagePlot, addScatterPlot
import matplotlib.pyplot as plt
'''
Description:


Author: Mikko Auvinen
        mikko.auvinen@fmi.fi 
        Finnish Meteorological Institute
'''
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
def addBlocks( T, stride, Lx, h ):
  Tdims = T.shape
  sy = stride[0]; sx = stride[1]
  ly = Lx[0];     lx = Lx[1]
  
  for i in range( int(np.ceil(Tdims[1]/sx)+1) ):
    ix  = i*sx 
    ix2 = ix + lx
    for j in range( int( np.ceil(Tdims[0]/sy)+1) ):
      jy = j*sy + int(np.mod(i,2)*(sy/2))
      jy2 = jy+ly
      if( ix2 > Tdims[1] or jy2 > Tdims[0] ):
        break
      else:
        #print(' ix1: ix2 = {}:{}, jy1:jy2 = {}:{} '.format(ix,ix2,jy,jy2))
        T[jy:jy2, ix:ix2] += h

  return T
 # =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 

#==========================================================#
parser = argparse.ArgumentParser(prog='addBlockMargin.py')
parser.add_argument("-f", "--filename",type=str, help="Name of the comp domain data file.")
parser.add_argument("-fo", "--fileout",type=str, help="Name of output Palm topography file.")
parser.add_argument("-s","--stride", help="Stride lengths for the block arrangement. [N, E]",\
  type=int,nargs=2,default=[None,None])
parser.add_argument("-L","--Lblocks", help="Block dimensions. [W, L]",\
  type=int,nargs=2,default=[None,None])
parser.add_argument("-mw","--mrgnW", help="Zero or non-zero margin widths as ratios (0-1): [L,R,B,T]",\
  type=float,nargs=4,default=[None,None,None,None])
parser.add_argument("-mh","--mrgnH", help="Margins block heights: [L,R,B,T]. Default=0",\
  type=float,nargs=4,default=[0.,0.,0.,0.])
parser.add_argument("-wa", "--writeAscii", help="Write 'TOPOGRAPHY_DATA' ascii file.",\
  action="store_true", default=False)
parser.add_argument("-z", "--zero", help="Zero the raster file first.",\
  action="store_true", default=False)
parser.add_argument("-p", "--printOn", help="Print the resulting raster data.",\
  action="store_true", default=False)
parser.add_argument("-pp", "--printOnly", help="Only print the resulting data. Don't save.",\
  action="store_true", default=False)
args = parser.parse_args()
writeLog( parser, args, args.printOnly )
#==========================================================#

filename   = args.filename
fileout    = args.fileout
mw         = args.mrgnW
mh         = args.mrgnH
stride     = args.stride
Lb         = args.Lblocks
zeroAll    = args.zero
printOn    = args.printOn
printOnly  = args.printOnly
writeAscii = args.writeAscii

if( mw.count(None) != 0 ):
  sys.exit(' Error! One of the margins widths is None. Exiting ...')

if( stride.count(None) != 0 ):
  sys.exit(' Error! One of the stride lengths is None. Exiting ...')

if( Lb.count(None) != 0 ):
  sys.exit(' Error! One of the block dimensions is None. Exiting ...')


# Read the raster tile to be processed.
Rdict = readNumpyZTile(filename)
R = Rdict['R']
Rdims = np.array(np.shape(R))
ROrig = Rdict['GlobOrig']
print(' Rdims = {} '.format(Rdims))
print(' ROrig = {} '.format(ROrig))

if( zeroAll ):
  R[:,:] = 0.


L12, R12, B12, T12 = marginIds( Rdims, mw )
L1 = L12[0]; L2 = L12[1]
R1 = R12[0]; R2 = R12[1]
B1 = B12[0]; B2 = B12[1]
T1 = T12[0]; T2 = T12[1]


if( not all( L12 == 0 ) ): R[:,L1:L2] = addBlocks( R[:,L1:L2], stride, Lb, mh[0] )
if( not all( R12 == 0 ) ): R[:,R1:R2] = addBlocks( R[:,R1:R2], stride, Lb, mh[1] )
if( not all( T12 == 0 ) ): R[T1:T2,:] = addBlocks( R[T1:T2,:], stride, Lb, mh[2] )
if( not all( B12 == 0 ) ): R[B1:B2,:] = addBlocks( R[B1:B2,:], stride, Lb, mh[3] )

if( not args.printOnly ):
  Rdict['R'] = R
  saveTileAsNumpyZ( fileout, Rdict )
  if( writeAscii ):
    fout= 'TOPOGRAPHY_DATA_BLOCK'
    np.savetxt(fout,np.round(R),fmt='%g')
  

if( args.printOn or args.printOnly ):
  
  figDims = 13.*(Rdims[::-1].astype(float)/np.max(Rdims))
  fig = plt.figure(num=1, figsize=figDims)
  fig = addImagePlot( fig, R, filename )

  plt.show()


