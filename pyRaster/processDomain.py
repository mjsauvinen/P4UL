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
#==========================================================#

parser = argparse.ArgumentParser(prog='processDomain.py')
parser.add_argument("-f", "--filename",type=str, help="Name of the comp domain data file.")
parser.add_argument("-fo", "--fileout",type=str, help="Name of output Palm topography file.")
parser.add_argument("-i0","--iZero", help="Pixel ids [N,E] for the zero level.",\
  type=int,nargs=2,default=[None,None])
parser.add_argument("-na", "--nansAbove", type=float, default=None,\
  help="Replace values above given threshold by <nans> (i.e. fill values). Default=None")
parser.add_argument("-nb", "--nansBelow", type=float, default=None,\
  help="Replace values below given threshold by <nans> (i.e. fill values). Default=None")
parser.add_argument("-mw","--mrgnW", help="Zero or non-zero margin widths as ratios (0-1): [L,R,B,T]",\
  type=float,nargs=4,default=[None,None,None,None])
parser.add_argument("-mr","--mrgnR", help="Margin ramp widths as ratios (0-1): [L,R,B,T]",\
  type=float,nargs=4,default=[None,None,None,None])
parser.add_argument("-mh","--mrgnH", help="Margins heights: [L,R,B,T]. Default=0",\
  type=float,nargs=4,default=[0.,0.,0.,0.])
helpFlt = ''' Filter type and its associated number. Available filters:
 median, percentile, rank, gaussian, local. Entering \"user, num\" allows the user
 to specify <num> different filters consecutively.
 Example entry: median 5'''
parser.add_argument("-ft","--filter",type=str,nargs=2,default=[None,None], help=helpFlt)
parser.add_argument("-rx","--rmax", type=float, default=None,\
  help="Recover peaks (after filtering) above given value.")
parser.add_argument("-hx","--hmax", type=float, default=None,\
  help="Maximum allowable height.")
parser.add_argument("-p", "--printOn", action="store_true", default=False,\
  help="Print the resulting raster data.")
parser.add_argument("-pp", "--printOnly", help="Only print the resulting data. Don't save.",\
  action="store_true", default=False)
args = parser.parse_args()
writeLog( parser, args, args.printOnly )
#==========================================================#

filename= args.filename
fileout = args.fileout
na      = args.nansAbove
nb      = args.nansBelow
i0      = args.iZero    # Rename
mw      = args.mrgnW
mr      = args.mrgnR
mh      = args.mrgnH
flt     = args.filter
hmax    = args.hmax
rmax    = args.rmax
printOn    = args.printOn
printOnly  = args.printOnly

# Test comment
# Another one

if( flt[0] == None):  fltStr  = ' '
else:                 fltStr  = flt[0]+'-filtered: '


# Read the raster tile to be processed.
Rdict = readNumpyZTile(filename)
R = Rdict['R']
Rdims = np.array(np.shape(R))
ROrig = Rdict['GlobOrig']
print(' Rdims = {} '.format(Rdims))
print(' ROrig = {} '.format(ROrig))

# Set the zero level according to the given pixel value.
if(i0.count(None) == 0):
  print(' Zero Level: {} '.format(R[i0[0],i0[1]]))
  R0 = R[i0[0],i0[1]]
  R -= R0
  R[R<0.] = 0.

R = applyMargins( R , mw, mr, mh )

# Apply desired filters.
Rf =  np.zeros( np.shape(R) , float)
Rf =  filterAndScale(Rf, R, flt )

# Apply nans where fill values are desired.
Rf =  replaceByNans( Rf, na, nb)

if( rmax is not None ):
  idv = (Rf > rmax)
  Rf[idv]  = np.maximum( Rf[idv], R[idv] )


if( hmax ):
  Rf = np.minimum( hmax , Rf )

Rdict['R'] = Rf; Rdict['GlobOrig'] = ROrig


if( not args.printOnly ):
  saveTileAsNumpyZ( fileout, Rdict )

if( args.printOn or args.printOnly ):
  figDims = 13.*(Rdims[::-1].astype(float)/np.max(Rdims))
  #print('Sum = {}'.format(np.sum(Rf)))
  fig = plt.figure(num=1, figsize=figDims)
  fig = addImagePlot( fig, Rf, fltStr+fileout )

  plt.show()


R = Rf = None
