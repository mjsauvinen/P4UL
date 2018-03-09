#!/usr/bin/env python
import sys
import argparse
import numpy as np
from mapTools import *
from utilities import filesFromList, writeLog
from plotTools import addImagePlot
import matplotlib.pyplot as plt
''' 
Description:


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='mergeGridTiles.py')
parser.add_argument("-f", "--filenames",type=str, help="List of ascii grid tile files.", \
   nargs='+')
parser.add_argument("-p", "--printOn", help="Print the extracted tile.",\
  action="store_true", default=False) 
parser.add_argument("-pp", "--printOnly", help="Only print the extracted tile. Don't save.",\
  action="store_true", default=False) 
parser.add_argument("-s", "--scale",type=float,\
   help="Scale factor for the output. Default=1.", default=1.)   
args = parser.parse_args() 
writeLog( parser, args )
#==========================================================#

# Check the file formats.
ascii = False; npz = False
for filename in args.filenames:
  suffix = filename[-3:]
  if  (suffix == 'asc'):
    ascii = True
  elif(suffix == 'npz'):
    npz   = True
  else:
    print ' File {} appears to be neither .asc or .npz format. Exiting.'.format(filename)
    sys.exit(1)
    
if( ascii and npz ):
  print ' File formats (.asc/.npz) should not be mixed. Exiting.'
  sys.exit(1)


# Initialize a list for directories
dictList = []
idc = 0
fileout = str()

for filename in args.filenames:
  if( ascii ):
    dataDict, idc = readAsciiGridHeader( filename, idc )
  elif( npz ):
    dataDict, idc = readNumpyZGridData( filename, idc )
    
  dictList.append( dataDict )
  dataDict = None
  fstr = filename.split('.')[0]
  fileout += fstr[-1]   # Gather the last character/number from the name.

fstr = filename.split('.')[0]
fileout = fstr[:-1]+'_'+fileout  # Complete the filename.
dPx = resolutionFromDicts( dictList )

#print ' dictList = {} '.format(dictList)  

gIJ, XOrig, Mrows, Mcols = arrangeTileGrid( dictList, [ascii,npz] )
print 'gIJ : {} '.format( gIJ )

Rdict = compileTileGrid( dictList, gIJ , Mrows, Mcols, [ascii, npz] )
Rdims = np.array(np.shape(Rdict['R']))
Rdict['GlobOrig'] = XOrig
Rdict['dPx'] = dPx

if(not args.printOnly ):
  saveTileAsNumpyZ( fileout, Rdict )

if( args.printOn or args.printOnly ):
  R = Rdict['R']
  figDims = 13.*(Rdims[::-1].astype(float)/np.max(Rdims))
  fig = plt.figure(num=1, figsize=figDims)
  fig = addImagePlot( fig, R, 'Combined Tiles: '+fileout )
  plt.show()
  R = None

