#!/usr/bin/env python
import argparse
import numpy as np
from mapTools import *
from plotTools import addImagePlot
from utilities import openStlFile, closeStlFile, writeStlFacet
import matplotlib.pyplot as plt
''' 
Description:


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='rasterToSTL.py')
parser.add_argument("-f", "--filename",type=str, help="Name of the .npz data file.")
parser.add_argument("-fo", "--fileout",type=str, help="Name of output STL file.",\
  default="topography.stl")
parser.add_argument("-p", "--printOn", help="Print the resulting raster data.",\
  action="store_true", default=False) 
parser.add_argument("-pp", "--printOnly", help="Only print the resulting data. Don't save.",\
  action="store_true", default=False) 
args = parser.parse_args() 
writeLog( parser, args )
#==========================================================#

filename  = args.filename
solidname = args.fileout


Rdict = readNumpyZTile( filename )
R = Rdict['R']
Rdims = np.array(np.shape(R))
ROrig = Rdict['GlobOrig']
dPx = Rdict['dPx']
Rdict = None
dPx = entry2Int( dPx )

print ' Rdims = {} '.format(Rdims)
print ' ROrig = {} '.format(ROrig)

if( ROrig[0] == 0 ):
  print ' Resetting y-origin. '
  ROrig[0] = Rdims[0]*dPx

nv = np.zeros(3,float)
v1 = np.zeros(3,float)
v2 = np.zeros(3,float)
v3 = np.zeros(3,float)



'''
   Upper
X******X
* *    *
*  *   *
*   *  * 
*    * * 
X******X
 Lower
'''

# Open the STL file and write the header.
fw = openStlFile( solidname )

for irow in xrange(Rdims[0]-1):
  dy = irow*dPx
  #print ' Processing row {} ... '.format(irow)
  for jcol in xrange(Rdims[1]-1):
    dx = jcol*dPx
    # Lower diagonal triangle.
    v1[0]=ROrig[1]+dx      ;  v1[1]=ROrig[0]-dy      ;  v1[2]=R[irow,jcol]
    v2[0]=ROrig[1]+dx      ;  v2[1]=ROrig[0]-(dy+dPx);  v2[2]=R[irow+1,jcol]
    v3[0]=ROrig[1]+(dx+dPx);  v3[1]=ROrig[0]-(dy+dPx);  v3[2]=R[irow+1,jcol+1]
    fw = writeStlFacet(fw, nv, v1, v2, v3 )
    # Upper diagonal triangle.
    # v1 -- same as above
    v2=v3.copy()
    v3[0]=ROrig[1]+(dx+dPx);  v3[1]=ROrig[0]-dy      ;  v3[2]=R[irow,jcol+1]
    fw = writeStlFacet(fw, nv, v1, v2, v3 )
    
closeStlFile( fw, solidname )

