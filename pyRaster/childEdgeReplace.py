#!/usr/bin/env python3
import sys
import argparse
import numpy as np
from mapTools import *
from utilities import filesFromList, writeLog
from plotTools import addImagePlot
import matplotlib.pyplot as plt
#import scipy.ndimage as sn
'''
Description:
Adjusts the terrain-topography, building and building id raster edges of a child domain such that the 
terrain+building contours in child edge match exactly with the contour in parent. This is to make sure 
that the child edge flow areas match exactly with the corresponding flow area on the parent side. 
Files must be processed one by one (use options -fc, -fp, and -fo). If there is no building information on 
the parent side, use the argument -npr (--noParent) and NaNs will be set to the child edge zones automatically.
The edge zone width can be set by the option -e, but the default value 1 parent cell width is recommended.
Tested 17.11.2022.

Author: Mia Aarnio
        mia.aarnio@fmi.fi
        Finnish Meteorological Institute
        
        Antti Hellsten
        antti.hellsten@fmi.fi
        Finnish Meteorological Institute
'''
#==========================================================#
parser = argparse.ArgumentParser(prog='childEdgeReplace.py',
      description='''Copies values from lower resolution parent raster and pastes them to higher resolution child raster outer edge.''')
parser.add_argument("-fc", "--child", metavar="CHILD", type=str,
                    help="Child domain raster data file (.npz).")
parser.add_argument("-fp", "--parent", metavar="PARENT", type=str, help="Parent domain raster data file (.npz).")
parser.add_argument("-fo", "--fileout", metavar="CHILD", type=str, help="Output modified child raster data file (.npz).")
parser.add_argument("-e", "--ebw", type=int, default = 1,help="Edge buffer width in parent-grid spacings (default = 1).")
parser.add_argument("-npr", "--noParent", action="store_true", default=False, help="No parent building related raster. Make NaN-filled edge to child raster. You must give parent terrain raster as -fp.")   
parser.add_argument("-p", "--printOn", action="store_true", default=False, help="Print the resulting raster data.")
parser.add_argument("-pp", "--printOnly", action="store_true", default=False,
                    help="Only print the resulting data. Don't save.")
parser.add_argument("-t", "--tolefactor", type=float, default = 0.000001,
                    help="Tolerance for grid matching tests, fraction of smallest grid spacing (default = 0.000001).") 
args = parser.parse_args()
writeLog(parser, args)

#==========================================================#
# Renaming ...
parentFile = args.parent 
childFile  = args.child
fileout    = args.fileout
ebw        = args.ebw
tolefactor = args.tolefactor
useNans    = args.noParent
printOn    = args.printOn
printOnly  = args.printOnly

# Read grid information from both files
#print('\n# - - - - - PARENT INFO - - - - - - - - #')
RdictParent  = readNumpyZTile(parentFile)   #pick up metadata t
Rparent      = RdictParent['R']    # pick up non-metadata e.g. the raster data
Rtype_parent = Rparent.dtype       # file datatype
nPxParent    = np.shape(Rparent)   # [j,i] grid size in both directions 
ROrigParent  = RdictParent['GlobOrig']  # Top left cell-center origin in global coordinates
dPxParent    = RdictParent['dPx']       #resolution [y, x]
RdictParent  = None                     #removing the original data file from memory 

#print(' Global origo (Top Left): [N,E] = [{}, {}]'.format(*ROrigParent))
print(' Size: [j,i] = [{}, {}]'.format(*nPxParent))
print(' Resolution: [dPy,dPx] = [{}, {}]'.format(*dPxParent))

#print('\n# - - - - - CHILD INFO - - - - - - - - #')  # all of these are python
# The following are python lists except Rchild (numpy-array) and Rtype
RdictChild  = readNumpyZTile(childFile) #pick up metadata 
Rchild      = RdictChild['R']           #pick up non-metadata e.g. the raster data 
Rtype_child = Rchild.dtype              # file  datatype
nPxChild    = np.shape(Rchild)          # [j,i] grid size in both direction
ROrigChild  = RdictChild['GlobOrig']    # Top left cell-center origin in global coordinates
dPxChild    = RdictChild['dPx']         #resolution [y, x]
RdictChild  = None                      #removing the original data file from memory 

# Check that both files are of same datatype:
if ( Rtype_child == Rtype_parent ):
  Rtype = Rtype_parent
else:
  sys.exit( ' Child and parent rasters are of different datatypes. Exiting...' )

#print(' Global origo (Top Left): [N,E] = [{}, {}]'.format(*ROrigChild))
print(' Size: [j,i] = [{}, {}]'.format(*nPxChild))
print(' Resolution: [dPy,dPx] = [{}, {}]'.format(*dPxChild))

# Check that parent and child are in correct order and resolutions: dPxChild & dPxParent
# Integer-divisibility test.  % is modulo 
tolerance = tolefactor * np.min(dPxChild)
modgsry0 = np.abs( dPxParent[0] %  dPxChild[0] )     
modgsry1 = np.abs( dPxChild[0] - dPxParent[0] % dPxChild[0] )
modgsrx0 = np.abs( dPxParent[1] %  dPxChild[1] )
modgsrx1 = np.abs( dPxChild[1] - dPxParent[1] % dPxChild[1] )

if ( modgsry0 > tolerance  and  modgsry1 > tolerance ):
  sys.exit(' Parent and Child grid resolutions in y are not compatible. Exiting...')
if ( modgsrx0 > tolerance  and  modgsrx1 > tolerance ):
  sys.exit(' Parent and Child grid resolutions in x are not compatible. Exiting...')

# Check that grid spacings are not the same
if ( np.abs( dPxParent[0] - dPxChild[0] ) < tolerance ): 
  sys.exit(' Parent and Child grid resolutions in y are the same. Exiting...')
if ( np.abs( dPxParent[1] - dPxChild[1] ) < tolerance ): 
  sys.exit(' Parent and Child grid resolutions in x are the same. Exiting...')  
  
# Ratio of grid spacings 0=j/y/N 1=i/x/E: warn if this is more than 3 indicating the
# rasters may have been given wrong.
# Note that gsr is of type integer.
gsr = np.zeros(2).astype(int)                            #create as a  numpy array
gsr = (np.round(dPxParent/dPxChild)).astype(int)
print(' Resolution ratio of Parent/Child :  [j,i]  = [{}, {}]'.format(*gsr))
if ( gsr[0] > 3.0 or gsr[1] > 3 ):
  print( ' Warning: Ratio of resolution of Parent to Child is large. Did you input correct files? ') 
  
# Child origo bottom left in Parent coordinates  
ChildOrigBL = np.zeros(2)
ChildOrigBL[0] =  ( ROrigChild[0] + dPxChild[0]/2.0 - nPxChild[0] * dPxChild[0] ) - \
                  ( ROrigParent[0] + dPxParent[0]/2.0 - nPxParent[0] * dPxParent[0] )
ChildOrigBL[1] =  ( ROrigChild[1] - dPxChild[1]/2.0 ) - (ROrigParent[1] - dPxParent[1]/2.0 )

# Child origo in Parent coordinates
print( ' Child origo bottom left corner:[j,i] = [{}, {}]'.format(*ChildOrigBL))

# Check that child is located inside the parent from left and south
if ( ChildOrigBL[0] < 0.0 or ChildOrigBL[1] < 0.0 ):
  sys.exit( ' Child not located inside the Parent. Exiting...')     

# Check that child is located all inside the parent from right and north
ChildURcorner = np.zeros(2)
ParentURcorner = np.zeros(2)
ChildURcorner = ChildOrigBL + dPxChild*nPxChild
ParentURcorner = dPxParent*nPxParent
print( 'ChildURcorner:[j,i] = [{}, {}]'.format(*ChildURcorner))
print( 'ParentURcorner:[j,i] = [{}, {}]'.format(*ParentURcorner))
if ( ChildURcorner[0] > ParentURcorner[0] or ChildURcorner[1] > ParentURcorner[1] ):        
  sys.exit( ' Child not fully located inside the Parent. Exiting...')
    
# If parentRLcorner = ChildRLcorner, print " parentRLcorner = ChildRLcorner, are you
# doing vertical nesting? if yes then this is right"
if ( ChildOrigBL[0] == 0.0 or ChildOrigBL[1] == 0.0 or
     ChildURcorner[0] == ParentURcorner[0] or ChildURcorner[1] == ParentURcorner[1] ):        
  print( ' Warning: Child edge located on Parent edge. Vertical nesting intended? ')

print(' Child Parent offset: [j,i] = [{}, {}]'.format(*ChildOrigBL))

# Child corners in parent indices
jCorner=np.zeros(2).astype(int)                   # oli .astype(np.int) jos tuli herjaa
iCorner=np.zeros(2).astype(int)                   # oli .astype(np.int) jos tuli herja

# jCorner[1] and iCorner[1] have +1 because of python slicing excludes the upper bounding
# element (creating RparentSlice1 below) 
jCorner[0] = np.round( ChildOrigBL[0]/dPxParent[0] ).astype(int)
jCorner[1] = np.round( jCorner[0] + nPxChild[0]/gsr[0] - 1 ).astype(int) + 1   
iCorner[0] = np.round( ChildOrigBL[1]/dPxParent[1] ).astype(int)
iCorner[1] = np.round( iCorner[0] + nPxChild[1]/gsr[1] - 1 ).astype(int) + 1
print( ' Child bottom left corner in parent indices: [j, i] = ', jCorner[0], iCorner[0] )
print( ' Child top right corner in parent indices: [j, i] = ', jCorner[1], iCorner[1] )

# Cut child area from parent at parent resolution
RparentSlice1 = Rparent[jCorner[0]:jCorner[1], iCorner[0]:iCorner[1]]
Rparent = None

# Create mapping from child indexing to RparentSlice indexing for refining
n1,e1 = np.ogrid[ 0:nPxChild[0] , 0:nPxChild[1] ]  # northing, easting 
n2,e2 = np.ogrid[ 0:nPxChild[0] , 0:nPxChild[1] ]  # northing, easting
# In the following, 0.5 is added to make sure that np.floor would never floor 
# an integer valued result represented by float down by one due to round-off error. 
n1 = np.floor((n1+0.5)/gsr[0]).astype(int)
e1 = np.floor((e1+0.5)/gsr[1]).astype(int)

# Create RparentSlice2 by refining RparentSlice1 to child resolution
#creating a new array in the size and resolution of the child
RparentSlice2 = np.zeros([nPxChild[0],nPxChild[1]], dtype=Rtype)
RparentSlice2[n2,e2] = RparentSlice1[n1,e1]

if ( useNans ):
   nanny  = np.nan    
   RparentSlice2.fill(nanny)   
   
# Substitute the interiors (all excluding the edge buffer zones) of Rchild into RparentSlice2
ebwChild = ebw * gsr
jsouth = ebwChild[0]
jnorth = nPxChild[0] - ebwChild[0]
ileft  = ebwChild[1]
iright = nPxChild[1] - ebwChild[1]
RparentSlice2[jsouth:jnorth,ileft:iright] = Rchild[jsouth:jnorth,ileft:iright] 

# Output 
Rdict = dict()
Rdict['GlobOrig'] = ROrigChild
Rdict['dPx']      = dPxChild
Rdict['R']        = RparentSlice2
    
if( not args.printOnly ):
  saveTileAsNumpyZ( fileout, Rdict )

if( args.printOn or args.printOnly ):
  fig = plt.figure(num=1, figsize=(13.5,13.5))
  fig = addImagePlot( fig, RparentSlice2, fileout )

  plt.show()
      


