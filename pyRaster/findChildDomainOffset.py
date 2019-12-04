#!/usr/bin/env python3
import sys
import argparse
import numpy as np
from mapTools import *
from utilities import writeLog

'''
Description:
Calculates child domain's offset to parent domain.

NOTE: This script uses [N,E] coordinates for global coordinates and [X,Y] for
coordinates relative to domain's origo. This might be confusing at certain
parts and requires special care when editing the file.
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='findChildDomainOffset.py', description='''Calculates child domain's offset to parent domain in local and global coordinates.''')
parser.add_argument("-fc", "--child", metavar="CHILD", type=str, help="Child domain raster data file (.npz).")
parser.add_argument("-fp", "--parent", metavar="PARENT", type=str, help="Parent domain raster data file (.npz).")
args = parser.parse_args()
writeLog(parser, args)

#==========================================================#
# Renaming ...
parentFile = args.parent 
childFile  = args.child



# Read grid information from both files
RdictParent = readNumpyZTile(parentFile)
nPxParent   = np.shape(RdictParent['R'])
ROrigParent = RdictParent['GlobOrig']  # Top left cell-center origin
dPxParent   = RdictParent['dPx']
# Grid rotation may not be present in the dictionary
if( 'gridRot' in RdictParent.keys()): gridRot = RdictParent['gridRot'] 
RdictParent = None

print('\n# - - - - - PARENT INFO - - - - - - - - #')
print(' Global origo (Top Left): [N,E] = [{}, {}]'.format(*ROrigParent))
print(' Size: [N,E] = [{}, {}]'.format(*nPxParent))
print(' Resolution: [dPy,dPx] = [{}, {}]'.format(*dPxParent))
print(' Grid rotation: [deg] = {}\n'.format(gridRot / (np.pi / 180.)))

RdictChild = readNumpyZTile(childFile)
nPxChild   = np.shape(RdictChild['R'])
ROrigChild = RdictChild['GlobOrig']
dPxChild   = RdictChild['dPx']
if( 'gridRot' in RdictChild.keys()): gridRotChild = RdictChild['gridRot']
RdictChild = None

print('\n# - - - - - CHILD INFO - - - - - - - - #')
print(' Global origo (Top Left): [N,E] = [{}, {}]'.format(*ROrigChild))
print(' Size: [N,E] = [{}, {}]'.format(*nPxChild))
print(' Resolution: [dPy,dPx] = [{}, {}]'.format(*dPxChild))
print(' Grid rotation: [deg] = {}'.format(gridRotChild / (np.pi / 180.))); print('')

if (gridRot != gridRotChild):
  sys.exit('Rotations of parent and child domain don\'t match! Exiting...')

# Calculate bottom left corner origos
ROrigParentBL = ROrigParent.copy()
ROrigParentBL[0] -= ((nPxParent[0]-1) * dPxParent[0] + dPxParent[0]/2.)
ROrigParentBL[1] -= dPxParent[1]/2.


# We have to transform the child origin to align with the parents local N,E coords.
# We effectively undo the rotation to XRTM, YTRM in extractDomainFromTile.py.
ROrigChildT = rotatePoint(ROrigParent, ROrigChild, -gridRot)
print(' Rotated ROrigChild = {}'.format( ROrigChildT )) 

ROrigChildBL = ROrigChildT.copy()
ROrigChildBL[0] -= ((nPxChild[0]-1) * dPxChild[0] + dPxChild[0]/2.)
ROrigChildBL[1] -= dPxChild[1]/2.
print('\n Bottom Left Origos.\n Parent = {} vs.\n Child = {}'.format(ROrigParentBL,ROrigChildBL))

# Offset of global bottom left origo coordinates
OrigOffset = np.round( ROrigChildBL - ROrigParentBL , decimals=1 )

print('\n# - - - - - OFFSETS - - - - - - - - #') 
print(' Bottom left origo offsets:')
OrigOffsetLocal = np.round( OrigOffset / dPxParent , decimals=1 )
print(' Parent domain\'s grid: [dN,dE] = [{}, {}]  \n'.format(*OrigOffset))
print(' Enter this into PARIN -->    {:.1f}, {:.1f} \n'.format(OrigOffset[1],OrigOffset[0]))
print(' Pixels in parent domain\'s grid: [jN,iE]= [{}, {}]'.format(*OrigOffsetLocal))

# Help the user to move the child domain to match the parent's grid
if (not(OrigOffsetLocal[0].is_integer()) and not(OrigOffsetLocal[1].is_integer()) ):
  print(' WARNING: The offset in pixel count is not an integer: {} '.format(OrigOffsetLocal))

else:
  # Check if the grid dimensions match, i.e. the edges align with the parent grid
  xRatio = np.round( nPxChild[1] * dPxChild[1] / dPxParent[1] , decimals=1)
  yRatio = np.round( nPxChild[0] * dPxChild[0] / dPxParent[0] , decimals=1)
  if (not(xRatio.is_integer() and yRatio.is_integer())):
    print(' WARNING: Child domain\'s grid edges don\'t align with the parent. Check your resolutions and dimensions.')
  else:
    print(' Child\'s grid aligns with the parent\'s grid.')
