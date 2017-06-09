#!/usr/bin/env python
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
parser = argparse.ArgumentParser(prog='distributeValuesToAreas.py', description='''Calculates child domain's offset to parent domain in local and global coordinates.''')
parser.add_argument("-fc", "--child", metavar="CHILD" ,type=str, help="Child domain raster data file (.npz).")
parser.add_argument("-fp", "--parent", metavar="PARENT" ,type=str, help="Parent domain raster data file (.npz).")
args = parser.parse_args()
writeLog( parser, args )

#==========================================================#

# Read grid information from both files
RdictParent = readNumpyZTile(args.parent)
nPxParent = np.shape(RdictParent['R'])
ROrigParent = RdictParent['GlobOrig']
dPxParent = RdictParent['dPx']
gridRot = RdictParent['gridRot']
RdictParent = None
print(' Global origo: [N,E] = [{}, {}]'.format(*ROrigParent))
print(' Size: [X,Y] = [{}, {}]'.format(*nPxParent))
print(' Resolution: [dPx,dPy] = [{}, {}]'.format(*dPxParent))
print(' Grid rotation: [deg] = {}'.format(gridRot/(np.pi/180.)));print('')

RdictChild = readNumpyZTile(args.child)
nPxChild = np.shape(RdictChild['R'])
ROrigChild = RdictChild['GlobOrig']
dPxChild = RdictChild['dPx']
gridRotChild = RdictChild['gridRot']
RdictChild = None

if (gridRot != gridRotChild):
  sys.exit('Rotations of parent and child domain don\'t match! Exiting...')

print(' Global origo: [N,E] = [{}, {}]'.format(*ROrigChild))
print(' Size: [X,Y] = [{}, {}]'.format(*nPxChild))
print(' Resolution: [dPx,dPy] = [{}, {}]'.format(*dPxChild))
print(' Grid rotation: [deg] = {}'.format(gridRotChild/(np.pi/180.)));print('')

# Calculate bottom left origos
ROrigParentBL = np.zeros(2);ROrigChildBL = np.zeros(2)
ROrigParentBL[0] = ROrigParent[0] - nPxParent[1]*dPxParent[1]*np.cos(gridRot)
ROrigParentBL[1] = ROrigParent[1] + nPxParent[0]*dPxParent[0]*np.sin(gridRot)

ROrigChildBL[0] = ROrigChild[0] - nPxChild[1]*dPxChild[1]*np.cos(gridRot)
ROrigChildBL[1] = ROrigChild[1] + nPxChild[0]*dPxChild[0]*np.sin(gridRot)

# Offset of global origo coordinates
OrigOffsetGlobal = ROrigChildBL - ROrigParentBL
print(' Bottom left origo offsets:')
print(' Global: [N,E] = [{}, {}]'.format(*OrigOffsetGlobal))

# Calculate local origo offset in the PALM grid of the parent domain by rotation
# Note that the offset is given in domain's grid [X,Y]!
OrigOffset = rotatePoint(OrigOffsetGlobal[::-1], -gridRot, 1/dPxParent)
print(' Parent domain\'s grid: [X,Y] = [{}, {}]'.format(*OrigOffset))

# Help the user to move the child domain to match the parent's grid
if (not(OrigOffset[0].is_integer() and OrigOffset[1].is_integer())):
  NewOffset = map(int, OrigOffset)
  # Rotate and scale back to original global top left origo
  NewOffset = rotatePoint(NewOffset, gridRot, dPxParent)[::-1]
  NewOrigo = NewOffset + ROrigParentBL
  NewOrigo[0] += nPxChild[1]*dPxChild[1]*np.cos(gridRot)
  NewOrigo[1] -= nPxChild[0]*dPxChild[0]*np.sin(gridRot)
  print(NewOrigo)
  childAdj = NewOrigo-ROrigChild
  
  print(' WARNING: Child\'s origo doesn\'t match to the parent\'s grid.')
  print(' Move the child domain area by [{}, {}] in N,E grid to align the origo.'.format(*childAdj))
else:
  # Check if the grid dimensions match, i.e. the edges align with the parent grid
  xRatio=nPxChild[0]*dPxChild[0]/dPxParent[0]
  yRatio=nPxChild[1]*dPxChild[1]/dPxParent[1]
  if (not(xRatio.is_integer() and yRatio.is_integer())):
    print(' WARNING: Child domain\'s grid edges don\'t align with the parent. Check your resolutions and dimensions.')