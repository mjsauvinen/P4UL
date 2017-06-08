#!/usr/bin/env python
import sys
import argparse
import numpy as np
from mapTools import *
from utilities import writeLog

''' 
Description:
Calculates child domain's offset to parent domain.
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
ROrigParent[0] -= nPxParent[1]*dPxParent[1]*np.cos(gridRot)
ROrigParent[1] += nPxParent[1]*dPxParent[1]*np.sin(gridRot)

ROrigChild[0] -= nPxChild[1]*dPxChild[1]*np.cos(gridRotChild)
ROrigChild[1] += nPxChild[1]*dPxChild[1]*np.sin(gridRotChild)

# Offset of global origo coordinates
OrigOffsetGlobal = ROrigChild - ROrigParent
print(' Bottom left origo offsets:')
print(' Global: [N,E] = [{}, {}]'.format(*OrigOffsetGlobal))

# Calculate local origo offset in the PALM grid of the parent domain by rotation
# Note that the offset is given in domain's grid [X,Y]!
OrigOffset = [(OrigOffsetGlobal[1]*np.cos(-gridRot)-OrigOffsetGlobal[0]*np.sin(-gridRot))/dPxParent[0],\
              (OrigOffsetGlobal[1]*np.sin(-gridRot)+OrigOffsetGlobal[0]*np.cos(-gridRot))/dPxParent[1]]

print(' Parent domain\'s grid: [X,Y] = [{}, {}]'.format(*OrigOffset)))
