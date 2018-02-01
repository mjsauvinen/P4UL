#!/usr/bin/env python
import sys
import argparse
import ConfigParser
import numpy as np

from netcdfTools import *
from mapTools import readNumpyZTile


'''
Description:
This grid generates input files for PALM simulations following PALM Input Data Standard (PIDS) v1.7.

Author:
Sasu Karttunen
sasu.karttunen@helsinki.fi
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='mergeNetCdfFiles.py')
parser.add_argument("config", type=str, default=None, help="Name of the input config file.")
args = parser.parse_args()
#==========================================================#

config = ConfigParser.ConfigParser()
config.read(args.config)

'''
Read the configuration file, perform some basic checks on input files
and print the input parameters into standard output.
'''
print("Input configuration file {}:".format(args.config))

print("\nGlobal attributes:")
globalAttributes = config._sections['Global']
del globalAttributes['__name__']
for attr, val in globalAttributes.iteritems():
  print("{}: {}".format(attr,val))

'''
Topography configuration data
'''

print("\nTopography configuration:")
topoConfig = config._sections['Topography']
del topoConfig['__name__']

# Read topography data
oroDict = readNumpyZTile(config.get('Topography', 'orography'),verbose=False)
oroR = oroDict['R']
oroDPx = oroDict['dPx']
oroNPx = np.shape(oroR)

buildDict = readNumpyZTile(config.get('Topography', 'buildings'),verbose=False)
buildR = buildDict['R']
buildDPx = buildDict['dPx']
buildNPx = np.shape(buildR)
buildLOD = len(buildNPx)-1 # 1=2D height field, 2=3D mask
if (buildLOD==2):
  bMask = True
else:
  bMask = False

for attr, val in topoConfig.iteritems():
  print("{}: {}".format(attr,val))

print("buildings_lod: {}".format(buildLOD))

if (buildNPx[0:2] != oroNPx):
  raise ValueError("Horizontal dimensionality of buildings array doesn't match the orography: the shape of orography array is {} while the shape of buildings array is {}.".format(oroNPx, buildNPx))

'''
Surface classification data
'''

print("\nSurface configuration:")
surfConfig = config._sections['Surface']
del surfConfig['__name__']

buildIDDict = readNumpyZTile(config.get('Surface', 'building_id'),verbose=False)
buildIDR = buildIDDict['R']
buildIDdPx = buildIDDict['dPx']
buildIDNPx = np.shape(buildIDR)

for attr, val in surfConfig.iteritems():
  print("{}: {}".format(attr,val))

if (buildIDNPx != oroNPx):
  raise ValueError("Dimensionality of building ID array doesn't match the orography: the shape of the orography array is {} while the shape of the building ID array is {}.".format(oroNPx, buildIDNPx))



'''
Open the output netCDF files and write global attributes to them.
'''
print("\nGenerating output datasets...")
pidsStaticDS = netcdfOutputDataset("PIDS_STATIC")
setPIDSGlobalAtrributes(pidsStaticDS, globalAttributes)

'''
Write dimensions into PIDS_STATIC
'''
xvStatic = createCoordinateAxis(pidsStaticDS, oroNPx, oroDPx, 1, 'x', 'f4', 'm', True, False)
xvStatic.long_name = "distance to origin in x-direction"
yvStatic = createCoordinateAxis(pidsStaticDS, oroNPx, oroDPx, 0, 'y', 'f4', 'm', True, False)
yvStatic.long_name = "distance to origin in y-direction"

# Write z axis if 3D masks are used
if (bMask):
  zvStatic = createCoordinateAxis(pidsStaticDS, buildNPx, buildDPx, 2, 'z', 'f4', 'm', True, False)
  zvStatic.long_name = "distance to origin in z-direction"

'''
Write minimum required data into PIDS_STATIC
'''
oroR = np.swapaxes(oroR,0,1) # [x,y] -> [y,x]
oroNCVar = createNetcdfVariable(pidsStaticDS, oroR, 'orography_2D', 0, 'm', 'f4', ('y','x'), False, False, fill_value=-9999.9)
oroNCVar.long_name= "orography"

if (bMask):
  buildR = np.swapaxes(buildR,0,2) # [x,y,z] -> [z,y,x]
  buildNCVar = createNetcdfVariable(pidsStaticDS, buildR, 'buildings_3D', 0, 'm', 'b', ('z','y','x'), False, False, fill_value=-127)
else:
  buildNCVar = createNetcdfVariable(pidsStaticDS, buildR, 'buildings_2D', 0, 'm', 'f4', ('y','x'), False, False, fill_value=-9999.9)

buildNCVar.long_name= "buildings"
buildNCVar.lod = int(buildLOD)

buildIDNCVar = createNetcdfVariable(pidsStaticDS, buildIDR, 'building_id', 0, 'm', 'i4', ('y','x'), False, False, fill_value=-9999)
buildIDNCVar.long_name= "building_id"

'''
Write optional variables into PIDS_STATIC
'''



'''
Close all datasets.
'''
netcdfWriteAndClose(pidsStaticDS)
