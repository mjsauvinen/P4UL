#!/usr/bin/env python
import sys
import argparse
import ConfigParser
import numpy as np
import os.path

from netcdfTools import *
from pidsTools import *


'''
Description:
This grid generates input files for PALM simulations following PALM Input Data Standard (PIDS)

PIDS version: 1.9

TODOS:
- Extend to cover more variables
- Check that old and new dimensions match each other

Author:
Sasu Karttunen
sasu.karttunen@helsinki.fi
Institute for Atmospheric and Earth System Research (INAR) / Physics
University of Helsinki

'''

#==========================================================#
parser = argparse.ArgumentParser(prog='genPIDSInput.py')
parser.add_argument("config", type=str, default=None, help="Name of the input config file.")
args = parser.parse_args()
#==========================================================#

config = ConfigParser.ConfigParser()
config.read(args.config)

'''
Read the configuration file, perform some basic checks on input files
and print the input parameters into standard output.
'''
print("===== Configuration file {} =====".format(args.config))

print("\n== Global attributes ==")
newGlobalAttributes = True
globalAttributes = readConfigSection(config, 'Global')
if(globalAttributes is None):
  print("No global attributes specified.")
  globalAttributes={}
  newGlobalAttributes = False

print("\n== Topography ==")
topoConfig = readConfigSection(config, 'Topography')
if(topoConfig is None):
  print("No topography specified.")

print("\n== Surface classification ==")
surfConfig = readConfigSection(config, 'Surface')
if(surfConfig is None):
  print("No surface classification specified.")

print("\n== Vegetation ==")
vegConfig = readConfigSection(config, 'Vegetation')
if(vegConfig is None):
  print("No vegetation specified.")

'''
Open the output netCDF files and write global attributes to them.

Intialize a dictionaries for existing dimensions and variables. This is important
because netCDF4-python doesn't allow direct overwrite for these. For now we trust
that the user doesn't try to append new data that has different dimensions than the
what already exist in the original PIDS_STATIC.

TODO: Do not trust the user and write checks for every variable.
'''

print("\n===== Static input file PIDS_STATIC =====")
staticAppend = os.path.isfile("PIDS_STATIC")
if(staticAppend):
  pidsStaticDS = netcdfOutputDataset("PIDS_STATIC", mode="r+")
  staticDims = pidsStaticDS.dimensions.keys()
  staticVars = pidsStaticDS.variables.keys()
  # Remove dims from vars with a filter
  staticVars = filter(lambda key: key not in staticDims, staticVars)
  print("Existing PIDS_STATIC file found, using append/update mode.")

  if(len(staticDims)>0):
    print("Dimensions: " + ', '.join(str(key) for key in staticDims))
    print("Variables: " + ', '.join(str(key) for key in staticVars))
    print("Warning: using old dimensions for new variables")
  else:
    print("No existing dimensions or variables found.")
else:
  pidsStaticDS = netcdfOutputDataset("PIDS_STATIC", mode="w")
  staticDims = []
  staticVars = []
  print("No existing PIDS_STATIC found, creating a new one.")


setPIDSGlobalAtrributes(pidsStaticDS, globalAttributes)

'''
Write topography with it's dimensions into PIDS_STATIC
'''

if(topoConfig is not None):
  oroFile = readConfigVariable(config, 'Topography', 'orography')
  if(oroFile is not None and oroFile!=""):
    oroVar = processOrography(oroFile,pidsStaticDS,staticVars,staticDims)

  buildFile = readConfigVariable(config, 'Topography', 'buildings')
  if(buildFile is not None and buildFile!=""):
    buildVar = processBuildings(buildFile,pidsStaticDS,staticVars,staticDims)

'''
Write surface classification data into PIDS_STATIC
'''

if(surfConfig is not None):
  buildIDFile = readConfigVariable(config, 'Surface', 'building_id')
  if(buildIDFile is not None and buildIDFile!=""):
    buildIDVar = processBuildingIDs(buildIDFile,pidsStaticDS,staticVars,staticDims)

  pavementTypeFile = readConfigVariable(config, 'Surface', 'pavement_type')
  if(pavementTypeFile is not None and pavementTypeFile!=""):
    pavementTypeVar = processPavementType(pavementTypeFile,pidsStaticDS,staticVars,staticDims)


if(vegConfig is not None):
  ladFile = readConfigVariable(config, 'Vegetation', 'lad')
  if(ladFile is not None and ladFile!=""):
    ladVar = processLAD(ladFile,pidsStaticDS,staticVars,staticDims)

  vegetationTypeFile = readConfigVariable(config, 'Vegetation', 'vegetation_type')
  if(vegetationTypeFile is not None and vegetationTypeFile!=""):
    vegetationTypeVar = processVegetationType(vegetationTypeFile,pidsStaticDS,staticVars,staticDims)

'''
Close PIDS_STATIC
'''
netcdfWriteAndClose(pidsStaticDS, verbose=False)
print("...static input file closed.")
