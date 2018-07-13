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

WARNING: I haven't had time to test this thoroughly so please check your PIDS files before using them.

TODOS:
- Extend to cover more variables
- Check that old and new dimensions match each other
- Skip file if nothing to do
- Write more comments to code, especially pidsTools.py.

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

print("\n== Chemistry ==")
chemConfig = readConfigSection(config, 'Chemistry')
if(chemConfig is None):
  print("No chemistry specified.")

print("\n== Aerosols ==")
aeroConfig = readConfigSection(config, 'Aerosols')
if(aeroConfig is None):
  print("No aerosols specified.")

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
print("Static output file PIDS_STATIC closed.")

'''
Move on to PIDS_CHEM
'''
if(chemConfig is not None):
  print("\n===== Chemistry input file PIDS_CHEM =====")
  chemAppend = os.path.isfile("PIDS_CHEM")
  if(chemAppend):
    pidsChemDS = netcdfOutputDataset("PIDS_CHEM", mode="r+")
    chemDims = pidsChemDS.dimensions.keys()
    chemVars = pidsChemDS.variables.keys()
    # Remove dims from vars with a filter
    chemVars = filter(lambda key: key not in chemDims, chemVars)
    print("Existing PIDS_CHEM file found, using append/update mode.")

    if(len(chemDims)>0):
      print("Dimensions: " + ', '.join(str(key) for key in chemDims))
      print("Variables: " + ', '.join(str(key) for key in chemVars))
      print("Warning: using old dimensions for new variables")
    else:
      print("No existing dimensions or variables found.")
  else:
    pidsChemDS = netcdfOutputDataset("PIDS_CHEM", mode="w")
    chemDims = []
    chemVars = []
    print("No existing PIDS_CHEM found, creating a new one.")

  setPIDSGlobalAtrributes(pidsChemDS, globalAttributes)

  emiCatInds = readConfigVariable(config, 'Chemistry', 'emission_category_index')
  if(emiCatInds is not None and emiCatInds!=""):
    emiCatIndsVar = processEmissionCategoryIndices(emiCatInds, pidsChemDS, chemVars, chemDims)

  emiInds = readConfigVariable(config, 'Chemistry', 'emission_index')
  if(emiInds is not None and emiInds!=""):
    emiInds = processEmissionIndices(emiInds, pidsChemDS, chemVars, chemDims)

  emiCatName = readConfigVariable(config, 'Chemistry', 'emission_category_name')
  if(emiCatName is not None and emiCatName!=""):
    emiCatName = processEmissionCategoryNames(emiCatName, pidsChemDS, chemVars, chemDims)

  emiSpeName = readConfigVariable(config, 'Chemistry', 'emission_species_name')
  if(emiSpeName is not None and emiSpeName!=""):
    emiSpeName = processEmissionSpeciesNames(emiSpeName, pidsChemDS, chemVars, chemDims)

  emiTimeFactors = readConfigVariable(config, 'Chemistry', 'emission_time_factors')
  emiTimeFactorsLod = readConfigVariable(config, 'Chemistry', 'emission_time_factors_lod')
  if(emiTimeFactors is not None and emiTimeFactors!=""):
    if(emiTimeFactorsLod==""):
      emiTimeFactorsLod = None
    emiTimeFactors = processEmissionTimeFactors(emiTimeFactors, emiTimeFactorsLod, pidsChemDS, chemVars, chemDims)


  netcdfWriteAndClose(pidsChemDS, verbose=False)
  print("Chemistry output file PIDS_CHEM closed.")

'''
Move on to PIDS_AERO
'''

if(aeroConfig is not None):
  print("\n===== Aerosol input file PIDS_AERO =====")
  aeroAppend = os.path.isfile("PIDS_AERO")
  if(aeroAppend):
    pidsAeroDS = netcdfOutputDataset("PIDS_AERO", mode="r+")
    aeroDims = pidsAeroDS.dimensions.keys()
    aeroVars = pidsAeroDS.variables.keys()
    # Remove dims from vars with a filter
    aeroVars = filter(lambda key: key not in aeroDims, aeroVars)
    print("Existing PIDS_AERO file found, using append/update mode.")

    if(len(aeroDims)>0):
      print("Dimensions: " + ', '.join(str(key) for key in aeroDims))
      print("Variables: " + ', '.join(str(key) for key in aeroVars))
      print("Warning: using old dimensions for new variables")
    else:
      print("No existing dimensions or variables found.")
  else:
    pidsAeroDS = netcdfOutputDataset("PIDS_AERO", mode="w")
    aeroDims = []
    aeroVars = []
    print("No existing PIDS_AERO found, creating a new one.")

  setPIDSGlobalAtrributes(pidsAeroDS, globalAttributes)

  compAeroStr = readConfigVariable(config, 'Aerosols', 'composition_aerosol')
  if(compAeroStr is not None and compAeroStr!=""):
    compAeroVar = processCompositionAerosol(compAeroStr, pidsAeroDS, aeroVars, aeroDims)



  netcdfWriteAndClose(pidsAeroDS, verbose=False)
  print("Aerosol output file PIDS_AERO closed.")

'''
Close PIDS_AERO
'''
