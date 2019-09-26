#!/usr/bin/env python3
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

print("\n== General configuration ==")
outputConfig = readConfigSection(config, 'General')
if (outputConfig is None):
  print("No output configuration specified, using defaults.")

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
Format output file names.
'''
nest_id = ""
if (outputConfig is not None):
  nestID = readConfigVariable(config, 'General', 'nest_id')
  if(nestID is not None and nestID!=""):
    nest_id = "_N{:0>2}".format(int(nestID))

'''
Open the output netCDF files and write global attributes to them.

Intialize a dictionaries for existing dimensions and variables. This is important
because netCDF4-python3 doesn't allow direct overwrite for these. For now we trust
that the user doesn't try to append new data that has different dimensions than the
what already exist in the original PIDS_STATIC.

TODO: Do not trust the user and write checks for every variable.
'''
# Check if there is any reason to write PIDS_STATIC.
# This is a bit ugly.
if(not all(v is None for v in [topoConfig, surfConfig, vegConfig])):
  pidsStaticFN = "PIDS_STATIC" + nest_id
  print("\n===== Static input file "+ pidsStaticFN +" =====")
  skipStaticFile = False
  staticAppend=False
  # Check what the user really wants if an existing output file is found
  if (os.path.isfile(pidsStaticFN)):
    while True:
      outputMode = str(raw_input("Existing "+ pidsStaticFN +" file found. Overwrite, append, skip or exit? [o/a/s/e] "))
      if (outputMode in ["e","E"]):
        print("Exit.")
        exit()
      elif (outputMode in ["s","S"]):
        skipStaticFile=True
        print("Skipping static input.")
        break
      elif (outputMode in ["o","O"]):
        staticAppend=False
        print("Overwriting an existing "+ pidsStaticFN +" file.")
        break
      elif (outputMode in ["a","A"]):
        staticAppend=True
        break
      else:
        print("Invalid selection.")
        continue
  else:
    print("No existing "+ pidsStaticFN +" found, creating a new one.")

  if (not skipStaticFile):
    if(staticAppend):
      print("Appending into existing "+ pidsStaticFN +" file.")
      pidsStaticDS = netcdfOutputDataset(pidsStaticFN, mode="r+")
      staticDims = pidsStaticDS.dimensions.keys()
      staticVars = pidsStaticDS.variables.keys()
      # Remove dims from vars with a filter
      staticVars = filter(lambda key: key not in staticDims, staticVars)

      if(len(staticDims)>0):
        print("Dimensions: " + ', '.join(str(key) for key in staticDims))
        print("Variables: " + ', '.join(str(key) for key in staticVars))
        print("Warning: using old dimensions for new variables")
      else:
        print("No existing dimensions or variables found.")
    else:
      pidsStaticDS = netcdfOutputDataset(pidsStaticFN, mode="w")
      staticDims = []
      staticVars = []



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

      waterTypeFile = readConfigVariable(config, 'Surface', 'water_type')
      if(waterTypeFile is not None and waterTypeFile!=""):
        waterTypeVar = processWaterType(waterTypeFile,pidsStaticDS,staticVars,staticDims)

      soilTypeFile = readConfigVariable(config, 'Surface', 'soil_type')
      if(soilTypeFile is not None and soilTypeFile!=""):
        soilTypeVar = processSoilType(soilTypeFile,pidsStaticDS,staticVars,staticDims)

      streetTypeFile = readConfigVariable(config, 'Surface', 'street_type')
      if(streetTypeFile is not None and streetTypeFile!=""):
        streetTypeVar = processStreetType(streetTypeFile,pidsStaticDS,staticVars,staticDims)


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
    print("Static output file "+ pidsStaticFN +" closed.")

'''
Move on to PIDS_CHEM
'''
if(chemConfig is not None):
  pidsChemFN = "PIDS_CHEM" + nest_id
  print("\n===== Chemistry input file "+ pidsChemFN +" =====")
  skipChemFile=False
  chemAppend=False
  if (os.path.isfile(pidsChemFN)):
    while True:
      outputMode = str(raw_input("Existing "+ pidsChemFN +" file found. Overwrite, append, skip or exit? [o/a/s/e] "))
      if (outputMode in ["e","E"]):
        print("Exit.")
        exit()
      elif (outputMode in ["s","S"]):
        skipChemFile=True
        print("Skipping chemistry input.")
        break
      elif (outputMode in ["o","O"]):
        chemAppend=False
        print("Overwriting an existing "+ pidsChemFN +" file.")
        break
      elif (outputMode in ["a","A"]):
        chemAppend=True
        break
      else:
        print("Invalid selection.")
        continue
  else:
    print("No existing "+ pidsChemFN +" found, creating a new one.")

  if (not skipChemFile):
    if(chemAppend):
      print("Appending into existing "+ pidsChemFN +" file.")
      pidsChemDS = netcdfOutputDataset(pidsChemFN, mode="r+")
      chemDims = pidsChemDS.dimensions.keys()
      chemVars = pidsChemDS.variables.keys()
      # Remove dims from vars with a filter
      chemVars = filter(lambda key: key not in chemDims, chemVars)

      if(len(chemDims)>0):
        print("Dimensions: " + ', '.join(str(key) for key in chemDims))
        print("Variables: " + ', '.join(str(key) for key in chemVars))
        print("Warning: using old dimensions for new variables")
      else:
        print("No existing dimensions or variables found.")
    else:
      pidsChemDS = netcdfOutputDataset(pidsChemFN, mode="w")
      chemDims = []
      chemVars = []

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
    print("Chemistry output file "+ pidsChemFN +" closed.")

'''
Move on to PIDS_AERO
'''

if(aeroConfig is not None):
  pidsAeroFN = "PIDS_SALSA" + nest_id
  print("\n===== Aerosol input file "+ pidsAeroFN +" =====")
  skipAeroFile = False
  aeroAppend=False
  if (os.path.isfile(pidsAeroFN)):
    while True:
      outputMode = str(raw_input("Existing "+ pidsAeroFN +" file found. Overwrite, append, skip or exit? [o/a/s/e] "))
      if (outputMode in ["e","E"]):
        print("Exit.")
        exit()
      elif (outputMode in ["s","S"]):
        skipAeroFile=True
        print("Skipping aerosol input.")
        break
      elif (outputMode in ["o","O"]):
        aeroAppend=False
        print("Overwriting an existing "+ pidsAeroFN +" file.")
        break
      elif (outputMode in ["a","A"]):
        aeroAppend=True
        break
      else:
        print("Invalid selection.")
        continue
  else:
    print("No existing "+ pidsAeroFN +" found, creating a new one.")

  if (not skipAeroFile):
    if(aeroAppend):
      pidsAeroDS = netcdfOutputDataset(pidsAeroFN, mode="r+")
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
      pidsAeroDS = netcdfOutputDataset(pidsAeroFN, mode="w")
      aeroDims = []
      aeroVars = []

    setPIDSGlobalAtrributes(pidsAeroDS, globalAttributes)

    compAeroStr = readConfigVariable(config, 'Aerosols', 'composition_aerosol')
    if(compAeroStr is not None and compAeroStr!=""):
      compAeroVar = processCompositionAerosol(compAeroStr, pidsAeroDS, aeroVars, aeroDims)

    aeroEmiStr = readConfigVariable(config, 'Aerosols', 'aerosol_emission_values')
    if (aeroEmiStr is not None and aeroEmiStr!=""):
      aeroSourceArea = readConfigVariable(config, 'Aerosols', 'aerosol_source_area')
      if(aeroSourceArea is not None and aeroSourceArea!=""):
        aeroEmissionDmid= readConfigVariable(config, 'Aerosols', 'aerosol_emission_dmid')
        if(aeroEmissionDmid is not None and aeroEmissionDmid!=""):
          aeroEminVar = processAerosolEmissionValues(aeroEmiStr, aeroSourceArea, aeroEmissionDmid, pidsAeroDS, aeroVars, aeroDims)

    emiMassFracsVar = processEmissionMassFractions(pidsAeroDS)

    emiCatInds = readConfigVariable(config, 'Chemistry', 'emission_category_index')
    if(emiCatInds is not None and emiCatInds!=""):
      emiCatIndsVar = processEmissionCategoryIndices(emiCatInds, pidsAeroDS, aeroVars, aeroDims)

    emiInds = readConfigVariable(config, 'Chemistry', 'emission_index')
    if(emiInds is not None and emiInds!=""):
      emiInds = processEmissionIndices(emiInds, pidsAeroDS, aeroVars, aeroDims)

    emiCatName = readConfigVariable(config, 'Chemistry', 'emission_category_name')
    if(emiCatName is not None and emiCatName!=""):
      emiCatName = processEmissionCategoryNames(emiCatName, pidsAeroDS, aeroVars, aeroDims)

    emiSpeName = readConfigVariable(config, 'Chemistry', 'emission_species_name')
    if(emiSpeName is not None and emiSpeName!=""):
      emiSpeName = processEmissionSpeciesNames(emiSpeName, pidsAeroDS, aeroVars, aeroDims)



    netcdfWriteAndClose(pidsAeroDS, verbose=False)
    print("Aerosol output file "+ pidsAeroFN +" closed.")

'''
Close PIDS_AERO
'''
