#!/usr/bin/env python3
import sys
import argparse
import configparser
import numpy as np
import os.path

from netcdfTools import *
from pidsTools import *


'''
Description:
This grid generates input files for PALM simulations following PALM Input Data Standard (PIDS)

PIDS version: 1.9

See the documentation: https://palm.muk.uni-hannover.de/trac/wiki/doc/app/iofiles/pids

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

UPDATE:
- Mikko: In python3 the ConfigParser library is renamed configparser
- Jani Stromberg: Added building type and surface fractions
- Mona Kurppa: Added and modified chemistry and salsa variables
- Jukka-Pekka Keskinen: Added support for crs information
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='genPIDSInput.py')
parser.add_argument("config", type=str, default=None, help="Name of the input config file.")
args = parser.parse_args()
#==========================================================#

config = configparser.ConfigParser()
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

print("\n== CRS ==")
crsConfig = readConfigSection(config, 'CRS')
if(crsConfig is None):
  print("No crs specified. Default will be used (ETRS89, UTM zone from given latitude).")

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
salsaConfig = readConfigSection(config, 'Aerosols')
if (salsaConfig is None):
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
      outputMode = str(input("Existing "+ pidsStaticFN +" file found. Overwrite, append, skip or exit? [o/a/s/e] "))
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
      dso_s = netcdfOutputDataset(pidsStaticFN, mode="r+")
      dims_s = dso_s.dimensions.keys()
      vars_s = dso_s.variables.keys()
      # Remove dims from vars with a filter
      vars_s = filter(lambda key: key not in dims_s, vars_s)

      if(len(dims_s)>0):
        print("Dimensions: " + ', '.join(str(key) for key in dims_s))
        print("Variables: " + ', '.join(str(key) for key in vars_s))
        print("Warning: using old dimensions for new variables")
      else:
        print("No existing dimensions or variables found.")
    else:
      dso_s = netcdfOutputDataset(pidsStaticFN, mode="w")
      dims_s = []
      vars_s = []



    setPIDSGlobalAtrributes(dso_s, globalAttributes)

    if (crsConfig is not None):
      processCRS(dso_s, crsConfig)

    '''
    Write topography with its dimensions into PIDS_STATIC
    '''

    if(topoConfig is not None):
      fname = readConfigVariable(config, 'Topography', 'orography')
      if(fname is not None and fname!=""):
        oroVar = processOrography(fname,dso_s,vars_s,dims_s)
      fname = None  
      

      fname = readConfigVariable(config, 'Topography', 'buildings')
      if(fname is not None and fname!=""):
        buildVar = processBuildings(fname,dso_s,vars_s,dims_s)
      fname = None

    '''
    Write surface classification data into PIDS_STATIC
    '''

    if(surfConfig is not None):
      fname = readConfigVariable(config, 'Surface', 'building_id')
      if(fname is not None and fname!=""):
        buildIDVar = processBuildingIDs(fname,dso_s,vars_s,dims_s)
      fname = None

      fname = readConfigVariable(config, 'Surface', 'building_type')
      if (fname is not None and fname != ""):
        buildingTypeVar = processBuildingType(fname, dso_s, vars_s, dims_s)
      fname = None

      fname = readConfigVariable(config, 'Surface', 'pavement_type')
      if(fname is not None and fname!=""):
        pavementTypeVar = processPavementType(fname,dso_s,vars_s,dims_s)
      fname = None

      fname = readConfigVariable(config, 'Surface', 'water_type')
      if(fname is not None and fname!=""):
        waterTypeVar = processWaterType(fname,dso_s,vars_s,dims_s)
      fname = None

      fname = readConfigVariable(config, 'Surface', 'soil_type')
      if(fname is not None and fname!=""):
        soilTypeVar = processSoilType(fname,dso_s,vars_s,dims_s)
      fname = None
      
      fname = readConfigVariable(config, 'Surface', 'street_type')
      if(fname is not None and fname!=""):
        streetTypeVar = processStreetType(fname,dso_s,vars_s,dims_s)
      fname = None

      fname = readConfigVariable(config, 'Surface', 'surface_fraction')
      if (fname is not None and fname != ""):
        surfaceFractionVar = processSurfaceFraction(fname, dso_s, vars_s,dims_s)
      fname = None
      
      fname = readConfigVariable(config, 'Surface', 'surface_temperature')
      if(fname is not None and fname !=""):
        surfTempTypeVar = processSurfaceTemperature(fname,dso_s,vars_s,dims_s)
      fname = None
      
      fname = readConfigVariable(config, 'Surface', 'ssws')
      if(fname is not None and fname !=""):
        sswsVar = processSurfaceScalarFlux(fname,dso_s,vars_s,dims_s)
      fname = None    
      

    if(vegConfig is not None):
      fname = readConfigVariable(config, 'Vegetation', 'lad')
      if(fname is not None and fname!=""):
        ladVar = processLAD(fname,dso_s,vars_s,dims_s)
      fname = None

      fname = readConfigVariable(config, 'Vegetation', 'vegetation_type')
      if(fname is not None and fname!=""):
        vegetationTypeVar = processVegetationType(fname,dso_s,vars_s,dims_s)
      fname = None

    '''
    Close PIDS_STATIC
    '''
    netcdfWriteAndClose(dso_s, verbose=False)
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
      outputMode = str(input("Existing "+ pidsChemFN +" file found. Overwrite, append, skip or exit? [o/a/s/e] "))
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
      dso_c = netcdfOutputDataset(pidsChemFN, mode="r+")
      dims_c = dso_c.dimensions.keys()
      vars_c = dso_c.variables.keys()
      # Remove dims from vars with a filter
      vars_c = filter(lambda key: key not in dims_c, vars_c)

      if(len(dims_c)>0):
        print("Dimensions: " + ', '.join(str(key) for key in dims_c))
        print("Variables: " + ', '.join(str(key) for key in vars_c))
        print("Warning: using old dimensions for new variables")
      else:
        print("No existing dimensions or variables found.")
    else:
      dso_c = netcdfOutputDataset(pidsChemFN, mode="w")
      dims_c = []
      vars_c = []

    setPIDSGlobalAtrributes(dso_c, globalAttributes)

    # Level-of-detail (LOD) of the emissions
    emissionLod = readConfigVariable(config, 'Chemistry', 'emission_values_lod')

    # Emission category index
    emiCatInds = readConfigVariable(config, 'Chemistry', 'emission_category_index')
    if(emiCatInds is not None and emiCatInds!=""):
      emiCatIndsVar = processEmissionCategoryIndices(emiCatInds, dso_c, vars_c, dims_c)

   # Names of the emission categories
    emiCatName = readConfigVariable(config, 'Chemistry', 'emission_category_name')
    if(emiCatName is not None and emiCatName!=""):
      emiCatName = processEmissionCategoryNames(emiCatName, dso_c, vars_c, dims_c)

    # Emission index per emitted chemical species
    emiInds = readConfigVariable(config, 'Chemistry', 'emission_index')
    if(emiInds is not None and emiInds!=""):
      emiInds = processEmissionIndices(emiInds, dso_c, vars_c, dims_c)

    # Names of the chemical species
    emiName = readConfigVariable(config, 'Chemistry', 'emission_name')
    if (emiName is not None and emiName!=""):
      emiName = processEmissionNames(emiName, dso_c, vars_c, dims_c)

    # If emission_values_lod=1, emission_values are rescaled using emission time factors
    if (emissionLod == '1'):
      emiTimeFactorsFile = readConfigVariable(config, 'Chemistry', 'emission_time_factors')
      emiTimeFactorsLod = readConfigVariable(config, 'Chemistry', 'emission_time_factors_lod')
      if(emiTimeFactorsFile is not None and emiTimeFactorsFile!=""):
        if(emiTimeFactorsLod==""):
          emiTimeFactorsLod = None
        emiTimeFactors = processEmissionTimeFactors(emiTimeFactorsFile, emiTimeFactorsLod,
                                                    dso_c, vars_c, dims_c)

    # If emission_values_lod=2, emissions can vary with time
    elif (emissionLod == '2'):
      emiTimes = readConfigVariable(config, 'Chemistry', 'emission_time')
      if (emiTimes is not None and emiTimes != ""):
        time_dim = createTimeDim(dso_c, emiTimes, dims_c)
      emiTimestamp = readConfigVariable(config, 'Chemistry', 'emission_timestamp')
      if ('True' in emiTimestamp):
        print('create timestamp')
        originTimeStr = readConfigVariable(config, 'Global', 'origin_time')
        timestamp = processEmissionTimestamp(dso_c, originTimeStr, vars_c, dims_c)

    # Emission value map
    emiStr = readConfigVariable(config, 'Chemistry', 'emission_values')
    if (emiStr is not None and emiStr != ""):
      sourceArea = readConfigVariable(config, 'Chemistry', 'source_area')
      if (sourceArea is not None and sourceArea != ""):
        emissionUnit = readConfigVariable(config, 'Chemistry', 'emission_unit')
        aeroEminVar = processEmissionValues(emiStr, sourceArea, emissionUnit, emissionLod,
                                            dso_c, vars_c, dims_c)


    netcdfWriteAndClose(dso_c, verbose=False)
    print("Chemistry output file "+ pidsChemFN +" closed.")

'''
Move on to PIDS_SALSA
'''

if (salsaConfig is not None):
  pidsSalsaFN = "PIDS_SALSA" + nest_id
  print("\n===== Aerosol input file "+ pidsSalsaFN +" =====")
  skipSalsaFile = False
  salsaAppend = False
  if (os.path.isfile(pidsSalsaFN)):
    while True:
      outputMode = str(input("Existing "+ pidsSalsaFN + " file found. " +
                             "Overwrite, append, skip or exit? [o/a/s/e] "))
      if (outputMode in ["e","E"]):
        print("Exit.")
        exit()
      elif (outputMode in ["s","S"]):
        skipSalsaFile = True
        print("Skipping aerosol input.")
        break
      elif (outputMode in ["o","O"]):
        salsaAppend = False
        print("Overwriting an existing " + pidsSalsaFN + " file.")
        break
      elif (outputMode in ["a","A"]):
        salsaAppend=True
        break
      else:
        print("Invalid selection.")
        continue
  else:
    print("No existing " + pidsSalsaFN + " found, creating a new one.")

  if (not skipSalsaFile):
    if (salsaAppend):
      pidsSalsaDS = netcdfOutputDataset(pidsSalsaFN, mode = "r+")
      salsaDims = pidsSalsaDS.dimensions.keys()
      salsaVars = pidsSalsaDS.variables.keys()
      # Remove dims from vars with a filter
      salsaVars = filter(lambda key: key not in salsaDims, salsaVars)
      print("Existing PIDS_SALSA file found, using append/update mode.")

      if (len(salsaDims)>0):
        print("Dimensions: " + ', '.join(str(key) for key in salsaDims))
        print("Variables: " + ', '.join(str(key) for key in salsaVars))
        print("Warning: using old dimensions for new variables")
      else:
        print("No existing dimensions or variables found.")
    else:
      pidsSalsaDS = netcdfOutputDataset(pidsSalsaFN, mode = "w")
      salsaDims = []
      salsaVars = []

    setPIDSGlobalAtrributes(pidsSalsaDS, globalAttributes)

    # Level-of-detail (LOD) of the emissions
    aeroEmissionLod = readConfigVariable(config, 'Aerosols', 'aerosol_emission_values_lod')

    # Emission category index
    emiCatInds = readConfigVariable(config, 'Aerosols', 'aerosol_emission_category_index')
    if(emiCatInds is not None and emiCatInds!=""):
      emiCatIndsVar = processEmissionCategoryIndices(emiCatInds, pidsSalsaDS, salsaVars, salsaDims)

    # Names of the emission categories
    emiCatName = readConfigVariable(config, 'Aerosols', 'aerosol_emission_category_name')
    if(emiCatName is not None and emiCatName!=""):
      emiCatName = processEmissionCategoryNames(emiCatName, pidsSalsaDS, salsaVars, salsaDims)

    # Names of the chemical components
    compNameStr = readConfigVariable(config, 'Aerosols', 'composition_name')
    if (compNameStr is not None and compNameStr != ""):
      compNameVar = processCompositionNames(compNameStr, pidsSalsaDS, salsaVars, salsaDims)

    # Mass fractions of chemical components in aerosol emissions
    emiMassFracsStr = readConfigVariable(config, 'Aerosols', 'emission_mass_fracs')
    if (emiMassFracsStr is not None and emiMassFracsStr != ""):
      emiMassFracsVar = processEmissionMassFracs(emiMassFracsStr, pidsSalsaDS, salsaVars, salsaDims)

    # If aerosol_emission_values_lod=1, aerosol_emission_values are rescaled using emission time
    # factors
    if (aeroEmissionLod == '1'):
      emiTimeFactorsFile = readConfigVariable(config, 'Aerosols', 'aerosol_emission_time_factors')
      emiTimeFactorsLod = readConfigVariable(config, 'Aerosols', 'aerosol_emission_time_factors_lod')
      if (emiTimeFactorsFile is not None and emiTimeFactorsFile!=""):
        if (emiTimeFactorsLod==""):
          emiTimeFactorsLod = None
        emiTimeFactors = processEmissionTimeFactors(emiTimeFactorsFile, emiTimeFactorsLod,
                                                    pidsSalsaDS, salsaVars, salsaDims)

    # If aerosol_emission_values_lod=2, emissions can vary with time. Also, the aerosol size
    # distribution of the emission is given using emission_number_fracs
    elif (aeroEmissionLod == '2'):
      emiTimes = readConfigVariable(config, 'Aerosols', 'aerosol_emission_time')
      if (emiTimes is not None and emiTimes != ""):
        time_dim = createTimeDim(pidsSalsaDS, emiTimes, salsaDims)

      # Number fractions of aerosol size bins in aerosol emissions
      emiNumberFracsStr = readConfigVariable(config, 'Aerosols', 'emission_number_fracs')
      if (emiNumberFracsStr is not None and emiNumberFracsStr != ""):
        aeroEmissionDmid = readConfigVariable(config, 'Aerosols', 'Dmid')
        if (aeroEmissionDmid is not None and aeroEmissionDmid!=""):
          emiNumberFracsVar = processEmissionNumberFracs(emiNumberFracsStr, aeroEmissionDmid,
                                                         pidsSalsaDS, salsaVars, salsaDims)
    # Aerosol emission value map
    aeroEmiStr = readConfigVariable(config, 'Aerosols', 'aerosol_emission_values')
    if (aeroEmiStr is not None and aeroEmiStr != ""):
      aeroSourceArea = readConfigVariable(config, 'Aerosols', 'aerosol_source_area')
      if (aeroSourceArea is not None and aeroSourceArea != ""):
        aeroEmissionUnit = readConfigVariable(config, 'Aerosols', 'aerosol_emission_unit')
        aeroEminVar = processAerosolEmissionValues(aeroEmiStr, aeroSourceArea, aeroEmissionUnit,
                                                   aeroEmissionLod, pidsSalsaDS, salsaVars,
                                                   salsaDims)

    # Close the output file PIDS_SALSA
    netcdfWriteAndClose(pidsSalsaDS, verbose=False)
    print("Aerosol output file "+ pidsSalsaFN +" closed.")

'''
Close PIDS_SALSA
'''
