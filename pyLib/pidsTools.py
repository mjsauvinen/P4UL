import netCDF4 as nc
import sys
import argparse
import numpy as np
import configparser
from itertools import islice, chain, repeat

from netcdfTools import *
from mapTools import readNumpyZTile


'''
Tools for genPIDSInput.py

Author:
Sasu Karttunen
sasu.karttunen@helsinki.fi
Institute for Atmospheric and Earth System Research (INAR) / Physics
University of Helsinki

UPDATE:
- Jani Stromberg: Added building type and surface fraction type/dimension. Fixed a bug in building
  ids.
- Mona Kurppa: Added and modified chemistry and salsa variables
               + changes in the function "map". Does not return a list in Python 3.
'''

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def setPIDSGlobalAtrributes(ds, globalAttributes):
  ''' Set PIDS global attributes to data set file '''
  strAttributeKeys = ['acronym', 'author', 'campaign', 'contact_person', 'creation_time',\
                      'comment', 'Conventions', 'data_content', 'dependencies', 'history'\
                      'keywords', 'license', 'location', 'origin_time', 'references', 'site',\
                      'source', 'title']
  floatAttributeKeys = ['origin_x', 'origin_y', 'origin_z', 'origin_lat', 'origin_lon',\
                        'rotation_angle']
  intAttributeKeys = ['version']
  for key in strAttributeKeys:
    try:
      setattr(ds, key, globalAttributes[key])
    except KeyError:
      if(getattr(ds, key, "")==""):
        setattr(ds, key, "")
  # Mandatory
  setattr(ds, 'Conventions', "CF-1.7")

  for key in floatAttributeKeys:
    try:
      setattr(ds, key, float(globalAttributes[key]))
    except KeyError:
      if(getattr(ds, key, 0.0)==0.0):
        setattr(ds, key, 0.0)

  for key in intAttributeKeys:
    try:
      setattr(ds, key, int(globalAttributes[key]))
    except KeyError:
      if(getattr(ds, key, 0)==0):
        setattr(ds, key, 0)

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def readConfigSection(config, name):
  ''' Prints and returns variables found from given section '''
  try:
    configSection = config._sections[name]
  except KeyError:
    return None

  for attr, val in configSection.items():
    print("{}: {}".format(attr,val))

  return configSection

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def readConfigVariable(config, section, name):
  try:
    var = config.get(section, name)
  except configparser.NoOptionError:
    return None

  return var


#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def parseStringArrayInput(input_str, dtype):
  # Return string type input as 2D numpy array
  # Example input: "2,5,6\n3,6,7\n6,7,9"
  # Example output np.array([[2,5,6],[3,6,7],[6,7,9]])
  rows = input_str.split("\\n")
  if(len(rows[0].split(","))==1):
    rows=rows[1:]
  arr = np.zeros((len(rows),len(rows[0].split(","))),dtype=dtype)
  for i in range(len(rows)):
    items = rows[i].split(",")
    try:
      arr[i,:] = np.array(list(map(dtype, items)))
    except ValueError:
      continue

  return arr

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def parseCharacterArray(input_str, maxstrlen):
  # Return string type input as 2D character array
  # Example input: "abba,cd,acdc"
  # Example output np.array([['a','b','b','a'],['c','d','',''],['a','c','d','c']])
  items = input_str.split(",")
  charr = map(list, items)
  # Fill missing elements with empty char and construct a 2d array
  def pad_array(charr):
    return list(islice(chain(charr, repeat('')),maxstrlen))
  charr = np.array(list(map(pad_array, charr)))
  return charr

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
#                    UNIVERSAL                      #
#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def createXDim(ds, nPx, dPx, dims):
  # Creates a new x-axis unless it already exists
  if('x' not in dims):
    x_dim = createCoordinateAxis(ds, nPx, dPx, 0, 'x', 'f4', 'm', True, False, verbose=False)
    x_dim.long_name = "distance to origin in x-direction"
    dims.append('x')
    return x_dim
  else:
    x_dim = ds.variables['x']
    return x_dim

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def createYDim(ds, nPx, dPx, dims):
  # Creates a new y-axis unless it already exists
  if('y' not in dims):
    y_dim = createCoordinateAxis(ds, nPx, dPx, 1, 'y', 'f4', 'm', True, False, verbose=False)
    y_dim.long_name = "distance to origin in y-direction"
    dims.append('y')
    return y_dim
  else:
    y_dim = ds.variables['y']
    return y_dim

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def createZDim(ds, nPx, dPx, dims):
  # Creates a new z-axis unless it already exists
  if('z' not in dims):
    z_dim = createCoordinateAxis(ds, nPx, dPx, 2, 'z', 'f4', 'm', True, False, verbose=False)
    z_dim.long_name = "height above origin"
    dims.append('z')
    return z_dim
  else:
    z_dim = ds.variables['z']
    return z_dim

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
#                   PIDS_STATIC                     #
#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def createZLADDim(ds, nPx, dPx, dims):
  # Creates a new zlad-axis unless it already exists
  if('zlad' not in dims):
    zlad_dim = createCoordinateAxis(ds, nPx, dPx, 2, 'zlad', 'f4', 'm', True, False, verbose=False)
    zlad_dim.long_name = "height above origin"
    dims.append('zlad')
    return zlad_dim
  else:
    zlad_dim = ds.variables['zlad']
    return zlad_dim

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def createZnsurfaceFractionDim(ds, nPx, dPx, dims):
  # Creates a new nsurface_fraction-axis unless it already exists
  if('nsurface_fraction' not in dims):
    znsurfaceFraction_dim = createCoordinateAxis(ds, nPx, dPx, 2, 'nsurface_fraction', 'i4', '',
                                                 True, False, verbose=False)
    znsurfaceFraction_dim.long_name = "height above origin"
    dims.append('nsurface_fraction')
    return znsurfaceFraction_dim
  else:
    znsurfaceFraction_dim = ds.variables['nsurface_fraction']
    return znsurfaceFraction_dim

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processOrography(fname,ds,vars,dims):
  # Write orography data to given ds
  oroDict = readNumpyZTile(fname,verbose=False)
  oroR = oroDict['R'][::-1,:]
  oroDPx = oroDict['dPx']
  oroNPx = np.shape(oroR)

  if('zt' in vars):
    ds.variables['zt'][:]=oroR
    return ds.variables['zt']
  else:
    # Create new dimensions or check if the orography matches old ones
    x_dim = createXDim(ds, oroNPx, oroDPx, dims)
    y_dim = createYDim(ds, oroNPx, oroDPx, dims)

    oroNCVar = createNetcdfVariable(ds, oroR, 'zt', 0, 'm', 'f4', ('y','x'), False, False,
                                    fill_value=-9999., verbose=False)
    oroNCVar.long_name= "terrain_height"

    return oroNCVar

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processBuildings(fname,ds,vars,dims):
  buildDict = readNumpyZTile(fname,verbose=False)
  buildR = buildDict['R'][::-1,:]
  buildDPx = buildDict['dPx']
  buildNPx = np.shape(buildR)
  buildLOD = len(buildNPx)-1 # 1=2D height field, 2=3D mask

  if(buildLOD==1):
    # Save as a 2D building height array
    if('buildings_2d' in vars):
      ds.variables['buildings_2d'][:]=buildR
      return ds.variables['buildings_2d']
    else:
      x_dim = createXDim(ds, buildNPx, buildDPx, dims)
      y_dim = createYDim(ds, buildNPx, buildDPx, dims)
      buildNCVar = createNetcdfVariable(ds, buildR, 'buildings_2d', 0, 'm', 'f4', ('y','x'), False,
                                        False, fill_value=-9999., verbose=False)
      buildNCVar.long_name = "building_height"
      buildNCVar.lod = int(buildLOD)
      return buildNCVar

  elif(buildLOD==2):
    # Save as a 3D boolean mask array
    #topo=fillTopographyArray(buildR, buildNPx, buildDPx, int)
    # I'm quite not sure why axes swapping has to be done but at least it works
    topo = np.swapaxes(buildR,0,2)
    if('buildings_3d' in vars):
      ds.variables['buildings_3d'][:]=topo
      return ds.variables['buildings_3d']
    else:
      x_dim = createXDim(ds, buildNPx, buildDPx, dims)
      y_dim = createYDim(ds, buildNPx, buildDPx, dims)
      z_dim = createZDim(ds, buildNPx, buildDPx, dims)
      buildNCVar = createNetcdfVariable(ds, topo, 'buildings_3d', 0, 'm', 'b', ('z','y','x'), False,
                                        False, fill_value=-127, verbose=False)
      buildNCVar.long_name = "building_flag"
      buildNCVar.lod = int(buildLOD)
      return buildNCVar
  else:
    raise ValueError("invalid number of dimensions in buildings array: {}".format(buildLOD+1))

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processBuildingIDs(fname,ds,vars,dims):
  buildIDDict = readNumpyZTile(fname,verbose=False)
  buildIDR = buildIDDict['R'][::-1,:]
  buildIDDPx = buildIDDict['dPx']
  buildIDNPx = np.shape(buildIDR)

  if('building_id' in vars):
    ds.variables['building_id'][:]=buildIDR
    return ds.variables['building_id']
  else:
    x_dim = createXDim(ds, buildIDNPx, buildIDDPx, dims)
    y_dim = createYDim(ds, buildIDNPx, buildIDDPx, dims)

    buildIDNCVar = createNetcdfVariable(ds, buildIDR, 'building_id', 0, 'm', 'i4', ('y','x'), False,
                                        False, fill_value=-9999, verbose=False)
    buildIDNCVar.long_name = "building id numbers"

    return buildIDNCVar

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processPavementType(fname,ds,vars,dims):
  pavementTypeDict = readNumpyZTile(fname,verbose=False)
  pavementTypeR = pavementTypeDict['R'][::-1,:]
  pavementTypeDPx = pavementTypeDict['dPx']
  pavementTypeNPx = np.shape(pavementTypeR)

  if('pavement_type' in vars):
    ds.variables['pavement_type'][:]=pavementTypeR
    return ds.variables['pavement_type']
  else:
    x_dim = createXDim(ds, pavementTypeNPx, pavementTypeDPx, dims)
    y_dim = createYDim(ds, pavementTypeNPx, pavementTypeDPx, dims)

    pavementTypeNCVar = createNetcdfVariable(ds, pavementTypeR, 'pavement_type', 0, 'm', 'b',
                                             ('y','x'), False, False, fill_value=-127, verbose=False)
    pavementTypeNCVar.long_name = "pavement type classification"

    return pavementTypeNCVar

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processWaterType(fname,ds,vars,dims):
  waterTypeDict = readNumpyZTile(fname,verbose=False)
  waterTypeR = waterTypeDict['R'][::-1,:]
  waterTypeDPx = waterTypeDict['dPx']
  waterTypeNPx = np.shape(waterTypeR)

  if('water_type' in vars):
    ds.variables['water_type'][:]=waterTypeR
    return ds.variables['water_type']
  else:
    x_dim = createXDim(ds, waterTypeNPx, waterTypeDPx, dims)
    y_dim = createYDim(ds, waterTypeNPx, waterTypeDPx, dims)

    waterTypeNCVar = createNetcdfVariable(ds, waterTypeR, 'water_type', 0, 'm', 'b', ('y','x'),
                                          False, False, fill_value=-127, verbose=False)
    waterTypeNCVar.long_name = "water type classification"

    return waterTypeNCVar

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processSoilType(fname,ds,vars,dims):
  soilTypeDict = readNumpyZTile(fname,verbose=False)
  soilTypeR = soilTypeDict['R'][::-1,:]
  soilTypeDPx = soilTypeDict['dPx']
  soilTypeNPx = np.shape(soilTypeR)

  if('soil_type' in vars):
    ds.variables['soil_type'][:]=soilTypeR
    return ds.variables['soil_type']
  else:
    x_dim = createXDim(ds, soilTypeNPx, soilTypeDPx, dims)
    y_dim = createYDim(ds, soilTypeNPx, soilTypeDPx, dims)

    soilTypeNCVar = createNetcdfVariable(ds, soilTypeR, 'soil_type', 0, 'm', 'b', ('y','x'), False,
                                         False, fill_value=-127, verbose=False)
    soilTypeNCVar.long_name = "soil type classification"

    return soilTypeNCVar

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processStreetType(fname,ds,vars,dims):
  streetTypeDict = readNumpyZTile(fname,verbose=False)
  streetTypeR = streetTypeDict['R'][::-1,:]
  streetTypeDPx = streetTypeDict['dPx']
  streetTypeNPx = np.shape(streetTypeR)

  if('street_type' in vars):
    ds.variables['street_type'][:]=streetTypeR
    return ds.variables['street_type']
  else:
    x_dim = createXDim(ds, streetTypeNPx, streetTypeDPx, dims)
    y_dim = createYDim(ds, streetTypeNPx, streetTypeDPx, dims)

    streetTypeNCVar = createNetcdfVariable(ds, streetTypeR, 'street_type', 0, 'm', 'b', ('y','x'),
                                           False, False, fill_value=-127, verbose=False)
    streetTypeNCVar.long_name = "street type classification"

    return streetTypeNCVar

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processLAD(fname,ds,vars,dims):
  ladDict = readNumpyZTile(fname,verbose=False)
  ladR = ladDict['R']
  ladDPx = ladDict['dPx']
  ladNPx = np.shape(ladR)

  # Same here as in buildings_3d, idk why this has to be done for 3D arrays
  print('1 lad shape = {}'.format(np.shape(ladR)))
  ladR = np.swapaxes(ladR,0,2)
  print('2 lad shape = {}'.format(np.shape(ladR)))

  if('lad' in vars):
    print(' lad is in vars ')
    ds.variables['lad'][:]=ladR
    return ds.variables['lad']
  else:
    print(' lad is NOT in vars ')
    x_dim = createXDim(ds, ladNPx, ladDPx, dims)
    y_dim = createYDim(ds, ladNPx, ladDPx, dims)
    zlad_dim = createZLADDim(ds, ladNPx, ladDPx, dims)

    ladNCVar = createNetcdfVariable(ds, ladR, 'lad', 0, 'm', 'f4', ('zlad','y','x'), False, False,
                                    fill_value=-9999., verbose=False)
    ladNCVar.long_name = "basal area density"

    return ladNCVar

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processVegetationType(fname,ds,vars,dims):
  vegetationTypeDict = readNumpyZTile(fname,verbose=False)
  vegetationTypeR = vegetationTypeDict['R'][::-1,:]
  vegetationTypeDPx = vegetationTypeDict['dPx']
  vegetationTypeNPx = np.shape(vegetationTypeR)

  if('vegetation_type' in vars):
    ds.variables['vegetation_type'][:]=vegetationTypeR
    return ds.variables['vegetation_type']
  else:
    x_dim = createXDim(ds, vegetationTypeNPx, vegetationTypeDPx, dims)
    y_dim = createYDim(ds, vegetationTypeNPx, vegetationTypeDPx, dims)

    vegetationTypeNCVar = createNetcdfVariable(ds, vegetationTypeR, 'vegetation_type', 0, 'm', 'b',
                                               ('y','x'), False, False, fill_value=-127,
                                               verbose=False)
    vegetationTypeNCVar.long_name = "vegetation type classification"

    return vegetationTypeNCVar


# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processBuildingType(fname, ds, vars, dims):
  buildingTypeDict = readNumpyZTile(fname, verbose=False)
  buildingTypeR = buildingTypeDict['R'][::-1, :]
  buildingTypeDPx = buildingTypeDict['dPx']
  buildingTypeNPx = np.shape(buildingTypeR)

  if ('building_type' in vars):
    ds.variables['building_type'][:] = buildingTypeR
    return ds.variables['building_type']
  else:
    x_dim = createXDim(ds, buildingTypeNPx, buildingTypeDPx, dims)
    y_dim = createYDim(ds, buildingTypeNPx, buildingTypeDPx, dims)

    buildingTypeNCVar = createNetcdfVariable(ds, buildingTypeR, 'building_type', 0, 'm', 'b',
                                             ('y', 'x'), False, False, fill_value=-127,
                                             verbose=False)
    buildingTypeNCVar.long_name = "building type classification"

    return buildingTypeNCVar


# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


def processSurfaceFraction(fname, ds, vars, dims):
  surfaceFractionDict = readNumpyZTile(fname, verbose=False)
  surfaceFractionR = surfaceFractionDict['R'][::-1, :, :]
  surfaceFractionDPx = surfaceFractionDict['dPx']
  surfaceFractionNPx = np.shape(surfaceFractionR)

  # Same here as in buildings_3d, idk why this has to be done for 3D arrays
  surfaceFractionR = np.swapaxes(surfaceFractionR, 0, 2)
  surfaceFractionR = np.swapaxes(surfaceFractionR, 2, 1)

  if ('surface_fraction' in vars):
    ds.variables['surface_fraction'][:] = surfaceFractionR
    return ds.variables['surface_fraction']
  else:
    x_dim = createXDim(ds, surfaceFractionNPx, surfaceFractionDPx, dims)
    y_dim = createYDim(ds, surfaceFractionNPx, surfaceFractionDPx, dims)
    znsurface_fraction_dim = createZnsurfaceFractionDim(ds, surfaceFractionNPx, surfaceFractionDPx,
                                                        dims)

    surfaceFractionNCVar = createNetcdfVariable(ds, surfaceFractionR, 'surface_fraction', 0, 'm',
                                                'f4', ('nsurface_fraction', 'y', 'x'), False, False,
                                                fill_value=-9999., verbose=False)
    surfaceFractionNCVar.long_name = "surface fraction"

    return surfaceFractionNCVar


#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
#           Both PIDS_CHEM and PIDS_SALSA           #
#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def createNcatDim(ds, ncat, dims):
  # Creates a new ncat dim unless it already exists

  if ('ncat' not in dims):
    ncat_dim = createNetcdfVariable(ds, ncat, 'ncat', len(ncat), '', 'i4', ('ncat',),
                                    parameter=True, verbose=False)
    ncat_dim.long_name = "number of emission categories"
    dims.append('ncat')
    return ncat_dim
  else:
    return ds.variables['ncat']

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def createStrlenDim(ds, strlen, dims):
  # Creates a new strlen dim unless it already exists

  if ('strlen' not in dims):
    strlen_dim = createNetcdfVariable(ds, strlen, 'strlen', len(strlen), '', 'i4', ('strlen',),
                                      parameter=True, verbose=False)
    dims.append('strlen')
    return strlen_dim
  else:
    return ds.variables['strlen']

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def createNhoursyearDim(ds, nhoursyear, dims):
  # Creates a new nhoursyear dim unless it already exists

  if ('nhoursyear' not in dims):
    nhoursyear_dim = createNetcdfVariable(ds, nhoursyear, 'nhoursyear', len(nhoursyear), '', 'i4',
                                          ('nhoursyear',), parameter=True, verbose=False)
    dims.append('nhoursyear')
    return nhoursyear_dim
  else:
    return ds.variables['nhoursyear']

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def createNmonthdayhourDim(ds, nmonthdayhour, dims):
  # Creates a new nmonthdayhour dim unless it already exists

  if ('nmonthdayhour' not in dims):
    nmonthdayhour_dim = createNetcdfVariable(ds, nmonthdayhour, 'nmonthdayhour', len(nmonthdayhour),
                                             '', 'i4', ('nmonthdayhour',), parameter=True,
                                             verbose=False)
    dims.append('nmonthdayhour')
    return nmonthdayhour_dim
  else:
    return ds.variables['nmonthdayhour']

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def createTimeDim(ds, time, dims):
  # Creates a new ncat dim unless it already exists

  if ('time' not in dims):
    time = list(map(float, time.split(",")))
    time_dim = createNetcdfVariable(ds, time, 'time', len(time), 's', 'f4', ('time',),
                                    parameter=True, verbose=False)
    time_dim.long_name = "seconds since the beginning of the simulation"
    dims.append('time')
    return time_dim
  else:
    return ds.variables['time']

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processEmissionTimestamp(ds, originTime, vars, dims):
  # Creates a new ncat dim unless it already exists

  from datetime import datetime, timedelta

  if ('field_len' not in dims):
    field_len = 64 # some default value, I don't know why (Mona)
    ds.createDimension('field_len', field_len)

  if ('timestamp' not in vars):
    time = ds.variables['time'][:]
    dtOriginTime = datetime.strptime(originTime+"00", "%Y-%m-%d %H:%M:%S %z")
    timestamp = []
    emptyStr = " " * field_len
    for t in range(len(time)):
      dt = timedelta(seconds = int(time[t]))
      dtStr = dtOriginTime + dt
      timestampStr = dtStr.strftime("%Y-%m-%d %H:%M:%S %z")[:-2].ljust(field_len)
      timestamp.append(timestampStr)

    timestampVar = ds.createVariable('timestamp', 'S1', ('time','field_len',))
    timestampVar[:] = list(map(lambda x : list(x), timestamp))
    timestampVar.long_name = "timestamps since the beginning of the simulation"
    vars.append('timestamp')
    return timestampVar

  else:
    return ds.variables['timestamp']

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processEmissionCategoryIndices(emiCatInds, ds, vars, dims):
  # In my opinion ncat is completely redundant parameter, emission_category_index should
  # be a coordinate variable with a length of ncat instead. Current setup in PALM doesn't
  # make any sense (Sasu)

  try:
    emiCatInds = list(map(int, emiCatInds.split(",")))
  except TypeError:
    print("Error: invalid value for emission_category_index in configuration file, expected a"+
          " comma-delimited list")
    exit(1)

  if ('emission_category_index' in vars):
    ds.variables['emission_category_index'][:] = emiCatInds
    return ds.variables['emission_category_index']

  else:
    ncat_dim = createNcatDim(ds, np.arange(1, len(emiCatInds)+1, 1), dims)
    emiCatIndsVar = createNetcdfVariable(ds, emiCatInds, 'emission_category_index', 0, '', 'i1',
                                         ('ncat',), parameter=False, verbose=False)
    emiCatIndsVar.long_name = "emission category index"
    emiCatIndsVar.standard_name = 'emission_cat_index'
    return emiCatIndsVar

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processEmissionCategoryNames(emiCatName, ds, vars, dims):
  # Creates emission_category_name unless it already exists

  try:
    # Max string length is chosen quite arbitralily here
    maxstrlen = 25
    emiCatName = parseCharacterArray(emiCatName, maxstrlen)
  except TypeError:
    print("Error: invalid value for emission_category_name in configuration file, expected a " +
          "comma delimited array")
    exit(1)

  if ('emission_category_name' in vars):
    ds.variables['emission_category_name'][:] = emiCatName
    return ds.variables['emission_category_name']

  else:
    ncat_dim = createNcatDim(ds, np.arange(1, np.shape(emiCatName)[0]+1, 1), dims)
    strlen_dim = createStrlenDim(ds, np.arange(1, maxstrlen+1, 1), dims)
    emiCatNameVar = createNetcdfVariable(ds, np.array(emiCatName), 'emission_category_name', 0, '',
                                         'S1', ('ncat','strlen',), parameter=False, verbose=False)
    emiCatNameVar.long_name = 'emission category name'
    return emiCatNameVar

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processEmissionTimeFactors(fname, lod, ds, vars, dims):
  # Create emission_time_factors unless it already exists

  # Try to read 'emission_time_factors' array from a given input file
  constant_factor = False
  try:
    emiTimeFactors = np.genfromtxt(fname, delimiter=",")
    ncat = np.shape(emiTimeFactors)[0]
    ncat_dim = createNcatDim(ds, np.arange(1, np.shape(emiTimeFactors)[0]+1, 1), dims)
  except:
    print("Cannot read emission_time_factors data from file \"{}\"".format(fname) +
          ". Use a constant factor emission_time_factors=1.0 and ncat=1.")
    factor = 1.0
    constant_factor = True
    if ('ncat' in dims):
      ncat = len(ds.variables['ncat'][:])
    else:
      ncat = 1
      ncat_dim = createNcatDim(ds, [ncat], dims)

  # Check the level of detail
  if (lod is None or lod=='2'):

    # Using a weighting factors given separately for each hour of the year
    nhoursyear = np.arange(1, 8760+1, 1) # 24*365
    if constant_factor:
      emiTimeFactors = np.zeros([ncat, len(nhoursyear)], dtype=float) + factor

    if (np.shape(emiTimeFactors)[-1]!=8760):
      raise ValueError("emission_time_factors data must contain exactly 8760 datapoints for "+
                       "every emission category when emission_time_factors_lod = 2")

    if ('emission_time_factors' in vars):
      if (np.shape(ds.variables['emission_time_factors'])[-1]!=8760):
        raise ValueError("The dimensions of the existing emission_time_factors does not match "+
                         "with the new emission_time_factors data")
      ds.variables['emission_time_factors'][:] = emiTimeFactors
      return ds.variables['emission_time_factors']

    else:
      nhoursyear_dim = createNhoursyearDim(ds, nhoursyear, dims)
      emiTimeFactorsVar = createNetcdfVariable(ds, emiTimeFactors, 'emission_time_factors', 0, '',
                                               'f4', ('ncat','nhoursyear',), parameter=False,
                                               verbose=False)
      emiTimeFactorsVar.long_name = "emission time scaling factors"
      emiTimeFactorsVar.lod = 2
      return emiTimeFactorsVar

  elif (lod=='1'):
    # Using a weighting factors based on the month, day of week and hour
    nmonthdayhour = np.arange(1, 91+1, 1) # see the documentation
    if constant_factor:
      emiTimeFactors = np.zeros([ncat, len(monthdayhour)], dtype=float) + factor

    if (np.shape(emiTimeFactors)[-1]!=91):
      raise ValueError("emission_time_factors data must contain exactly 90 datapoints for "+
                       "every emission category when emission_time_factors_lod = 1")

    if ('emission_time_factors' in vars):
      if (np.shape(ds.variables['emission_time_factors'])[-1]!=90):
        raise ValueError("The dimensions of the existing emission_time_factors does not match "+
                         "with the new emission_time_factors data")
      ds.variables['emission_time_factors'][:] = emiTimeFactors
      return ds.variables['emission_time_factors']

    else:
      nhoursyear_dim = createNmonthdayhourDim(ds, nmonthdayhour, dims)
      emiTimeFactorsVar = createNetcdfVariable(ds, emiTimeFactors, 'emission_time_factors', 0, '',
                                               'f4', ('ncat','nmonthdayhour',), parameter=False,
                                               verbose=False)
      emiTimeFactorsVar.long_name = "emission time scaling factors"
      emiTimeFactorsVar.lod = 1
      return emiTimeFactorsVar

  else:
    raise ValueError("invalid value for emission_time_factors_lod: {}".format(lod))

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
#                    PIDS_CHEM                      #
#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def createZChemDim(ds, z, dims):
  # Creates a new z dim unless it already exists

  if ('z' not in dims):
    zchem_dim = createNetcdfVariable(ds, z, 'z', len(z), '', 'i4', ('z',), parameter=True,
                                     verbose=False)
    zchem_dim.long_name = "distance to origin in z-direction"
    dims.append('z')
    return zchem_dim
  else:
    return ds.variables['z']

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def createNspeciesDim(ds, nspecies, dims):
  # Creates a new nspecies dim unless it already exists

  if ('nspecies' not in dims):
    nspecies_dim = createNetcdfVariable(ds, nspecies, 'nspecies', len(nspecies), '', 'i4',
                                        ('nspecies',), parameter=True, verbose=False)
    nspecies_dim.long_name = "number of emission species"
    dims.append('nspecies')
    return nspecies_dim
  else:
    return ds.variables['nspecies']

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processEmissionIndices(emiInds, ds, vars, dims):
   # Creates emission_index unless it already exists. Again, nspecies is redundant (Sasu)

  try:
    emiInds = list(map(int, emiInds.split(",")))
  except TypeError:
    print("Error: invalid value for emission_index in configuration file, expected a comma " +
          "delimited list")
    exit(1)

  if ('emission_index' in vars):
    ds.variables['emission_index'][:] = emiInds
    return ds.variables['emission_index']

  else:
    nspecies_dim = createNspeciesDim(ds, np.arange(1, len(emiInds)+1, 1), dims)
    emiIndsVar = createNetcdfVariable(ds, emiInds, 'emission_index', 0, '', 'u2', ('nspecies',),
                                      parameter=False, verbose=False)
    emiIndsVar.long_name = "emission species index"
    return emiIndsVar

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processEmissionNames(emiName, ds, vars, dims):
  # Creates emission_name unless it already exists

  try:
    # Max string length is chosen quite arbitralily here
    maxstrlen = 25
    emiName = parseCharacterArray(emiName, maxstrlen)
  except TypeError:
    print("Error: invalid value for emission_name in configuration file, expected a "+
          "comma delimited array")
    exit(1)

  if ('emission_name' in vars):
    ds.variables['emission_name'][:] = emiName
    return ds.variables['emission_name']

  else:
    strlen_dim = createStrlenDim(ds, np.arange(1, maxstrlen+1, 1), dims)
    emiNameVar = createNetcdfVariable(ds, np.array(emiName), 'emission_name', 0, '', 'S1',
                                         ('nspecies','strlen',), parameter=False, verbose=False)
    emiNameVar.long_name = "emission species name"
    return emiNameVar

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processEmissionValues(emiStr, fnameSource, unit, lod, ds, vars, dims):
  # Creates aerosol_emission_values unless it already exists

  # Read in a source map: contains emission category (ncat) values
  sourceDict = readNumpyZTile(fnameSource, verbose=False)
  sourceR = sourceDict['R'][::-1,:,:]
  sourceDPx = sourceDict['dPx']
  sourceNPx = np.shape(sourceR)

  try:
    emiVals = parseStringArrayInput(emiStr, float) # emission value (per time and) per category
  except TypeError:
    print("Error: invalid value for emission_values in configuration file, expected a " +
          "comma-delimited list")
    exit(1)

  nspecies = len(ds.variables['nspecies'][:]) # number of chemical species
  ncat = len(ds.variables['ncat'][:]) # number of emission categories

  if (not 'emission_values' in vars):
    x_dim = createXDim(ds, sourceNPx, sourceDPx, dims)
    y_dim = createYDim(ds, sourceNPx, sourceDPx, dims)
    z_dim = createZChemDim(ds, np.array([1]), dims)

  if (lod == '1'): # emission_values(z,y,x,nspecies,ncat) where z=1
    emission_values = np.zeros([1, sourceNPx[0], sourceNPx[1], nspecies, ncat],
                               dtype=float) - 9999.
    emiVals = np.squeeze(emiVals)
    for ispec in range(nspecies):
      for n in range(ncat):
        emission_values[0,sourceR[:,:,n]==1,ispec,n] = emiVals[ispec,n]

    if ('emission_values' in vars):
      ds.variables['emission_values'][:] = emission_values
      return ds.variables['emission_values']

    else:
      emiValVar = createNetcdfVariable(ds, emission_values, 'emission_values', 0, unit, 'f4',
                                       ('z','y','x','nspecies','ncat',), False, False,
                                        fill_value=-9999., verbose=False)
      emiValVar.long_name= 'emission values'
      emiValVar.lod = 1
      return emiValVar

  elif (lod is None or lod == '2'): # emission_values(time,z,y,x,nspecies) where z=1
    times = ds.variables['time']
    emission_values = np.zeros([len(times), 1, sourceNPx[0], sourceNPx[1], nspecies], dtype=float)
    for n in range(ncat):
      for t in range(len(times)):
        for ispec in range(nspecies):
          emission_values[t,0,sourceR[:,:,n]==1,ispec] += emiVals[ispec,t]
    emission_values[emission_values < 1e-30] = -9999. # fill value

    if ('emission_values' in vars):
      ds.variables['emission_values'][:] = emission_values
      return ds.variables['emission_values']

    else:
      emiValVar = createNetcdfVariable(ds, emission_values, 'emission_values', 0, unit, 'f4',
                                       ('time','z','y','x','nspecies',), False, False,
                                       fill_value=-9999., verbose=False)
      emiValVar.long_name= 'emission values'
      emiValVar.lod = 2
      return emiValVar

  else:
    raise ValueError("invalid value for emission_lod: {}".format(lod))

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
#                    PIDS_SALSA                     #
#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def createCompositionIndexDim(ds, composition_index, dims):
  # Creates a new composition_index dim unless it already exists

  if ('composition_index' not in dims):
    composition_index_dim = createNetcdfVariable(ds, composition_index, 'composition_index',
                                                 len(composition_index), '', 'i4',
                                                 ('composition_index',), parameter=True, verbose=False)
    composition_index_dim.long_name = "aerosol composition index"
    dims.append('composition_index')
    return composition_index_dim
  else:
    return ds.variables['composition_index']

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processCompositionNames(compName, ds, vars, dims):
  # Creates a new composition_name variable unless it already exists
  try:
    # Max string length is chosen quite arbitralily here
    maxstrlen = 25
    compName = parseCharacterArray(compName, maxstrlen)
  except TypeError:
    print("Error: invalid value for composition_name in configuration file, expected a "+
          "comma delimited array")
    exit(1)

  if ('composition_name' in vars):
    ds.variables['composition_name'][:] = compName
    return ds.variables['composition_name']

  else:
    composition_index =  np.arange(1, np.shape(np.array(compName))[0]+1, 1)
    composition_index_dim = createCompositionIndexDim(ds, composition_index, dims)
    compNameVar = ds.createVariable('composition_name', 'S1', ('composition_index','strlen',))
    compNameVar[:] = list(map(lambda x : list(x), compName))
    compNameVar.long_name = 'aerosol composition name'
    return compNameVar

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processAerosolEmissionValues(emiStr, fnameSource, unit, lod, ds, vars, dims):
  # Creates aerosol_emission_values unless it already exists

  # Read in a source map: contains emission category (ncat) values
  sourceDict = readNumpyZTile(fnameSource, verbose=False)
  sourceR = sourceDict['R'][::-1,:,:]
  sourceDPx = sourceDict['dPx']
  sourceNPx = np.shape(sourceR)

  try:
    emiVals = parseStringArrayInput(emiStr, float) # emission value (per time and) per category
  except TypeError:
    print("Error: invalid value for aerosol_emission_values in configuration file, expected a " +
          "comma-delimited list")
    exit(1)

  if ('aerosol_emission_values' not in vars):
    x_dim = createXDim(ds, sourceNPx, sourceDPx, dims)
    y_dim = createYDim(ds, sourceNPx, sourceDPx, dims)
    ncat = np.shape(emiVals)[0]
    ncat_dim = createNcatDim(ds, np.arange(1, ncat+1, 1), dims)

  if (lod == '1'):
    aerosol_emission_values = np.zeros([sourceNPx[0], sourceNPx[1], ncat], dtype=float) - 9999.
    emiVals = np.squeeze(emiVals)
    for n in range(ncat):
      aerosol_emission_values[sourceR[:,:,n]==1,n] = emiVals[n]

    if ('aerosol_emission_values' in vars):
      ds.variables['aerosol_emission_values'][:] = aerosol_emission_values
      return ds.variables['aerosol_emission_values']

    else:
      emiValVar = createNetcdfVariable(ds, aerosol_emission_values, 'aerosol_emission_values', 0,
                                       unit, 'f4', ('y','x','ncat',), False, False, fill_value=-9999.,
                                       verbose=False)
      emiValVar.long_name= 'aerosol emission values'
      emiValVar.lod = 1
      return emiValVar

  elif (lod is None or lod == '2'):
    times = ds.variables['time']
    aerosol_emission_values = np.zeros([len(times), sourceNPx[0], sourceNPx[1], ncat],
                                       dtype=float) - 9999.
    for t in range(len(times)):
      for n in range(ncat):
        aerosol_emission_values[t,sourceR[:,:,n]==1,n] = emiVals[n,t]

    if ('aerosol_emission_values' in vars):
      ds.variables['aerosol_emission_values'][:] = aerosol_emission_values
      return ds.variables['aerosol_emission_values']

    else:
      emiValVar = createNetcdfVariable(ds, aerosol_emission_values, 'aerosol_emission_values', 0,
                                       unit, 'f4', ('time','y','x','ncat',), False, False,
                                       fill_value=-9999., verbose=False)
      emiValVar.long_name= 'aerosol emission values'
      emiValVar.lod = 2
      return emiValVar

  else:
    raise ValueError("invalid value for aerosol_emission_lod: {}".format(lod))

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processEmissionMassFracs(emiMassFracsD, ds, vars, dims):
  # Creates emission_mass_fracs unless it already exists

  ncat = len(ds.variables['ncat'])
  try:
    emiMassFracs = parseStringArrayInput(emiMassFracsD, float)
  except TypeError:
    print("Error: invalid value for emission_mass_fracs in configuration file, expected a " +
          "newline and comma delimited matrix")
    exit(1)

  if (np.shape(emiMassFracs)[0] != ncat):
    raise ValueError("Not correct dimensions of emission_mass_fracs(ncat,composition_index)")
  if (np.shape(emiMassFracs)[1] != len(ds.variables['composition_index'])):
    raise ValueError("Not correct dimensions of emission_mass_fracs(ncat,composition_index)")

  if ('emission_mass_fracs' in vars):
    ds.variables['emission_mass_fracs'][:] = emiMassFracs
    return ds.variables['emission_mass_fracs']

  else:
    emiMassFracsVar = createNetcdfVariable(ds, emiMassFracs, 'emission_mass_fracs', 0, '', 'f4',
                                           ('ncat','composition_index',), False, False,
                                           fill_value=-9999., verbose=False)
    emiMassFracsVar.long_name = "mass fractions of chemical components in aerosol emissions"
    emiMassFracsVar.units = ""
    return emiMassFracsVar

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processEmissionNumberFracs(emiNumberFracsD, emiDmid, ds, vars, dims):
  # Creates emission_numbs_fracs unless it already exists

  ncat = len(ds.variables['ncat'])
  if ('Dmid' in dims):
    emiDmids = ds.variables['Dmid']
  else:
    try:
      emiDmids = list(map(float, emiDmid.split(",")))
      dmid_Dim = createNetcdfVariable(ds, emiDmids, 'Dmid', len(emiDmids), 'm', 'f4', ('Dmid',),
                                      parameter=True, verbose=False)
      dims.append('Dmid')
    except TypeError:
      print("Error: invalid value for aerosol_emission_dmid in configuration file, expected a " +
            "comma-delimited list")
      exit(1)
  nbins = len(emiDmids)

  try:
    emiNumberFracs = parseStringArrayInput(emiNumberFracsD, float)
  except TypeError:
    print("Error: invalid value for composition_aerosol in configuration file, expected a " +
        "newline and comma delimited matrix")
    exit(1)

  if (np.shape(emiNumberFracs)[0] != ncat):
    raise ValueError("Incorrect 0 dimension in emission_number_fracs(ncat,Dmid)")
  if (np.shape(emiNumberFracs)[1] != nbins):
    raise ValueError("Incorrect 1 dimension in emission_number_fracs(ncat,Dmid)")

  if ('emission_number_fracs' in vars):
    ds.variables['emission_number_fracs'][:] = emiNumberFracs
    return ds.variables['emission_number_fracs']

  else:
    emiNumberFracsVar = createNetcdfVariable(ds, emiNumberFracs, 'emission_number_fracs', 0, '',
                                             'f4', ('ncat','Dmid',), False, False,
                                             fill_value=-9999., verbose=False)
    emiNumberFracsVar.long_name = "number fractions of aerosol size bins in aerosol emissions"
    emiNumberFracsVar.units = ""
    return emiNumberFracsVar
