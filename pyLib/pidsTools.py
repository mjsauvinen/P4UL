import netCDF4 as nc
import sys
import argparse
import numpy as np
import ConfigParser

from netcdfTools import *
from mapTools import readNumpyZTile


'''
Author:
Sasu Karttunen
sasu.karttunen@helsinki.fi
Institute for Atmospheric and Earth System Research (INAR) / Physics
University of Helsinki
'''

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def setPIDSGlobalAtrributes(ds, globalAttributes):
  ''' Set PIDS global attributes to data set file '''
  strAttributeKeys = ['Conventions', 'palm_version', 'title', 'acronym', 'campaign',\
                   'institution', 'author', 'contact_person', 'licence', 'history',\
                   'keywords', 'references', 'comment', 'data_content', 'source',\
                   'dependencies', 'location', 'site']
  floatAttributeKeys = ['origin_x', 'origin_y', 'origin_z', 'origin_lat', 'origin_lon',\
                        'rotation_angle', 'origin_time']
  for key in strAttributeKeys:
    try:
      setattr(ds, key, globalAttributes[key])
    except KeyError:
      if(getattr(ds, key, "")==""):
        setattr(ds, key, "")

  for key in floatAttributeKeys:
    try:
      setattr(ds, key, float(globalAttributes[key]))
    except KeyError:
      if(getattr(ds, key, 0.0)==0.0):
        setattr(ds, key, 0.0)

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def readConfigSection(config, name):
  ''' Prints and returns variables found from given section '''
  try:
    configSection = config._sections[name]
  except KeyError:
    return None

  del configSection['__name__']
  for attr, val in configSection.iteritems():
    print("{}: {}".format(attr,val))

  return configSection

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def readConfigVariable(config, section, name):
  try:
    var = config.get(section, name)
  except ConfigParser.NoOptionError:
    return None

  return var

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

def processOrography(fname,ds,vars,dims):
  # Write orography data to given ds
  oroDict = readNumpyZTile(fname,verbose=False)
  oroR = oroDict['R']
  oroDPx = oroDict['dPx']
  oroNPx = np.shape(oroR)

  if('zt' in vars):
    ds.variables['zt'][:]=oroR
    return ds.variables['zt']
  else:
    # Create new dimensions or check if the orography matches old ones
    x_dim = createXDim(ds, oroNPx, oroDPx, dims)
    y_dim = createYDim(ds, oroNPx, oroDPx, dims)

    oroNCVar = createNetcdfVariable(ds, oroR, 'zt', 0, 'm', 'f4', ('y','x'), False, False, fill_value=-9999.9, verbose=False)
    oroNCVar.long_name= "terrain_height"

    return oroNCVar

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processBuildings(fname,ds,vars,dims):
  buildDict = readNumpyZTile(fname,verbose=False)
  buildR = buildDict['R']
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
      buildNCVar = createNetcdfVariable(ds, buildR, 'buildings_2d', 0, 'm', 'f4', ('y','x'), False, False, fill_value=-9999.9, verbose=False)
      buildNCVar.long_name= "building_height"
      buildNCVar.lod = int(buildLOD)
      return buildNCVar

  elif(buildLOD==2):
    # Save as a 3D boolean mask array

    # I'm quite not sure why axes swapping has to be done but at least it works
    buildR = np.swapaxes(buildR,0,2)
    if('buildings_3d' in vars):
      ds.variables['buildings_3d'][:]=buildR
      return ds.variables['buildings_3d']
    else:
      x_dim = createXDim(ds, buildNPx, buildDPx, dims)
      y_dim = createYDim(ds, buildNPx, buildDPx, dims)
      z_dim = createZDim(ds, buildNPx, buildDPx, dims)
      buildNCVar = createNetcdfVariable(ds, buildR, 'buildings_3d', 0, 'm', 'b', ('z','y','x'), False, False, fill_value=-127, verbose=False)
      buildNCVar.long_name= "building_flag"
      buildNCVar.lod = int(buildLOD)
      return buildNCVar
  else:
    raise ValueError("invalid number of dimensions in buildings array: {}".format(buildLOD+1))

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processBuildingIDs(fname,ds,vars,dims):
  buildIDDict = readNumpyZTile(fname,verbose=False)
  buildIDR = buildIDDict['R']
  buildIDDPx = buildIDDict['dPx']
  buildIDNPx = np.shape(buildIDR)

  if('building_id' in vars):
    ds.variables['building_id'][:]=buildIDR
    return ds.variables['building_id']
  else:
    x_dim = createXDim(ds, buildIDNPx, buildIDDPx, dims)
    y_dim = createYDim(ds, buildIDNPx, buildIDDPx, dims)

    buildIDNCVar = createNetcdfVariable(ds, buildIDR, 'building_id', 0, 'm', 'i4', ('y','x'), False, False, fill_value=-9999, verbose=False)
    buildIDNCVar.long_name= "building id numbers"

    return buildIDNCVar

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processPavementType(fname,ds,vars,dims):
  pavementTypeDict = readNumpyZTile(fname,verbose=False)
  pavementTypeR = pavementTypeDict['R']
  pavementTypeDPx = pavementTypeDict['dPx']
  pavementTypeNPx = np.shape(pavementTypeR)

  if('pavement_type' in vars):
    ds.variables['pavement_type'][:]=pavementTypeR
    return ds.variables['pavement_type']
  else:
    x_dim = createXDim(ds, pavementTypeNPx, pavementTypeDPx, dims)
    y_dim = createYDim(ds, pavementTypeNPx, pavementTypeDPx, dims)

    pavementTypeNCVar = createNetcdfVariable(ds, pavementTypeR, 'pavement_type', 0, 'm', 'i4', ('y','x'), False, False, fill_value=-9999, verbose=False)
    pavementTypeNCVar.long_name= "pavement type classification"

    return pavementTypeNCVar

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processLAD(fname,ds,vars,dims):
  ladDict = readNumpyZTile(fname,verbose=False)
  ladR = ladDict['R']
  ladDPx = ladDict['dPx']
  ladNPx = np.shape(ladR)

  # Same here as in buildings_3d, idk why this has to be done for 3D arrays
  ladR = np.swapaxes(ladR,0,2)

  if('lad' in vars):
    ds.variables['lad'][:]=ladR
    return ds.variables['lad']
  else:
    x_dim = createXDim(ds, ladNPx, ladDPx, dims)
    y_dim = createYDim(ds, ladNPx, ladDPx, dims)
    zlad_dim = createZLADDim(ds, ladNPx, ladDPx, dims)

    ladNCVar = createNetcdfVariable(ds, ladR, 'lad', 0, 'm', 'f4', ('z','y','x'), False, False, fill_value=-9999.9, verbose=False)
    ladNCVar.long_name= "leaf area density"

    return ladNCVar

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def processVegetationType(fname,ds,vars,dims):
  vegetationTypeDict = readNumpyZTile(fname,verbose=False)
  vegetationTypeR = vegetationTypeDict['R']
  vegetationTypeDPx = vegetationTypeDict['dPx']
  vegetationTypeNPx = np.shape(vegetationTypeR)

  if('vegetation_type' in vars):
    ds.variables['vegetation_type'][:]=vegetationTypeR
    return ds.variables['vegetation_type']
  else:
    x_dim = createXDim(ds, vegetationTypeNPx, vegetationTypeDPx, dims)
    y_dim = createYDim(ds, vegetationTypeNPx, vegetationTypeDPx, dims)

    vegetationTypeNCVar = createNetcdfVariable(ds, vegetationTypeR, 'vegetation_type', 0, 'm', 'i4', ('y','x'), False, False, fill_value=-9999, verbose=False)
    vegetationTypeNCVar.long_name= "vegetation type classification"

    return vegetationTypeNCVar

#=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
