#!/usr/bin/env python3

import netCDF4 as nc
import sys
import argparse
import numpy as np
from utilities import partialMatchFromList, popKeyFromDict

debug = True

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def asciiEncode(uList, uStr):
  n = len(uList)
  if(n > 0):
    uList = list(uList)  # This might be a tuple coming in
    for i in range(len(uList)):
      if isinstance(uList[i], bytes): uList[i] = uList[i].decode()
  else:
    print(' Dictionary {} has zero length. Exiting ...'.format(uStr))
    sys.exit(1)

  return uList

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def netcdfDataset(filename, verbose=True):
  # Create Dataset
  ds = nc.Dataset(filename)

  # Generate a list of variables and independent variables contained in the file.
  varList = asciiEncode(ds.variables.keys(), 'Variables')
  dimList = asciiEncode(ds.dimensions.keys(), 'Dimensions')

  if(verbose):
    print(' Variable List : {} '.format(varList))
    print(' Dimension List : {} '.format(dimList))

  return ds, varList, dimList

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def netcdfDataset2(filename, verbose=True):
  # Create Dataset
  ds = nc.Dataset(filename)
  
  varList = ds.variables.keys()
  dimList = ds.dimensions.keys()
  # Generate a dictionary that contains an individual listing of dimensions for each variable
  vD = dict() # Variable dimensions
  uD = dict() # All units
  for vn in varList:
    if( vn not in dimList ):
      vD[vn] = ds.variables[vn].dimensions
      uD[vn] = ds.variables[vn].units
    else:
      uD[vn] = ds.variables[vn].units
      
  return ds, vD, uD

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def netcdfOutputDataset(filename, mode='w'):
  if( isinstance( filename, bytes ) ):
    filename = filename.decode()
  dso = nc.Dataset(filename, mode, format='NETCDF4')
  return dso

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def netcdfWriteAndClose(dso, verbose=True):
  if(verbose):
    print('Writing of output data  .... ')
  dso.close()
  if(verbose):
    print(' ... done. File closed.')

  dso = None

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def read1DVariableFromDataset( dimStr, varStr, ds, iLOff=0, iROff=0, cl=1):
  # iLOff: left offset
  # iROff: right offset
  # cl   : coarsening level
  if(varStr in ds.variables.keys()):
    vs = ds.variables[varStr]
    dimList = vs.dimensions   # return a list of variable dimensions ('time', 'x', 'y', etc.)
    print(' dimList = {} '.format( dimList ))
    vdim = partialMatchFromList( dimStr, dimList )
    print(' Reading variable {} ... '.format(vdim))
    try:     
      var = ds.variables[vdim][:]
    except:  
      print(' Cannot read the array of variable: {}.'.format(varStr))
      sys.exit(1)
    
    if(iROff == 0 or (iROff is None) ):
      var = var[(0 + iLOff):]
    else:
      var = var[(0 + iLOff):-abs(iROff)]
    
    if( cl > 1 ): var = var[::cl]
    
  else:
    print(' Variable {} not in list {}.'.format(varStr, ds.variables.keys()))
    sys.exit(1)

  return var, np.shape(var)

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


def readVariableFromDataset(varStr, ds, cl=1 ):
  if( varStr in ds.variables.keys() ):

    vdims = asciiEncode(ds.variables[varStr].dimensions, ' Variable dimensions ')
    var   = ds.variables[varStr][:]

    if( cl>1 ):
      var = coarsenVariable( var, vdims, cl )

    # Load the independent variables and wrap them into a dict
    dDict = dict()
    for dname in vdims:
      dData = ds.variables[dname][:]
      if( ('time' in dname) or ('t+' in dname) ): 
        dDict[dname] = dData
      else:
        dDict[dname] = dData[::cl]
      dData = None

  else:
    sys.exit(' Variable {} not in list {}.'.format(varStr, ds.variables.keys()))

  return var, dDict

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def coarsenVariable( var, vardims, cl ):
  
  if( cl == 1 ): return var

  nd = len(vardims)

  tL = ['time', 't+', 'time+']; timefound = False
  for ts in tL:
    if( ts in vardims ): timefound = True

  if( nd==4 ):
    var = var[:,::cl,::cl,::cl]
  elif(nd==3 and not timefound ):
    var = var[::cl,::cl,::cl]
  elif(nd==3 and timefound     ):
    var = var[:,::cl,::cl]
  elif(nd==2 and not timefound ):
    var = var[::cl,::cl]
  elif(nd==2 and timefound     ):
    var = var[:,::cl]
  elif(nd==1 and not timefound ):
    var = var[::cl]
  else:
    pass
  
  return var 

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def read3DVariableFromDataset(varStr, ds, iTOff=0, iLOff=0, iROff=0, cl=1, meanOn=False):
  
  # iLOff: left offset
  # iROff: right offset
  # cl   : coarsening level
  varStr = partialMatchFromList( varStr , ds.variables.keys() )
  print(' Reading variable {} ... '.format(varStr))
  var, dDict = readVariableFromDataset(varStr, ds, cl=1 )
  print(' ... done.')


  iL = 0 + int(iLOff/cl)
  iR = int(abs(iROff/cl))
  iT = 0 + int(iTOff)
  if(iR == 0):
    # Param list (time, z, y, x )
    if(meanOn):
      vo = var[iL:, iL:, iL:]
    else:
      vo = var[iT:, iL:, iL:, iL:]
  else:
    if(meanOn):
      vo = var[iL:-iR, iL:-iR, iL:-iR]
    else:
      vo = var[iT:, iL:-iR, iL:-iR, iL:-iR]

  var = None

  return vo, np.array(vo.shape)

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def read3dDataFromNetCDF( fname, varNames, cl=1, zeroNans=True ):
  '''
  Establish two boolean variables which indicate whether the created variable is an
  independent or dependent variable in function createNetcdfVariable().
  '''
  parameter = True;  variable  = False

  '''
  Create a NETCDF input dataset (ds), and its associated lists of dependent (varList)
  and independent (dimList) variables.
  '''
  ds, varList, paramList = netcdfDataset(fname)
  vDict = dict()
  
  for vn in varNames:
    vname = partialMatchFromList( vn , varList ) # Obtain the correct name 
    print(' Extracting {} from dataset in {} ... '.format( vname, fname ))
    var, dDict = readVariableFromDataset(vname, ds, cl )
    print(' {}_dims = {}\n Done!'.format(vn, var.shape ))
    
    vDict[vn] = var   # Store the variable under the provided name 
    
    # Rename the keys in dDict to simplify the future postprocessing
    if( 'time' not in vDict.keys() ):
      vDict['time'] = popKeyFromDict('time', dDict)
      
    if( 'x' not in vDict.keys() ):
      vDict['x'] = popKeyFromDict('x', dDict, True) # firstLetterOnly = True
      
    if( 'y' not in vDict.keys() ):
      vDict['y'] = popKeyFromDict('y', dDict, True) # firstLetterOnly = True
    
    if( 'z' not in vDict.keys() ):
      vDict['z'] = popKeyFromDict('z', dDict, True) # firstLetterOnly = True

  dDict = None
    
  for dn in vDict.keys():
    if( zeroNans ):
      idNan = np.isnan(vDict[dn])
      vDict[dn][idNan] = 0.

  return vDict

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


def interpolatePalmVectors(v0, vc_dims, cmpStr, meanOn=False):

  icmp = int()
  iOn = False
  jOn = False
  kOn = False
  kCopy = False
  if(cmpStr == 'i'):
    icmp = 3
    iOn = True
  elif(cmpStr == 'j'):
    icmp = 2
    jOn = True
  elif(cmpStr == 'k'):
    icmp = 1
    kOn = True
  elif(cmpStr == 'kc'):
    icmp = 1
    kCopy = True
  else:
    print('Invalid component string: {}. Exiting ...'.format(cmpStr))
    sys.exit(1)
    
  vMasked = False
  
  if( np.ma.is_masked(v0) ):
    vMasked = True
    vc = np.ma.zeros( vc_dims )
    mc = np.zeros( vc_dims[1:], bool ) # Skip time axis 
  else:
    vc = np.zeros(vc_dims)
  
  if( meanOn ):
    if( vMasked ): vm = np.ma.zeros(vc_dims[1:])
    else:          vm = np.zeros(vc_dims[1:])
  else:
    vm = np.array([])

  # Create index arrays for interpolation.
  jl = np.arange(0, vc_dims[icmp]); jr = jl + 1  # x,y,z: left < right

  nTo,    nzo, nyo, nxo = np.shape(v0)
  nTimes, nz,  ny,  nx  = vc_dims
  
  if( nz == nzo ): k1 = 0
  else:            k1 = 1
  
  for i in range(nTimes):
    if( vMasked ): tmp0 = v0[i, :, :, :].data.view()
    else:          tmp0 = v0[i, :, :, :].view()

    if(iOn):
      tmp1 = (tmp0[:, :, jl] + tmp0[:, :, jr]) * 0.5
      tmp2 = tmp1[k1:, 0:-1, :].view()
    if(jOn):
      tmp1 = (tmp0[:, jl, :] + tmp0[:, jr, :]) * 0.5
      tmp2 = tmp1[k1:, :, 0:-1].view()
    if(kOn):
      tmp1 = (tmp0[jl, :, :] + tmp0[jr, :, :]) * 0.5
      tmp2 = tmp1[:, 0:-1, 0:-1].view()
    if( kCopy ):
      tmp1 = tmp0[jl, :, :].view()
      tmp2 = tmp1[:, 0:-1, 0:-1].view()

    vc[i, :, :, :] = tmp2

    if(meanOn):
      vm += tmp2

  # Clear memory.
  tmp0 = None; tmp1 = None; tmp2 = None

  if( vMasked ):
    tmp0 = v0[nTo//2,:,:,:].mask.view()
    
    if( iOn ):
      tmp1 = (tmp0[:, :, jl] * tmp0[:, :, jr])
      tmp2 = tmp1[k1:,0:-1, :].view()
    if(jOn):
      tmp1 = (tmp0[:, jl, :] * tmp0[:, jr, :])
      tmp2 = tmp1[k1:, :, 0:-1].view()
    if(kOn):
      tmp1 = (tmp0[jl, :, :] * tmp0[jr, :, :])
      tmp2 = tmp1[:, 0:-1, 0:-1].view()
    if( kCopy ):
      tmp1 = tmp0[jl, :, :].view()
      tmp2 = tmp1[:, 0:-1, 0:-1].view()
    
    mc[:,:,:] = tmp2
    
    # Clear memory again
    tmp0 = None; tmp1 = None; tmp2 = None
    
    mc = np.reshape( mc, (1, nz, ny, nx))
    if( np.abs( v0.fill_value ) > 1e3 ): 
      mc += ( np.abs(vc[-1,:,:,:]) > 0.4*np.abs(v0.fill_value) )
    
    vc.mask = np.repeat( mc, nTimes, axis=0 ) # Seems wasteful but appears to be mandatory
    np.ma.set_fill_value( vc , v0.fill_value )
    mc = None  # Clear 


  if(meanOn):
    vm *= float(nTimes)**(-1)
    if( vMasked ):
      vm.mask = vc.mask[0,:,:,:]

  print(' Interpolation along the {}^th direction completed.'.format(cmpStr))

  return vc, vm    # vm is empty if meanOn=False.

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def vectorPrimeComponent(vc, vm):
  vc_dims = np.shape(vc)
  vp = np.zeros_like( vc )

  nTimes = vc_dims[0]
  print(' Computing primes for {} times ... '.format(nTimes))

  for i in range(nTimes):
    vp[i, :, :, :] = vc[i, :, :, :] - vm[:, :, :]

  print(' ... done.')

  return vp

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# NOTE below: variables after *args are interpreted as keyword arguments with default values.

def createNetcdfVariable(dso, v, vName, vLen, vUnits, vType, vTuple, parameter, *args,\
  zlib=False, fill_value=-9999., verbose=True, mask_value=-9999.):
  
  if(parameter):
    dso.createDimension(vName, vLen)
  
    
  var = dso.createVariable(vName, vType, vTuple, zlib=zlib, fill_value=fill_value)
  var.units = vUnits

  if(not np.ma.is_masked(v) and v is not None):
    # Convert variable into masked array. Use [mask_value] to define the mask. 
    # This ensures the fill_values are inserted correctly in netcdf4 function createVariable().
    v = v.view( np.ma.MaskedArray )
    if( (mask_value is not np.nan)  and (mask_value is not None)): 
      v.mask = (v==mask_value)
    elif( mask_value is np.nan ):
      v.mask = np.isnan(v)
    else:
      v.mask = False

  if v is not None:
    var[:] = v
    v = None

  if(parameter): pStr = 'parameter'
  else:          pStr = 'variable'

  if(verbose):
    print(' NetCDF {} {} successfully created. '.format(pStr, vName))

  return var

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


def createCoordinateAxis(dso, Rdims, Rdpx, axis, varname, formatstr, unit, parameter, zlib=False, verbose=True, offset=0.0):
  arr = np.empty(Rdims[axis])
  for i in range(Rdims[axis]):
    # dpx is in [N,E], see getGeoTransform() in gdalTools.py
    arr[i] = np.maximum(0.0, i + offset) * Rdpx[axis]
  axvar = createNetcdfVariable( \
    dso, arr, varname, len(arr), unit, formatstr, (varname,), parameter, zlib, verbose=verbose )
  arr = None
  return axvar

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


def fillTopographyArray(Rtopo, Rdims, Rdpx, datatype):
  topodims = np.array([Rdims[2], Rdims[0], Rdims[1]])
  topo = np.zeros(topodims, dtype=datatype)
  print(' \n Filling 3D array from topography data...')
  print(' Dimensions [z,y,x]: [{}, {}, {}]'.format(*topodims))
  print(' Total number of data points: {}'.format(np.prod(topodims)))
  for x in range(Rdims[1]):
    for y in range(Rdims[0]):
      # Reverse the y-axis because of the top-left origo in raster
      maxind = int(round(Rtopo[-y - 1][x] / Rdpx[2]))+1
      if(maxind>1):
        topo[0:maxind, y, x] = 1
  print(' ...done. \n')
  return topo

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def read3dDictVarFromNetCDF( fname, nameDict, cl=1 ):
  '''
  Establish two boolean variables which indicate whether the created variable is an
  independent or dependent variable in function createNetcdfVariable().
  '''
  parameter = True;  variable  = False

  '''
  Create a NETCDF input dataset (ds), and its associated lists of dependent (varList)
  and independent (dimList) variables.
  '''
  ds, varList, paramList = netcdfDataset(fname)

  '''
  Read cell center coordinates and time.
  Create the output independent variables right away and empty memory.
  '''
  vn = nameDict['varname']
  time, time_dims = read1DVariableFromDataset('time', vn, ds, 0, 0, 1 ) # All values.
  x, x_dims = read1DVariableFromDataset(nameDict['xname'], vn, ds, 0, 0, cl )
  y, y_dims = read1DVariableFromDataset(nameDict['yname'], vn, ds, 0, 0, cl )
  z, z_dims = read1DVariableFromDataset(nameDict['zname'], vn, ds, 0, 0, cl )
  x[np.isnan(x)] = 0.  # Clear away NaNs
  y[np.isnan(y)] = 0.  #
  z[np.isnan(z)] = 0.  #
  '''
  Read in the velocity components.
  PALM netCDF4:
  u(time, zu_3d, y, xu)
  v(time, zu_3d, yv, x)
  w(time, zw_3d, y, x)
  '''
  print(' Extracting {} from dataset ... '.format( vn ))
  v, v_dims = read3DVariableFromDataset(vn, ds, 0, 0, 0, cl) # All values.
  print(' {}_dims = {}\n Done!'.format(vn, v_dims ))

  dataDict = dict()
  dataDict[vn] = v
  dataDict['x'] = x
  dataDict['y'] = y
  dataDict['z'] = z
  dataDict['time'] = time

  return dataDict

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
