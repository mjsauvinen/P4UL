#!/usr/bin/env python3
import argparse
import numpy as np
import sys
from netcdfTools import *
'''Calculates forces from given netCDF file. Currently only pressure drag from
mask output is implemented.

Jukka-Pekka Keskinen, 2021

'''

#==========================================================#
parser = argparse.ArgumentParser(prog='extractBoundaryForces.py')
parser.add_argument('-f', '--filenames',type=str, nargs='+', help=
                    'Input files. All should have the same time axis.')
parser.add_argument('-A', '--cellArea',type=float, default=1.0,
                    help='Cell area for pressure force calculations. '
                    'Default: 1m².')
args = parser.parse_args()

#==========================================================#
# Initial renaming operations and variable declarations

filenames = args.filenames
area = args.cellArea
#==========================================================#


ds, varList, paramList = netcdfDataset(filenames[0],False)

if 'p' in varList:
  print('  Pressure found from variables, calculating drag.')
  print('  Using cell area '+str(area)+' in pressure force calculation \n')
  aikakoko = ds['p'].shape[0]
  ds0 = None
  varList = None
  paramList = None

  for ti in range(aikakoko):
    psumx = 0.0
    psumy = 0.0
    for tiedosto in filenames:
      ds, varList, paramList = netcdfDataset(tiedosto,False)
      itr = ds['p'].shape
      for zi in range(itr[1]):
        for yi in range(itr[2]):
          p = ds['p'][ti,zi,yi,:].data
          pm = ds['p'][ti,zi,yi,:].mask
          if np.size(pm) > 1:
            raki = np.where(pm)[0][0]-1
            if raki > 0:
              psumx = psumx + p[raki]*area
            raki = np.where(pm)[0][-1]+1
            if raki < np.size(p):
              psumx = psumx - p[raki]*area
        for xi in range(itr[3]):
          p = ds['p'][ti,zi,:,xi].data
          pm = ds['p'][ti,zi,:,xi].mask
          if np.size(pm) > 1:
            raki = np.where(pm)[0][0]-1
            if raki > -1:
              psumy = psumy + p[raki]*area 
            raki = np.where(pm)[0][-1]+1
            if raki < np.size(p):
              psumy = psumy - p[raki]*area
    print('  * * * Drag * * *')
    print('  time:  '+np.array2string(ds['time'][ti],precision=0)+' s')
    print('  X:     '+np.array2string(psumx, precision=1)+' N')
    print('  Y:     '+np.array2string(psumy, precision=1)+' N')
    print('  total: '+np.array2string(np.sqrt(psumx**2+psumy**2), precision=1)+' N')
else:
  print('  Pressure not found from variables, nothing to do.')
