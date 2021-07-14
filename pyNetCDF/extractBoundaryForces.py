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
parser.add_argument('-o', '--outputFile',type=str, default=None, help=
                    'Save output forces to the given netCDF file.')
args = parser.parse_args()

#==========================================================#
# Initial renaming operations and variable declarations

filenames = args.filenames
outputfile = args.outputFile
area = args.cellArea
#==========================================================#


ds, varList, paramList = netcdfDataset(filenames[0],False)

if 'p' in varList:
  print('  Pressure found from variables, calculating drag.')
  print('  Using cell area '+str(area)+' in pressure force calculation \n')
  aikakoko = ds['p'].shape[0]
  aika = ds['time'][:].data
  ds.close()
  ds = None
  varList = None
  paramList = None

  Dx=-9999.0*np.ones(aikakoko)
  Dy=-9999.0*np.ones(aikakoko)

  print('  * * * Drag * * *')
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
      ds.close()      
    Dx[ti] = psumx
    Dy[ti] = psumy
    print('  time: '+np.array2string(aika[ti],precision=0)+' s')
    print('     X:     '+np.array2string(psumx, precision=0)+' N')
    print('     Y:     '+np.array2string(psumy, precision=0)+' N')
    print('     total: '+np.array2string(np.sqrt(psumx**2+psumy**2), precision=0)+' N')
  if not outputfile==None:
    dso = netcdfOutputDataset(outputfile)
    print('\n Outputting drag data to file '+outputfile)
    createNetcdfVariable(dso, aika, 'time', aikakoko, 'seconds', 'f8', ('time',), True)
    createNetcdfVariable(dso, Dx, 'Drag x', aikakoko, 'Newtons', 'f8', ('time',), False,
                         fill_value=-9999.0)
    createNetcdfVariable(dso, Dy, 'Drag y', aikakoko, 'Newtons', 'f8', ('time',), False,
                         fill_value=-9999.0)
    dso.close()

else:
  print('  Pressure not found from variables, nothing to do.')
