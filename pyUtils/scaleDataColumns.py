#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse
from utilities import filesFromList, inputIfNone
''' 
Description: A script to scale data in files with [x, y1, y2, ..., yn] format.
Scaling is done: vs = (y - v0)/v*


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''
#=======MAIN PROGRAM========================================#
parser = argparse.ArgumentParser()
parser.add_argument("strKey", nargs='?', default=None,\
  help="Search string for collecting files.")
parser.add_argument('-x', '--prefix', type=str, default='nond_',\
  help="Prefix for file names. Default='nond_'  .. as in non dimensional")
parser.add_argument('-ic', '--icols', type=int, nargs='+',\
  help='Columns to scale. Ex: 0,1,2,4')
parser.add_argument("-v0", "--vRef", type=float, nargs='+', default=[0.0],\
  help="Reference values 'v0' in v+ = (v - v0)/v* for each column in --icols. Default = [0.]")
parser.add_argument("-vs", "--vStar", type=float, nargs='+', default=[1.0],\
  help="Characteristic value 'v*' in v+ = (v - v0)/v* for each column in --icols. Default = [1.]")
args = parser.parse_args()
#==========================================================#
# Rename ...
strKey  = args.strKey
v0      = np.array( args.vRef  )  # Convert to numpy array
vs      = np.array( args.vStar )
icols   = args.icols
prefix  = args.prefix
#==========================================================#

vx0 = np.zeros( np.shape(icols), float )
vxS = np.ones(  np.shape(icols), float )

if( (len(v0) == 1) and (len(icols) != 1) ):
  vtmp = np.zeros( np.shape(icols), float )
  vtmp[:] = v0[0]
  v0 = vtmp.copy(); vtmp = None
elif( len(v0) != len(icols) ):
  sys.exit(' Error: incompatible array lengths for columns and v0 values. Exiting.')

if( (len(vs) == 1) and (len(icols) != 1) ):
  vtmp = np.ones( np.shape(icols), float )
  vtmp[:] = vs[0]
  vs = vtmp.copy(); vtmp = None
elif( len(vs) != len(icols) ):
  sys.exit(' Error: incompatible array lengths for columns and vs values. Exiting.')
  
#==========================================================#

strKey = inputIfNone( strKey , " Enter search string: " )

fileNos, fileList = filesFromList( strKey+"*")

for fn in fileNos:
  try:    dat = np.loadtxt(fileList[fn])
  except: dat = np.loadtxt(fileList[fn], delimiter=',')
  
  j = 0
  for i in icols:
    if i not in range(dat.shape[1]):
      sys.exit(' Error: {} not in range({}). Exiting ...'.format(i,dat.shape[1]) )
    
    dat[:,i] = ( dat[:,i] - v0[j] )/vs[j]
    j += 1
    
  fileout = prefix+fileList[fn]
  print(' Writing out file: {} '.format( fileout ) )
  
  headerStr = 'Scaling done with v+ = (v - v0)/v* with v0 = {} and v* = {} for columns = {}'.format(v0, vs, icols)
  np.savetxt(fileout, dat[:,:], fmt='%3.6e', header=headerStr)
  dat = None

print(' All done! ')