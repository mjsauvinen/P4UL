#!/usr/bin/env python3

import argparse
import numpy as np
import sys
from utilities import writeLog
from mapTools import readNumpyZTile, saveTileAsNumpyZ
#=============================================================================================================#
'''merge3DObjects

A script for combining two 3D objects in npz files. Can be used in e.g. mesh creation.

The merged array in the npz file is S. Perhaps more flexibility on the
matter later on.

Author: Jukka-Pekka Keskinen, FMI, 2020

'''
#==================================================================================================#
parser = argparse.ArgumentParser(prog='merge3DObjects.py',description="Merge two 3D objects within npz files.")
parser.add_argument("-f1", "--filename1", type=str, \
  help="Name of the file containing the first object. This is placed on top of the object in the second file.")
parser.add_argument("-f2", "--filename2", type=str, \
  help="Name of the file containing the second object. This acts as the background for the object in the first file.")
parser.add_argument("-mx", "--maxVal", type=float, default=None,\
  help="Force to maintain given maximum value.")
parser.add_argument("-fo", "--fileout", type=str, default='output.npz',\
  help="Name of the output file.")
parser.add_argument("-mloc", "--mergeloc", nargs=3 , type=int, default=[0,0,0], \
  help="Merge location. Indices in the second file where the object from first file will be placed.")
parser.add_argument("-a", "--add", action="store_true", default=False,\
  help="Add object in the first file on top of the second file.")
args = parser.parse_args()
writeLog( parser, args )
#==================================================================================================#

maxVal = args.maxVal


# Load files
F1=readNumpyZTile(args.filename1)
S1=F1['S']
F1 = None

F2=readNumpyZTile(args.filename2)
S2=F2['S']

# Calculate merge indices
ms=args.mergeloc
me=ms+np.asarray(S1.shape)
me=np.minimum(me,np.asarray(S2.shape))

if any(np.less_equal(me-ms,0)):
  sys.exit('Unable to merge objects. The first object is out of bounds.')

# Combine S1 with S2
if args.add:
  S2[ms[0]:me[0],ms[1]:me[1],ms[2]:me[2]]+=S1
else:
  S2[ms[0]:me[0],ms[1]:me[1],ms[2]:me[2]]=S1

S1 = None # clear memory

if( maxVal is not None ):
  S2 = np.minimum( S2, maxVal )

F2['S']=S2

#Save npz
saveTileAsNumpyZ(args.fileout,F2)
