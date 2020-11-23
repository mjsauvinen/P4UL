#!/usr/bin/env python3

import argparse
import numpy as np
import sys
from utilities import writeLog
from mapTools import readNumpyZTile, saveTileAsNumpyZ
#===================================================================================================#
'''merge3dNpz

A script for combining two 3D objects in npz files. Can be used in e.g. mesh creation.

The merged array in the npz file is S. Perhaps more flexibility on the
matter later on.

Author: Jukka-Pekka Keskinen, FMI, 2020

'''
#==================================================================================================#
parser = argparse.ArgumentParser(prog='merge3dNpz.py',\
  description="Merge two 3D objects within npz files.")
parser.add_argument("-f1", "--filename1", type=str, \
  help="Input file containing the first numpy array which is merged into -f2.")
parser.add_argument("-f2", "--filename2", type=str, \
  help="Name of the file containing the second numpy array onto which -f1 is merged.")
parser.add_argument("-mx", "--maxVal", type=float, default=None,\
  help="Force to maintain given maximum value.")
parser.add_argument("-fo", "--fileout", type=str, default='output.npz',\
  help="Name of the output file.")
parser.add_argument("-ml", "--mergeloc", nargs=3 , type=int, default=[0,0,0],\
  help="Merge location. Indices (jN, iE, kZ) in -f2 where the object from -f1 will be placed.")
parser.add_argument("-a", "--add", action="store_true", default=False,\
  help="Add object in the first file on top of the second file.")
args = parser.parse_args()
writeLog( parser, args )
#==================================================================================================#

f1      = args.filename1
f2      = args.filename2
fileout = args.fileout
maxVal  = args.maxVal
ml      = args.mergeloc
add     = args.add 

# = = = = = = = = = = = = = = = = = #

# Load files
F1=readNumpyZTile(f1)
S1=F1['S']
F1 = None


F2=readNumpyZTile(f2)
S2=F2['S']
dataType = S2.dtype

# Convert the type of S1 into the type of S2.
S1 = S1.astype(dataType)


# Calculate merge indices
me = ml + np.asarray(S1.shape)
me = np.minimum(me, np.asarray(S2.shape))
# Problem here! If me is clipped above, the slicing below fails as shapes will not agree.

if any(np.less_equal(me-ml,0)):
  sys.exit('Unable to merge objects. The first object is out of bounds.')

idn2 = np.isnan(S2)
Nnan = np.count_nonzero(idn2)

if( Nnan > 0 ):
  idn1 = np.isnan(S1)  
  ida  = ~idn2[ml[0]:me[0],ml[1]:me[1],ml[2]:me[2]] * ~idn1  # Where to add
  idr  =  idn2[ml[0]:me[0],ml[1]:me[1],ml[2]:me[2]] * ~idn1  # Where to insert



# Combine S1 with S2
if( add ):  
  if( Nnan > 0 ):
    S2[ml[0]:me[0],ml[1]:me[1],ml[2]:me[2]][ida]+= S1[ida]
    S2[ml[0]:me[0],ml[1]:me[1],ml[2]:me[2]][idr] = S1[idr]
  else:
    S2[ml[0]:me[0],ml[1]:me[1],ml[2]:me[2]]+= S1
else:
  S2[ml[0]:me[0],ml[1]:me[1],ml[2]:me[2]]=S1

S1 = None # clear memory

if( maxVal is not None ):
  S2 = np.minimum( S2, maxVal )

F2['S']=S2

#Save npz
saveTileAsNumpyZ(fileout,F2)
