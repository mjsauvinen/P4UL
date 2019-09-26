#!/usr/bin/env python3
from utilities import filesFromList, writeLog
from footprintTools import writeNumpyZFootprintRaw
import sys
import argparse
import numpy as np

''' 
Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''
# = # = # = # Function definitions # = # = # = # = # = # = # 

# = # = # = # End Function definitions  # = # = # = # = # = #

#========================================================== #
parser = argparse.ArgumentParser(prog='footprintRaw2Npz.py')
parser.add_argument("fileKey",type=str, default='TARGET', \
  help="Search string for raw footprint input files. Default=TARGET ")
parser.add_argument("-fo", "--fileout",type=str, default='FP', \
  help="Name of the .npz footprint output file. Default=FP")
args = parser.parse_args() 
writeLog( parser, args )
#========================================================== #

# Rename ... that's all.
fileKey = args.fileKey
fileout = args.fileout

# = = = = = = = = = = = = = = = = = = = = = = = = = = = =   #

fileNos, fileList = filesFromList( fileKey+'*' )

a = None # Open a storage list which will contain N np.arrays.

print(' Read the raw footprint data files ... \n ... and this may take awhile ... ')
for fn in fileNos:
  if( a == None ):
    a = np.loadtxt( fileList[fn] )
  else:
    a = np.concatenate( (a , np.loadtxt(fileList[fn])) )
  print ' ... array dimensions = {} '.format(a.shape)
print(' ... done !\n')

# Palm output file format:
# origin_x(0), origin_y(1), origin_z(2), x(3), y(4), z(5), speed_x(6), speed_y(7), speed_z(8) 

# Write .npz output file
writeNumpyZFootprintRaw( fileout , a )

a = None  # Clear memory.
