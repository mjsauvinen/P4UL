#!/usr/bin/env python
from utilities import filesFromList
from utilities import writeLog
from plotTools import addContourf, addToPlot
from footprintTools import *
from mapTools import readNumpyZTile, filterAndScale, farFieldIds
import sys
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
''' 
Description:


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''
# = # = # = # Function definitions # = # = # = # = # = # = # 
# = # = # = # = # = # = # = # = # = # = # = # = # = # = #
# = # = # = # End Function definitions # = # = # = # = # = #

#========================================================== #
parser = argparse.ArgumentParser(prog='footprintDiff.py')
parser.add_argument("-fr", "--fileRef", type=str,help="Ref Footprint file. (npz format)")
parser.add_argument("-fd", "--filesDiff", type=str,\
  help="Search string for diff footprint files. (npz format)", default='.npz')
help_x ="Percentage of farthest x-koords excluded from the diff. (Default=None)." 
parser.add_argument("-x","--excl", type=float, default=None, help=help_x)
parser.add_argument("-pp", "--printOnly", help="Only print the contour. Don't save.",\
  action="store_true", default=False)
parser.add_argument("--save", help="Save the figure right away.", action="store_true",\
  default=False)
args = parser.parse_args() 
writeLog( parser, args, args.printOnly )
#========================================================== #

# Rename ... that's all.
fileRef   = args.fileRef
filesDiff = args.filesDiff
excl      = args.excl
printOnly = args.printOnly
saveOn    = args.save            #  save the fig

# = = = = = = = = = = = = = = = = = = = = = = = = = = = =   #
# xO := origin coords. # xt := target coords. # ut := target speed
try:
  Fref, X, Y, Z, C = readNumpyZFootprint( fileRef ) # IdsOn=True
except:
  sys.exit(' Could not read the footprint file: {}'.format(fileRef))

# Gather footprint data files: 
fileNos, fileList = filesFromList( "*"+filesDiff.strip(".npz")+"*.npz" )

for fn in fileNos:
  try:
    Fi, X, Y, Z, Ci = readNumpyZFootprint( fileList[fn] )
  except:
    sys.exit(' Could not read the footprint file: {}'.format(fileList[fn]))

  # Resolution:
  dPx = np.array([ (X[0,1]-X[0,0]) , (Y[1,0]-Y[0,0]) ])

  # Normalize s.t. global integral becomes one.
  print(' Normalizing the footprints such that SUM(Fp) = 1 ...')
  C1 = 1./np.sum( Fref  * np.prod(dPx));  Fref *= C1
  C2 = 1./np.sum( Fi    * np.prod(dPx));  Fi   *= C2
  print('... done! C1_ref = {} and C2 = {}'.format(C1, C2))

  
  print(' DIFF: {} vs. {}'.format(fileRef, fileList[fn]))
  dF2   = np.zeros( Fref.shape )
  idx   = np.ones( Fref.shape, bool )
  if( excl ): idx -= farFieldIds( X, excl ) # Subtracting True values gives False
  
  idNz       = (Fref != 0.) * idx
  FrefMean   = np.mean( Fref[idNz] )
  dF2[idx]   = ( (Fref[idx] - Fi[idx])*np.prod(dPx) )**2 #/ FrefMean
  dnorm = np.sqrt( np.sum(dF2) )   #/np.prod( np.shape(Fref) )
  print(' norm of diff = {}\n'.format(dnorm))
  print('#--------------------#')
  
  

  