#!/usr/bin/env python
from utilities import filesFromList
from utilities import writeLog
from plotTools import addContourf, addToPlot
from footprintTools import *
from mapTools import readNumpyZTile, filterAndScale
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
parser = argparse.ArgumentParser(prog='footprintMaskOutput.py')
parser.add_argument("-f", "--filename", type=str,help="Footprint file. (npz format)")
parser.add_argument("-fm", "--filemask", type=str, help="Mask file. (npz format)")
parser.add_argument("-pp", "--printOnly", help="Only print the contour. Don't save.",\
  action="store_true", default=False)
parser.add_argument("--save", help="Save the figure right away.", action="store_true",\
  default=False)
args = parser.parse_args() 
writeLog( parser, args, args.printOnly )
#========================================================== #

# Rename ... that's all.
filename  = args.filename
filemask  = args.filemask
printOnly = args.printOnly
saveOn    = args.save            #  save the fig

# = = = = = = = = = = = = = = = = = = = = = = = = = = = =   #
# xO := origin coords. # xt := target coords. # ut := target speed
try:
  Fp, X, Y, Z, C, IDict = readNumpyZFootprint( filename, True ) # IdsOn=True
except:
  sys.exit(' Could not read the footprint file: {}'.format(filename))

try:
  Rm, Rmdims, RmOrig, dPx = readNumpyZTile( filemask )
except:
  sys.exit(' Could not read the mask file: {}'.format(filemask))


# To unify the treatment, let's add 100% to the IDict.
if(not IDict): IDict = {}  # In case IDict comes in as <None>
IDict[100] = np.ones( np.shape(Fp) , bool ) # All true


Nm = int(np.max(Rm))+1  # Number of different mask ID's (assuming first is zero)
dA = np.prod(dPx)

# Read the user specified source strengths
Qd  = np.ones(Nm, float)    # Default 
try:
  Q = input("Enter {} source strengths (Q) separated by commas: ".format(Nm))
except:
  print(" All source strengths are set to unity. Q[:] = 1.")
  Q = Qd
if( len(Q) != Nm ):
  sys.exit(" Error!  len(Q) = {}. It should be {}. Exiting ...".format(len(Qe),Nm))


for key in IDict.keys():
  # CSV header
  idx = IDict[key]
  print('{}%:\n \"Mask ID\",\"[%]\",\"SUM( Fp*M*dA )\",\"SUM( FP*dA )\" ,\" Q \"'.format(key))
  #Fptot = np.sum(Fp[idx]*dA)
  Fptot  = 0.
  FpM    = np.zeros( Nm )
  for im in xrange(Nm):
    M = (Rm == im).astype(int)
    #print(' sum(Fp*idx) = {}, min(M) = {}'.format(np.min(Fp*idx), np.min(M)))
    FpM[im] = Q[im] * np.sum(Fp*idx*M*dA)
    Fptot   += FpM[im]
  for im in xrange(Nm):
    pStr = '{}, {}, {}, {}, {}'.format(im,FpM[im]/Fptot*100.,FpM[im],Fptot,Q[im])
    print(pStr)
  print('----------')

mpl.rcParams['font.size'] = 18.0
plt.figure(num=1, figsize=(9.,6.));
lbl = np.array(['Buildings','Impervious','Grass',\
  'Low Vegitation','High Vegitation', 'Water', '','Road'])
dat = FpM/Fptot*100.
ix = (dat > 0.)  # Valid points. In order to isolate meaningless entries.
#colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral']
cs=np.array(['lightskyblue', 'g', 'r', 'c', 'm', 'y', 'k', 'w'])
expl = np.zeros(Nm)
expl[2:] = 0.1



plt.pie(dat[ix],explode=expl[ix],labels=lbl[ix],colors=cs[ix],\
  autopct='%1.1f%%',shadow=True, startangle=90)
#plt.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%', shadow=True, startangle=90)
# Set aspect ratio to be equal so that pie is drawn as a circle.
plt.axis('equal')

if(saveOn):
  plt.savefig( filename.strip('.npz')+'.jpg' )
plt.show()

#print(' Footprint: {} '.format(np.shape(Fp)))
#print(' Mask     : {} '.format(np.shape(Rm)))
#print(' FpM    : {} '.format(FpM)) 



