#!/usr/bin/env python3
import sys
import pylab as pl
import numpy as np
import argparse
from plotTools import extractFromRAW
from utilities import filesFromList, basicAnalysis


#==========================================================#

def reduceByRadius(x,y,a,Rlim=1.0,xo=0., yo=0.):
  b = []   # Empty list
  r = np.sqrt( (x-xo)**2 + (y-yo)**2 )
  for i in range(len(r)):
    if( r[i] <= Rlim ):
      b.append(a[i])
  #print ' len(a) = {}, len(b) = {}'.format(len(a), len(b))
  return np.array(b)

#==========================================================#
#==========================================================# 

'''
nargs:
  - N (an integer). N arguments from the command line will be gathered together into a list.
  -'?'. One argument will be consumed from the command line if possible, and produced as a single item. If no command-line argument is present, the value from default will be produced.
  -'*'. All command-line arguments present are gathered into a list.
  -'+'. Just like '*', all command-line args present are gathered into a list. Additionally, an error message will be generated if there was not at least one command-line argument present.
'''

sepStr = ' # = # = # = # = # = # = # = # = '
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filename", help="Name of the target .raw file",\
  default=None)
parser.add_argument("-v", "--var", help="Variable Name in RAW-file", nargs='+',\
  default=["U_z"])
parser.add_argument("-r","--ref", help="Reference value", nargs='?', type=float,\
  default=0.)

args = parser.parse_args()    

print type(args.var), args.var
varList = ["x","y","z"]
varList.extend(args.var)#; print varList
v_ref = args.ref#;         print v_ref


if( not args.filename ):
  fileNos, fileList = filesFromList( "*.raw" )
  filename = fileList[0]
else:
  fileList = [args.filename]

for fx in fileList:
  data = extractFromRAW( fx , varList )

  x=data[0]; y=data[1]; z=data[2]
  v=data[3] # Only one var allow at this point.

  RadiusLimit = 0.616   
  vd = reduceByRadius(x,y,v,Rlim=0.6)

  #strOut = 
  print " "; print sepStr
  print 'Analyzing {0} from file {1}.'.format( args.var[0] , fx )
  print 
  v_a = basicAnalysis(vd, args.var[0], v_ref, 1 ) 
  print " "; print sepStr
