#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt
import argparse
from utilities import filesFromList
from plotTools import userLabels, plotXX
''' 
Description: A script to plot files with data in [x, y1, y2, ..., yn] format.
Run the script in the directory where you have the data files:
$ uniplot.py <search string> [options]


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''
#=======MAIN PROGRAM========================================#
parser = argparse.ArgumentParser()
parser.add_argument("strKey", help="Search string for collecting files.",\
  nargs='?', default=None)
parser.add_argument("--log", help="Logarithmic y-axis.", action="store_true",\
  default=False)
parser.add_argument("--labels", help="User specified labels.", action="store_true",\
  default=False)
parser.add_argument("-fx", "--factorX", help="Multiplication factor for x-values: fx*x",\
  type=float, default=1.0)
parser.add_argument("-fy", "--factorY", help="Multiplication factor for y-values: fy*y",\
  type=float, default=1.0)
args = parser.parse_args()
#==========================================================#

if( not args.strKey ): 
  args.strKey = raw_input(" Enter search string: ")
  if( not args.strKey ): sys.exit(1)

plt.rc('xtick', labelsize=16); #plt.rc('ytick.major', size=10)
plt.rc('ytick', labelsize=20); #plt.rc('ytick.minor', size=6)

while 1:

  fileNos, fileList = filesFromList( "*"+args.strKey+"*" )

  pfig = plt.figure(num=1, figsize=(9.,9.));
  for fn in fileNos:
    pfig = plotXX( pfig, fileList[fn], args.log, args.factorX, args.factorY )

  if(args.labels):
    pfig = userLabels( pfig )
  plt.grid(True)
  plt.legend(loc=0)

  plt.show()
