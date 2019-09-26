#!/usr/bin/python3
import sys
import pylab as pl
import argparse
from utilities import filesFromList
from plotTools import userLabels, plotCSV
''' 
A script to extract and plot data from a CSV file.

Author: Mikko Auvinen, Department of Applied Mechanics, Aalto University 
'''

#==========================================================#
parser = argparse.ArgumentParser()
parser.add_argument("strKey", help="Search string for collecting files.",nargs='?',\
    default=".csv")
parser.add_argument("--magy", help="Magnitude of all variables.", action="store_true",\
    default=False)
parser.add_argument("--yx", help="Reverse axes: plot(x,y) --> plot(y,x)", action="store_true",\
    default=False)    
parser.add_argument("--labels", help="User specified labels.", action="store_true",\
    default=False)
parser.add_argument("--reuse", help="Reuse once specified variable selections.", action="store_true",\
    default=False)
#==========================================================#
args = parser.parse_args()    
strKey = args.strKey


while 1:

  fileNos, fileList = filesFromList( "*"+strKey+"*" )

  pfig = pl.figure(num=1, figsize=(8.5,8.5));
  for fn in fileNos:
    pfig = plotCSV( pfig, fileList[fn], args.yx, args.magy, args.reuse )
    
  if( args.labels ): pfig = userLabels( pfig )
  pl.grid(True)
  pl.legend(loc=0)

  pl.show()


