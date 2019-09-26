#!/usr/bin/python3
# -*- coding: utf-8 -*-
''' 
A script to plot files with data in [x, y, z, v1, v2, ..., vn] format.

Author: Mikko Auvinen, Department of Applied Mechanics, Aalto University

Example usage: (1) Place in your $HOME/bin directory (with plotTools.py) and make it an executable.
  (2) Run the script in the directory where you have the data files with an optional file search string:
   cuniplot.py [optional str, default="."] 
'''

import sys
import matplotlib.pyplot as plt
import argparse
from utilities import filesFromList
from plotTools import plotDY, userLabels


# = = = MAIN PROGRAM = = = = =

parser = argparse.ArgumentParser()
parser.add_argument("strKey", help="Search string for collecting files.",\
    default=None)
parser.add_argument("-d","--dim", nargs=1, help="No. of spatial dimensions.", type=int,\
    action="store", default=3)
parser.add_argument("--labels", help="User specified labels.", action="store_true",\
    default=False)
parser.add_argument("--yx", help="Transpose x and y axes.", action="store_true",\
    default=False)

try: 
  args = parser.parse_args()
except: 
  args = parser.parse_args([0])   # Send strKey=0
  args.strKey = raw_input(" Enter search string: ")
  if( not args.strKey ): sys.exit(1)
  
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=16)


while 1:

  fileNos, fileList = filesFromList( args.strKey )

  pfig = plt.figure(num=1, figsize=(8.,8.));
  for fn in fileNos:
    pfig = plotDY( pfig, fileList[fn], args.dim, args.yx )
  
  if( args.labels):
    pfig = userLabels( pfig )
  plt.grid(True)
  plt.legend(loc=0)

  plt.show()



