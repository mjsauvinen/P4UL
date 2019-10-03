#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt
import argparse
from utilities import filesFromList, inputIfNone
from plotTools import userLabels, plotXX
''' 
Description: A script to plot multiple files with data in [x, y1, y2, ..., yn] format.
Run the script in the directory where you have the data files:
$ uniplot.py <search string> [options]


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''
#=======MAIN PROGRAM========================================#
parser = argparse.ArgumentParser()
parser.add_argument("strKey", nargs='?', default=None,\
  help="Search string for collecting files.")
parser.add_argument("--log", action="store_true", default=False,\
  help="Logarithmic y-axis.")
parser.add_argument("--labels", action="store_true", default=False,\
  help="User specified labels.")
parser.add_argument("-rx","--revAxes", action="store_true", default=False,\
  help="Reverse axes: X <--> Y.")
parser.add_argument("-fx", "--factorX", type=float, default=1.0,\
  help="Multiplication factor for x-values: fx*x")
parser.add_argument("-fy", "--factorY", type=float, default=1.0,\
  help="Multiplication factor for y-values: fy*y")
parser.add_argument("-lm", "--linemode", type=int, default=1,\
  help="Mode for displaying the color and type of lines. See the source code. Default=1")
parser.add_argument("-s", "--save", type=str, default=None, \
  help="Name of the saved figure. Default=None")
args = parser.parse_args()
#==========================================================#
# Rename ...
strKey  = args.strKey
factorX = args.factorX
factorY = args.factorY
linemode= args.linemode
revAxes = args.revAxes
logOn   = args.log
labelsOn= args.labels
saveFig = args.save

strKey = inputIfNone( strKey , " Enter search string: " )

styleStr = 'seaborn-white' # 'ggplot'  # 'seaborn-paper'
#plt.style.use(styleStr)



while 1:

  fileNos, fileList = filesFromList( strKey+"*" )

  pfig = plt.figure(num=1, figsize=(12.,9.5));
  for fn in fileNos:
    pfig = plotXX( pfig, fileList[fn], logOn, factorX, factorY, revAxes, linemode )

  if( labelsOn ):
    print(' userLabels ')
    pfig = userLabels( pfig )
  plt.grid(True)
  plt.legend(loc=0)
  
  if( saveFig ):
    pfig.savefig( saveFig, format='jpg', dpi=300)

  plt.show()
  pfig = None
