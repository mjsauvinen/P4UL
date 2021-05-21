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

plt.rcParams["legend.labelspacing"] = 1.2
plt.rcParams["legend.framealpha"]   = 1.
plt.rcParams["legend.edgecolor"]    = 'k'
plt.rcParams["legend.fontsize"] = 18
plt.rcParams["legend.handleheight"]   = 1.0  # default: 0.7
plt.rcParams["legend.handlelength"]   = 2.2  # default: 2.0

#=======MAIN PROGRAM========================================#
parser = argparse.ArgumentParser()
parser.add_argument("-fs", "--fileSearchStr", nargs='+', required=True,\
  help="File search strings for collecting files for each subplot. Required.")
parser.add_argument("-hz","--horizontal", action="store_true", default=False,\
  help="Horizontal stacking of subplots.")
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
parser.add_argument("-nl","--nolegend", action="store_true", default=False,\
  help="Do not draw a legend.")
parser.add_argument("-all", "--allfiles", help="Select all files automatically.",\
  action="store_true", default=False)
parser.add_argument("-yl", "--ylims", type=float, nargs=2, default=[None,None],\
  help="Bounds (limits) for the y axes")
parser.add_argument("-xl", "--xlims", type=float, nargs=2, default=[None,None],\
  help="Bounds (limits) for the x axes")
parser.add_argument("-lm", "--linemode", type=int, choices=[1,2], default=1,\
  help="Mode for displaying the color and type of lines. See the source code. Default=1")
parser.add_argument("-lw", "--linewidth", type=float, default=2.6,\
  help="Line width. Default=2.6")
parser.add_argument("-s", "--save", type=str, default=None, \
  help="Name of the saved figure. Default=None")
args = parser.parse_args()
#==========================================================#
# Rename ...
fsearch    = args.fileSearchStr
horizontal = args.horizontal
vertical   = ~horizontal
Cx         = args.factorX
Cy         = args.factorY
logOn      = args.log
labelsOn   = args.labels
allfiles   = args.allfiles
saveFig    = args.save
Nplots     = len( fsearch )

fsize = (11., Nplots*9.)

if( not horizontal ):
  pfig, axs = plt.subplots(Nplots, 1, sharex=True, sharey=True, figsize=fsize)
else:
  pfig, axs = plt.subplots(1, Nplots, sharex=True, sharey=True, figsize=fsize[::-1])

# Assemble plot dict
pD = dict()  
pD['logOn']    = logOn; 
pD['Cx']       = Cx
pD['Cy']       = Cy
pD['revAxes']  = args.revAxes
pD['lm']       = args.linemode
pD['ylims']    = args.ylims
pD['xlims']    = args.xlims
pD['lw']       = args.linewidth
pD['reset']    = True

# loop over each subplot
for i in range(Nplots):

  fileNos, fileList = filesFromList( fsearch[i]+"*", allfiles )

  pD['reset'] = True
  for fn in fileNos:
    pD['filename'] = fileList[fn]
    pfig = plotXX( pfig, pD, axs[i] )
    pD['reset'] = False

  axs[i].grid(True)
  axs[i].legend(loc=1) # upper right ... for now

if( labelsOn ):
  print(' userLabels ')
  pfig = userLabels( pfig )

plt.tight_layout()

if( saveFig is not None ):
  pfig.savefig( saveFig, format='jpg', dpi=300)
else:
  plt.show()

pfig = None
