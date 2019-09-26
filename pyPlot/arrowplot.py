#!/usr/bin/env python3
''' 
A script to make an arrow plot from files with data in [x, y, dx, dy ] format.

Author: Mikko Auvinen, Department of Applied Mechanics, Aalto University

Example usage: 
  (1) Place in your $HOME/bin directory (with plotTools.py) and make it an executable.
  (2) Run the script in the directory where you have the data files:
    >> arrowplot.py <search string>   
'''
import sys
import matplotlib.pyplot as plt
import numpy as np
from plotTools import filesFromList
from plotTools import userLabels, arrow2DPlot

plt.rc('xtick', labelsize=24)
plt.rc('ytick', labelsize=24)

try:
  strKey = sys.argv[1]
except:
  strKey = "."
  
while 1:
  fileNos, fileList = filesFromList( strKey )
  afig = plt.figure(num=1, figsize=(12.,11.));

  try:
    scale = float( input(' Enter scaling factor. Scale='))
  except:
    scale = 1.

  x  = np.loadtxt( fileList[fileNos[0]] )
  ax = afig.add_axes( [0.075, 0.075 , 0.85 , 0.85] )
  lines=ax.plot( x[:,0], x[:,1], '.k' )

  fillOn = True
  for i in range(len(fileNos)):
    fn = fileNos[i]
    afig = arrow2DPlot( afig, fileList[fn], scale, i , fillOn )
    if( fillOn ): 
      fillOn = not fillOn
  
  pl.show()
    
