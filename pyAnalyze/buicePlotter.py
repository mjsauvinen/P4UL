#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
import math
import glob
import numpy as np
import pylab as pl
from plotTools import filesFromList, plotDY

# = = = FUNCTION DEFINITIONS = = = = =

def plotBuice( fig, fileStr, logOn ):
  x = np.loadtxt(fileStr)
  y = x[:,1]
  ax  = fig.add_axes( [0.115, 0.075 , 0.85 , 0.81] ) #[left, up, width, height], fig.add_subplot(111)

  # Print each column separately
  for i in xrange((x.shape[1]-3)):
    if( logOn ):
      lines=ax.semilogy(np.abs(x[:,i+3]), y[:] , linewidth=1.1 , label=fileStr)
    else:
      lines=ax.plot(x[:,i+3], y[:], linewidth=1.1, label=fileStr )

  ax.set_xlabel(" X/H + 10*U/U_in ",fontsize=22)
  ax.set_ylabel(" Y ",fontsize=22)
  return fig

# - - - - - - -

def plotExprData( fig, fileStr ):
  x = np.loadtxt(fileStr)
  y = x[:,1]

  for xp in xrange(0,41):
    if( str(xp) in fileStr ):
      xL = float(xp)
      print 'x = ', xL
  
  if( '-' in fileStr ):
    xL *= -1.

  ax = fig.add_axes( [0.115, 0.075 , 0.85 , 0.81] ) #[left, up, width, height], fig.add_subplot(111)
  for i in xrange((x.shape[1]-3)):
    lines=ax.plot((xL+10.*x[:,i+3]), y[:],'o', linewidth=1.5, label=fileStr )
  ax.set_xlabel(" X/H + 10*U/U_in ",fontsize=22)
  ax.set_xticks([-6,6,14,20,30,40])
  ax.set_ylabel(" Y ",fontsize=22)
  return fig

# = = = MAIN PROGRAM = = = = =

try:
  strKey = sys.argv[1]
except:
  strKey = "."

pl.rc('xtick', labelsize=16)
pl.rc('ytick', labelsize=16)
pl.rc('ytick.major', size=10)
pl.rc('ytick.minor', size=6)


#os.chdir(path)

while 1:

  pfig = pl.figure(num=1, figsize=(8.5,8.5))

  # == Computational Data == #
  fileNos, fileList = filesFromList( strKey )

  for fn in fileNos:
    pfig = plotBuice( pfig, fileList[fn], 0 )

  # == Experimental Data == #
  fileNos, fileList = filesFromList( "Expr_x" )
  for fn in fileNos:
    pfig = plotExprData( pfig, fileList[fn] )

  pl.grid(True)
  pl.title(' Buice Eaton Diffuser \n Velocity Profile Comparison ', fontsize=20)
  pl.legend()

  pl.show()



