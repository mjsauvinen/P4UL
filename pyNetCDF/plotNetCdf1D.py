#!/usr/bin/env python
from netcdfTools import *
from utilities import selectFromList, removeEntriesFromList
from plotTools import addToPlot, userLabels
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
''' 
Description:


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='plotNetCdf1D.py')
parser.add_argument("-f", "--filename",type=str, help="Name of the input NETCDF file.")
parser.add_argument("--log", help="Logarithmic y-axis.", action="store_true",\
  default=False)
parser.add_argument("--labels", help="User specified labels.", action="store_true",\
  default=False)
args = parser.parse_args() 
#==========================================================#
# Initial renaming operations and variable declarations

filename = args.filename
labelsOn = args.labels
logOn    = args.log

'''
Establish two boolean variables which indicate whether the created variable is an
independent or dependent variable in function createNetcdfVariable().
'''
parameter = True;  variable = False

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
''' 
Create a NETCDF input dataset (ds), and its associated lists of dependent (varList)
and independent (paramList) variables. 
'''
ds, varList, paramList = netcdfDataset(filename, False )
varList = removeEntriesFromList(varList, paramList)


while 1:
  print(' Select variables for the y-axis.')
  iyList = selectFromList( varList )
  
  # Construct the plot.
  fig = plt.figure(num=1, figsize=(12.,9.5))
  for iy in iyList:
    # Read data
    y, xdict = readVariableFromDataset( varList[iy], ds, varList)
    
    if('time' not in xdict.keys()): 
      print(' time not in dict = {} '.format(xdict.keys()))
      print(' adding time=[None] ...')
      xdict['time'] = np.array([None])
      y = y.reshape((1,len(y)))  # to ensure len(np.shape(y) == 2
    
    # If time is the only parameter, use it on x-axis
    if(len(np.shape(y)) == 1):
      xStr = 'time'
      labelStr = ' {}({}) '.format(varList[iy],'time')
      plotTxt   = [labelStr, xStr, varList[iy]]
      fig = addToPlot(fig, xdict[xStr], y, labelStr, plotTxt, logOn)
    
    elif(len(np.shape(y)) == 2):
      time  = xdict['time']
      xList = xdict.keys(); xList.remove('time') # remove time ...
      xStr = xList[0]                            # and take what is left.
      # Plot each time instance of the profile
      for j in xrange(len(time)):
        if( time[j] == None ): 
          labelStr  = ' {}({})'.format(varList[iy],xStr)
        else:          
          labelStr  = ' {0}(time={1:.0f} s, {2})'.format(varList[iy],time[j],xStr)
        
        idy = (y[j,:] < 1.e+9)
        if( np.count_nonzero(idy) > (len(idy)/3)):
          plotTxt   = [labelStr, xStr, varList[iy]]
          fig = addToPlot(fig, xdict[xStr][idy], y[j,idy], labelStr, plotTxt, logOn)
        else:
          pass
    
    else:
      sys.exit(' Plotting of profiles with more than one spatial dimension not supported.')

  
  if(labelsOn):
    fig = userLabels( fig )
  plt.legend(loc=0)
  plt.show()

print(' Bye! ')