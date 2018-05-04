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
parser.add_argument("--tavg", help="Time-averaged values (if applicable).", action="store_true",\
  default=False)
parser.add_argument("-ts", "--timeskip",type=int, default=0, \
  help="Number of skipped time instances (if applicable).")
parser.add_argument("--labels", help="User specified labels.", action="store_true",\
  default=False)
parser.add_argument("-wa", "--writeAscii", help="Write X, Y[-1,:] (or time averaged Y) data to an ascii file.",\
  action="store_true", default=False)
parser.add_argument("-waa", "--writeAllAscii", help="Write X, Y data to an ascii file.",\
  action="store_true", default=False)
args = parser.parse_args() 
#==========================================================#
# Initial renaming operations and variable declarations

filename = args.filename
tskip    = args.timeskip
labelsOn = args.labels
logOn    = args.log
timeAverageOn = args.tavg
writeAscii = args.writeAscii or args.writeAllAscii
writeAll   = args.writeAllAscii

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
    y, xdict = readVariableFromDataset( varList[iy], ds, 1)
    
    if('time' not in xdict.keys()): 
      print(' time not in dict = {} '.format(xdict.keys()))
      print(' adding time=[None] ...')
      xdict['time'] = np.array([None])
      y = y.reshape((1,len(y)))  # to ensure len(np.shape(y)) == 2
    
    # If time is the only parameter, use it on x-axis
    if(len(np.shape(y)) == 1):
      xStr = 'time'
      labelStr = ' {}({}) '.format(varList[iy],'time')
      plotTxt   = [labelStr, xStr, varList[iy]]
      fig = addToPlot(fig, xdict[xStr], y, labelStr, plotTxt, logOn)
      if( writeAscii ):
        print(' (1) Writing data to ascii file: {}.dat'.format(varList[iy])) 
        np.savetxt(varList[iy]+'.dat', np.c_[xdict[xStr],y] )   # x,y,z equal sized 1D arrays
    
    elif(len(np.shape(y)) == 2):
      time  = xdict['time']
      xList = xdict.keys(); xList.remove('time') # remove time ...
      xStr = xList[0]                            # and take what is left ... z-something, typically
      tskip = min(tskip, len(time)-1)
      
      
      # Plot profiles for chosen time instances
      if( not timeAverageOn ):
        for j in xrange(tskip, len(time)):
          if( time[j] == None ): 
            labelStr  = ' {}({})'.format(varList[iy],xStr)
          else:          
            labelStr  = ' {0}(time={1:.0f} s, {2})'.format(varList[iy],time[j],xStr)
        
          ieff = (y[j,:] < 1.e+9); ieff[-1] = False  # Set last false to skip it
          if( np.count_nonzero(ieff) > (len(ieff)/3)):
            plotTxt   = [labelStr, varList[iy], xStr]
            fig = addToPlot(fig, y[j,ieff], xdict[xStr][ieff], labelStr, plotTxt, logOn)
          else:
            pass
      
      else: # timeAverageOn
        tskip = max( 1, tskip )  # The first should almost always be skipped.
        labelStr  = '  <{}>_t({})'.format(varList[iy],xStr)
        plotTxt   = [labelStr, varList[iy], xStr]
        fig = addToPlot(fig, np.nanmean(y[tskip:,:-1], axis=0), xdict[xStr][:-1], labelStr, plotTxt, logOn)
        
        
        
      if( writeAscii ):
        print(' (2) Writing data to ascii file: {}.dat'.format(varList[iy]))
        print(' x.shape = {} vs y.shape = {}'.format(np.shape(xdict[xStr]), np.shape(y)))
        if(writeAll):
          hStr = ' {} at times = {}'.format(varList[iy], time )
          np.savetxt(varList[iy]+'.dat', np.c_[xdict[xStr], np.transpose(y[tskip:,:-1])], header=hStr )
        else:
          hStr = ' {} at times = {}'.format(varList[iy], time[-1] )
          if( timeAverageOn):
            np.savetxt(varList[iy]+'_tavg.dat', np.c_[xdict[xStr][:-1], np.transpose(np.nanmean(y[tskip:,:-1], axis=0))], header=hStr)
          else:   
            np.savetxt(varList[iy]+'.dat', np.c_[xdict[xStr][:-1], np.transpose(y[-1,:-1])], header=hStr )
    
    else:
      sys.exit(' Plotting of profiles with more than one spatial dimension not supported.')

  
  if(labelsOn):
    fig = userLabels( fig )
  plt.legend(loc=0)
  plt.show()

print(' Bye! ')
