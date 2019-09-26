#!/usr/bin/env python3
import sys
import glob
import numpy as np
import argparse
import subprocess as sb
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filename", help="Name of the target .vtk files",\
  default='U_zPlane.vtk')
parser.add_argument("-s", "--str2num", help="String to be replaced by a number.",\
  default='zPlane')
parser.add_argument("-ts", "--timescale", help="Factor for scaling time in file name",\
  type=float, default=0.)

args = parser.parse_args()

factor = args.timescale
origName = args.filename
str2num  = args.str2num
timeDirs = glob.glob('[0-9]*')  # Extracts all time directories as str list.
timeValues = map( float, timeDirs )
tmax = np.max(timeValues)

# Determine timescale factor automatically if none is provided. 
if( factor == 0 ):
    if (tmax < 10.):
        factor = 1000
    elif(tmax < 100.):
        factor = 100
    elif(tmax < 1000.):
        factor = 10
    else:
        factor = 1

for i in xrange(len(timeDirs)):
    timeStr = '{:g}'.format(int(timeValues[i]*factor))
    newName = origName.replace(str2num, timeStr )
    cmd = 'mv {}/{} ./{}'.format(timeDirs[i], origName, newName )
    print cmd
    retcode = sb.call( cmd , shell=True)

print ' Completed !'

