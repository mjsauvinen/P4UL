#!/usr/bin/env python3
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import rfft, irfft, fftfreq
from utilities import filesFromList, writeLog
''' 
Description: Fourier filter (time,scalar) time series given via txt file.


Author: Mikko Auvinen
        mikko.auvinen@fmi.fi
        Finnish Meteorological Institute
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='fourierFilterData.py')
parser.add_argument("fileKey", default=None,\
  help="Search string for collecting files.")
parser.add_argument("-o", "--outstr",type=str, default="ff_",\
  help="Prefix for the output file. Default=fflt_.")
parser.add_argument("-lf", "--cutfreq", type=float, default=0.01,\
  help="Frequency cutoff. FFT coefs will be zeroed for frequecies above this value.")
parser.add_argument("-fx", "--factorX", type=float, default=1.0,\
  help="Multiplication factor for x-values: fx*x")
parser.add_argument("-fy", "--factorY", type=float, default=1.0,\
  help="Multiplication factor for y-values: fy*y")
parser.add_argument("-sr", "--skiprows", type=int, default=1,\
  help="Skip rows when reading the file. Default=1")
parser.add_argument("-mw", "--marginwidth", type=float, default=0.0,\
  help="Relative width [0-1] of optional zero-valued margin added to the front of the signal. Default=0")
parser.add_argument("-p", "--printOn", help="Print the resulting time series data.",\
  action="store_true", default=False)
args = parser.parse_args()
writeLog( parser, args )
#==========================================================#
# Initial renaming operations and variable declarations

fileKey    = args.fileKey 
outstr     = args.outstr
cutfreq    = args.cutfreq
fx         = args.factorX
fy         = args.factorY
sr         = args.skiprows
mw         = args.marginwidth
printOn    = args.printOn

# Obtain a list of files to include.
fileNos, fileList = filesFromList( fileKey+'*' )
for fn in fileNos:
  
  fileout = outstr+fileList[fn]

  try:    t, s = np.loadtxt( fileList[fn], skiprows=sr, usecols=(0,1), delimiter=',', unpack=True)  
  except: t, s = np.loadtxt( fileList[fn], skiprows=sr, usecols=(0,1), unpack=True )

  
  t *= fx; s *= fy
  N = len(s)
  m = int( mw*N )
  
  s2 = np.zeros( N*2 + m )
  s2[:m] = s[0]
  s2[m:(N+m)] = s[:]
  s2[(N+m):] = s[::-1]
  
  vfreq = fftfreq( len(s2) , d=np.mean( t[1:]-t[:-1] ) )
  Fs    = rfft( s2 ); s2 = None
  Fs[(np.abs(vfreq) > cutfreq)] = 0
  sf = irfft( Fs )
  sf = sf[m:(N+m)]
  sf[:4] = s[:4] # First four identical to original data
  idz = (sf < 0.)
  sf[idz] = 0.; idz = None
  #Fs = None
  
  print(' Writing out file: {} '.format( fileout ) )
  np.savetxt(fileout, np.c_[t/fx,sf], fmt='%3.6e')
  
  if(printOn):
    plt.figure(1)
    plt.plot(t/fx,s ,'b', t/fx, sf,'r')
    #plt.figure(2)
    #plt.semilogy(vfreq, Fs)
    plt.show()

print(' Done! ')
