#!/usr/bin/env python
import sys
import pylab as pl
import numpy as np
import argparse
from utilities import filesFromList, basicAnalysis
from plotTools import extractFromCSV
from wave import waveInformation

#==========================================================#

def waveAnalysis( U , xc ):
  t = 1.; Ulim = 0.2; fileout = file('U-bundles.dat','w')
  waveInformation(t, U , xc, Ulim, fileout )
  
#==========================================================#
sepStr = ' # = # = # = # = # = # = # = # = '
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filename", help="Name of the target .csv file", default=None)
parser.add_argument("-v", "--var", help="Variable Name in CSV-file", nargs='+', default=["U:1"])
parser.add_argument("-r", "--ref", help="Reference value", nargs='?', type=float, default=0.)
parser.add_argument("-w", "--waveMode", help="Discrete Wave Mode (overrides -v/--var option)",\
  action="store_true", default=False)

args = parser.parse_args()    
#==========================================================# 
#print type(args.var), args.var

if( not args.filename ):
    fileNos, fileList = filesFromList( "*.csv" )
    try: args.filename = fileList[0]
    except: sys.exit('Could not obtain a valid CSV file. Exiting ...')
    
if( args.waveMode ):
  varList = ["U:1","U:2", "arc_length"]
  dat = extractFromCSV( args.filename , varList )
  Umag = np.sqrt( dat[0]**2 + dat[1]**2 )
  xL   = dat[2]
  waveAnalysis( Umag, xL ) 
  
else:  
  varList = args.var
  v_ref = args.ref; print v_ref
  v = extractFromCSV( args.filename , varList )[0]

  print " "; print sepStr
  print 'Analyzing {0} from file {1}.'.format( args.var[0] , args.filename )
  print 
  v_a = basicAnalysis(v, args.var[0], v_ref, 1 ) 
  print " "; print sepStr

