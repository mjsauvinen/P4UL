#!/usr/bin/env python3
from utilities import sortTimes, extractMatchingTerms
from paraTools import extractParaviewLineData
from wave import waveInformation
import os, sys
import numpy as np
import argparse
    
# = # = # = # = # = Arguments and their parsing # = # = # = # = # = # 

parser = argparse.ArgumentParser()

parser.add_argument("-f","--filename", type=str, help="Name of the output file.", default="waveInfo")
parser.add_argument("-U","--UName", type=str, help="Name of liquid phase U. (default=U)", \
    default="U")
parser.add_argument("-a","--alphaName", type=str, help="Name of liquid phase alpha. (default=alpha.liquid)",\
    default="alpha.liquid")
parser.add_argument("-p1","--point1", type=float, action='append', nargs=3, \
    help="Coords. of point 1:  x1 y1 z1")
parser.add_argument("-p2","--point2", type=float, action='append', nargs=3, \
    help="Coords. of point 2:  x2 y2 z2")
parser.add_argument("-s","--screen", help="Output values ONLY to screen.", action="store_true",\
    default=False)
parser.add_argument("-t","--times", type=float, action='append', nargs='+', \
    help="Specific times: t1 t2 t3 ...", default=None)

args = parser.parse_args()


# = # = # = # = # = # = # = # = # = # = # = #= # = #= # = #
# = Extract and employ the data from argument parser.

UName = args.UName
alphaName = args.alphaName

p1 = args.point1[0]; p2 = args.point2[0]
print " p1 = {}, p2 = {}".format(p1,p2)

fo = None
if( not args.screen ): fo = file("{}.dat".format(args.filename), 'w')


# = # = # = # = # = # = # = # MAIN = # = # = # = #= # = #= # = #

baseDir = os.getcwd()
timeDirs = sortTimes() 
timeFloats = map( float , timeDirs ) # Back to numbers (floats)

if( args.times ):
    #print " times = {}".format(args.times)
    times = extractMatchingTerms( args.times[0] , timeFloats , 1 )
    if( len(times) == 0 ): 
        sys.exit('None of the times found! Exiting ...')
else:
    times = timeFloats


varList = [ UName, alphaName ]
csvVars = ["U.liq","Points:0"]
calcSettings = { "Function":UName+"_X*"+alphaName, "ResultArrayName": csvVars[0] } 

line = [p1, p2]



Uw_min = 1.e-3   # This threshold value may need finetuning.
for t in times: 
    print 'time = {}'.format(t)
    [Ul, x] = extractParaviewLineData( baseDir, varList, line, calcSettings, csvVars, t )
    #print ' size( Ul )= {}'.format(np.size(Ul))

    waveInformation(t, Ul, x, Uw_min, fo )


if( isinstance(fo, file) ): fo.close()
