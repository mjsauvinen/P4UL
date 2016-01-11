#!/usr/bin/env python
import sys
import glob
import numpy as np
import pylab as pl
import argparse
from utilities import filesFromList
from plotTools import addToPlot, userLabels



# =========== === FUNCTIONS === ======================== #
def isolateNumericalValues( s , stripChars ):
    for i in xrange(len(s)):
        for sc in stripChars:
            s[i] = s[i].strip(sc)
    
    return np.array(map( float, s )) # Turns it into column array

# ===================================================== #

def extractColumns( data , icols  ):
    ncols = len(icols); nrows = np.shape(data)[0]
    v = np.zeros( (ncols, nrows) , float) # cols and rows reversed here!
    
    for j in xrange(ncols):
        v[j] = isolateNumericalValues( data[:,icols[j]] , ['(',')'])
    
    return v

# ============ === ARGS === =========================== #

parser = argparse.ArgumentParser()
parser.add_argument("-t","--torgue", help="Plot Torgues", action="store_true",\
    default=False)
parser.add_argument("--labels", help="User specified labels.", action="store_true",\
    default=False)
args = parser.parse_args()

# ============== === MAIN === ======================= #

fileNos, fileList = filesFromList('./postProcessing/forces/*/forces*'  )

Ftot = [[],[],[]]; Ttot =[[],[],[]]
time = []

for fn in fileNos:
    raw = np.loadtxt( fileList[fn], dtype=str )
    ncols = len(raw[0,:])-1  # How many columns of Force data

    t  = np.array( extractColumns( raw , [0,] )[0] )
    Fp = np.array( extractColumns( raw , [1,2,3] ) )
    Fv = np.array( extractColumns( raw , [4,5,6] ) )
    Fo = np.array( extractColumns( raw , [7,8,9] ) )
    Tp = np.array( extractColumns( raw , [10,11,12] ) )
    Tv = np.array( extractColumns( raw , [13,14,15] ) )
    To = np.array( extractColumns( raw , [16,17,18] ) )

    time.extend(t)
    for n in xrange(3):
        Ftot[n].extend( Fp[n]+Fv[n]+Fo[n] )
        Ttot[n].extend( Tp[n]+Tv[n]+To[n] )

    #print np.shape(t), np.shape(time), np.shape(Fp[0]),np.shape(Ftot[0])     

fig = pl.figure(1, figsize=(8.5,8.5))
if( not args.torgue ):
    plotTxt = ["Forces", "Iter/Time (s)", "F (N)"]
    fig = addToPlot(fig, time, Ftot[0],'Fx', plotTxt, False)
    fig = addToPlot(fig, time, Ftot[1],'Fy', plotTxt, False)
    fig = addToPlot(fig, time, Ftot[2],'Fz', plotTxt, False)
else:
    plotTxt = ["Torques", "Iter/Time (s)", "T (Nm)"]
    fig = addToPlot(fig, time, Ttot[0],'Tx', plotTxt, False)
    fig = addToPlot(fig, time, Ttot[1],'Ty', plotTxt, False)
    fig = addToPlot(fig, time, Ttot[2],'Tz', plotTxt, False)

if(args.labels):
    fig = userLabels( fig )
pl.legend(loc=0)
pl.grid(True)


pl.show()
