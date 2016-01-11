#!/usr/bin/env python
import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
from plotTools import addToPlot
from spectraTools import *
from utilities import filesFromList


'''
Kaimal & Finnigan:
The requirements for averaging time T with T >> Tau_{\alpha} can then be expressed
in terms of \sigma^{2}_{\bar{\alpha}}, the variance of the measured time mean \bar{\alpha}
about the expected ensemple mean, and \sigma^{2}_{\alpha}, the ensemble variance of \alpha.
'''
#==========================================================#
parser = argparse.ArgumentParser()
parser.add_argument("fileKey", help="Search string for collecting files.", default=None)
parser.add_argument("-T", "--DeltaT", help="Length of period to be analyzed in seconds [s].",\
  type=float)
parser.add_argument("-iv", "--icolv", help="Column number for analyzed signal. First is 0.",\
  nargs='?', type=int, default=1)
parser.add_argument("-nb", "--nbins", help="Number of frequency bins.", type=int, default=76)
parser.add_argument("-n", "--normalize", help="Compute f*S/sigma^2.", action="store_true", \
  default=False) 
parser.add_argument("-p", "--printOn", help="Print the numpy array data.",\
  action="store_true", default=False) 
parser.add_argument("-pp", "--printOnly", help="Only print the numpy array data. Don't save.",\
  action="store_true", default=False)

# -- group for mutually exclusive entries -- #
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-it", "--icolt", type=int,\
  help="Column number in file for time (first column is 0).")
group.add_argument("-fs", "--samplingFreq", help="Sampling frequency in [Hz].", type=float)
# ---    end group              ------------ #
args = parser.parse_args()    
#==========================================================# 
# Rename ...
fileKey   = args.fileKey
DeltaT    = args.DeltaT
jcolV     = args.icolv
jcolT     = args.icolt
normalize = args.normalize
samplingFreqUser  = args.samplingFreq
Nbins     = args.nbins
#==========================================================# 

# Read the data from raw file:
jcols = []
if( jcolT != None ): jcols.append(jcolT)
jcols.append(jcolV)

# Obtain a list of files to include.
fileNos, fileList = filesFromList( '*'+fileKey+'*' )

# Read and gather the time and variable arrays.
time, v = timeSeriesFromFiles( fileNos, fileList, jcols )

nterms  = np.shape(v) # Number of variables in the file, number of entries in the data. 
print(' Number of terms in the data series, Nv = {}'.format(nterms))

# Determine the sampling frequency.
samplingFreq   = samplingFrequency( time, samplingFreqUser )
vw             = applyTapering( v , DeltaT , samplingFreq  )

# Evaluate, Power (P), power spectral energy (E), and power spectral density (S).
P, E, S, freqs = evalSpectra( vw, samplingFreq, normalize )

Sbin, fbin = frequencyBins( S , freqs, Nbins )
Ebin, fbin = frequencyBins( E , freqs, Nbins )
Pbin, fbin = frequencyBins( P , freqs, Nbins )

Refbin = 1.E-1 * fbin[Nbins/2:]**(-5./3.)


if( normalize ):
  labelStr = 'Normalized Power Spectral Density '
  plotStr  = [ labelStr ," Frequency [Hz] "," f*S/$sigma^2$ "] 
else:
  labelStr = ' Power Spectral Density: $\Phi(f)$ '
  plotStr  = [ labelStr ," Frequency [Hz] "," $\Phi(f)$ "] 

fig1 = plt.figure(num=1, figsize=(12.,10.))
fig1 = addToPlot(fig1, fbin, Sbin, labelStr, plotStr, logOn=True)
fig1 = addToPlot(fig1, fbin[Nbins/2:], np.nanmean(Sbin)*Refbin , ' Model -5/3 curve', plotStr, logOn=True)
plt.legend(loc=0)


fig2 = plt.figure(num=2, figsize=(12.,10.))
labelStr = ' Energy Spectrum: E(f) '
plotStr  = [" Energy Spectrum "," Frequency [Hz] "," E(f) "] 
fig2 = addToPlot(fig2, fbin, Ebin, labelStr, plotStr, logOn=True)
fig2 = addToPlot(fig2, fbin[Nbins/2:], np.nanmean(Ebin)*Refbin , ' Model -5/3 curve', plotStr, logOn=True)
plt.legend(loc=0)


fig3 = plt.figure(num=3, figsize=(12.,10.))
labelStr = ' Power Spectrum: P(f) '
plotStr  = [" Power Spectrum "," Frequency [Hz] "," P(f) "] 
fig3 = addToPlot(fig3, fbin, Pbin, labelStr, plotStr, logOn=True)
fig3 = addToPlot(fig3, fbin[Nbins/2:], np.nanmean(Pbin)*Refbin , ' Model -5/3 curve', plotStr, logOn=True)
plt.legend(loc=0)

plt.show()
