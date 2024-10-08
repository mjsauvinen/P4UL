#!/usr/bin/env python3
import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
from spectraTools import timeSeriesFromFiles, spectraAnalysis
from utilities import filesFromList

'''
Kaimal & Finnigan:
The requirements for averaging time T with T >> Tau_{\alpha} can then be expressed
in terms of \sigma^{2}_{\bar{\alpha}}, the variance of the measured time mean \bar{\alpha}
about the expected ensemple mean, and \sigma^{2}_{\alpha}, the ensemble variance of \alpha.
'''
#==========================================================#
parser = argparse.ArgumentParser()
parser.add_argument("fileKey", default=None,\
  help="Search string for collecting files.")
parser.add_argument("-iv", "--icolv", nargs='?', type=int, default=1,\
  help="Column index for analyzed signal (first = 0).")
parser.add_argument("-vn","--varname", type=str, default=None,\
  help="Variable name to appear in the figures.")
parser.add_argument("-m", "--mode", type=str, default='S', choices=['S', 'E', 'P'],\
  help="Mode: 'S': power spectral density, 'E': energy spectrum, 'P': power spectrum.")

# -- group for mutually exclusive entries -- #
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-it", "--icolt", type=int,\
  help="Opt 1: Column index for time in file (first = 0).")
group.add_argument("-fs", "--samplingFreq", type=float,\
  help="Opt 2: Sampling frequency in [Hz].")
# ---    end group              ------------ #

parser.add_argument("-nb", "--nbins", type=int, default=76,\
  help="Number of frequency bins. Default=76")
parser.add_argument("-n", "--normalize", action="store_true", default=False,\
  help="Compute f*S/$\sigma^2$.")
parser.add_argument("-p", "--printOn", action="store_true", default=False,\
  help="Print the numpy array data.")
parser.add_argument("-pp", "--printOnly", action="store_true", default=False,\
  help="Only print the numpy array data. Don't save.")
args = parser.parse_args()
#==========================================================# 
# Rename ...
fileKey   = args.fileKey
jcolV     = args.icolv
jcolT     = args.icolt
normalize = args.normalize
varname   = args.varname
samplingFreqUser  = args.samplingFreq
mode      = args.mode
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

fig = None
fig = spectraAnalysis(fig, v, time, varname, Nbins, mode, normalize)

plt.legend(loc=0)
plt.show()
