#!/home/mjsauvin/python/bin/python
import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
from plotTools import addToPlot
from spectraTools import *
from netcdfTools import *


'''
Kaimal & Finnigan:
The requirements for averaging time T with T >> Tau_{\alpha} can then be expressed
in terms of \sigma^{2}_{\bar{\alpha}}, the variance of the measured time mean \bar{\alpha}
about the expected ensemple mean, and \sigma^{2}_{\alpha}, the ensemble variance of \alpha.
'''
#==========================================================#
sepStr = ' # = # = # = # = # = # = # = # = '
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filename", type=str, help="Name of the input NETCDF file.")
parser.add_argument("-v", "--varname",  type=str, help="Name of the variable in NETCDF file.",\
  default='u')
#parser.add_argument("-T", "--DeltaT", help="Length of period to be analyzed in seconds [s].",\
#  type=float)
parser.add_argument("-nb", "--nbins", help="Number of frequency bins.", type=int, default=76)
parser.add_argument("-n", "--normalize", help="Compute f*S/sigma^2.", action="store_true", \
  default=False)
parser.add_argument("-xn", "--xname",type=str, help="Specify the x coordinate. e.g. xu or x",\
  default='x')
parser.add_argument("-yn", "--yname",type=str, help="Specify the y coordinate. e.g. yv or y",\
  default='y')
parser.add_argument("-zn", "--zname",type=str, help="Specify the z coordinate. e.g. zu_3d or zw_3d",\
  default='zu_3d')
parser.add_argument("-p", "--printOn", help="Print the numpy array data.",\
  action="store_true", default=False) 
parser.add_argument("-pp", "--printOnly", help="Only print the numpy array data. Don't save.",\
  action="store_true", default=False)
parser.add_argument("-c", "--coarse", type=int, help="Coarsening level. Int > 1.",\
  default=1)
args = parser.parse_args()    
#==========================================================# 
# Rename ...
filename  = args.filename
varName   = args.varname
normalize = args.normalize
Nbins     = args.nbins
cl        = abs(args.coarse)
xname     = args.xname
yname     = args.yname
zname     = args.zname
#==========================================================# 

'''
Establish two boolean variables which indicate whether the created variable is an
independent or dependent variable in function createNetcdfVariable().
'''
parameter = True;  variable  = False
''' 
Create a NETCDF input dataset (ds), and its associated lists of dependent (varList)
and independent (dimList) variables. 
'''
ds, varList, paramList = netcdfDataset(filename)
'''
Read cell center coordinates and time.
Create the output independent variables right away and empty memory.
'''
time, time_dims = read1DVariableFromDataset('time', ds, paramList, 0, 0, 1 ) # All values.
x, x_dims = read1DVariableFromDataset( xname,ds, paramList, 0, 0, cl )
y, y_dims = read1DVariableFromDataset( yname,ds, paramList, 0, 0, cl ); print(' y_dims = {} '.format(y_dims))
y[np.isnan(y)] = 0.  # Special treatment.
z, z_dims = read1DVariableFromDataset( zname ,ds, paramList, 0, 0, cl ) # Exclude the first value.

'''
Read in the velocity components.
PALM netCDF4: 
  u(time, zu_3d, y, xu)
  v(time, zu_3d, yv, x)
  w(time, zw_3d, y, x)
'''

v0, v0_dims = read3DVariableFromDataset( varName , ds, varList, 0, 0, 0, cl ) # All values.

infoStr = '''
 Coord. range:
 min(x)={0} ... max(x)={1}, nx = {2}
 min(y)={3} ... max(y)={4}, ny = {5}
 min(z)={6} ... max(z)={7}, nz = {8}
'''.format(np.min(x), np.max(x), len(x), \
  np.min(y), np.max(y), len(y),\
    np.min(z), np.max(z), len(z))

print(infoStr)
ixyz = input(" Enter ids: ix, iy, iz = ")
if( len(ixyz) != 3 ):
  sys.exit(' Error! You must provide 3 values for ix, iy and iz. Exiting ...')

ix = np.minimum( ixyz[0] , len(x)-1 ); ix = np.maximum( ix , 0 )
iy = np.minimum( ixyz[1] , len(y)-1 ); iy = np.maximum( iy , 0 )
iz = np.minimum( ixyz[2] , len(z)-1 ); iz = np.maximum( iz , 0 )

v = v0[:,int(iz),int(iy),int(ix)]

nterms  = np.shape(v) # Number of variables in the file, number of entries in the data. 
print(' Number of terms in the data series, Nv = {}'.format(nterms))

# Determine the sampling frequency.
samplingFreq   = samplingFrequency( time, None )
print(' sampling frequency = {}'.format(samplingFreq))

DeltaT = (time[-1]-time[0])
vw     = applyTapering( v , DeltaT , samplingFreq  )

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
