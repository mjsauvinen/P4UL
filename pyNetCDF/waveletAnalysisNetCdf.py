#!/usr/bin/env python3
import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
from plotTools import addToPlot, addImagePlotDict
from analysisTools import sensibleIds, groundOffset, discreteWaveletAnalysis
from waveletTools import wtDataset, Morlet
from netcdfTools import read3dDataFromNetCDF
from utilities import filesFromList
from txtTools import openIOFile
from scipy import signal, ndimage
''' 
Description: A script to perform wavelet analysis on data stored in a NETCDF file.

Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
        
        Simone Boi
        University of Helsinki
'''
#==========================================================#
sepStr = ' # = # = # = # = # = # = # = # = '
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filename", type=str,\
  help="Name of the input NETCDF file.")
parser.add_argument("-fo", "--fileout", type=str, default="out.nc", \
  help="Name of the output NETCDF file. Default=out.nc")
parser.add_argument("-v", "--varname",  type=str, default='u',\
  help="Name of the variable in NETCDF file. Default='u' ")
parser.add_argument("-m", "--mode", type=str, default='spectro', \
  choices=['spectrohist', 'scalohist', 'spectro', 'scalo'],\
  help="Mode: histo-, spectro-, or scalo-gram.")
parser.add_argument("-P", "--power", action="store_true", default=False,\
  help="Plot power histo-, spectro- or scalo-gram.")
parser.add_argument("-C", "--complex", action="store_true", default=False,\
  help="Complex spectro- or scalo-gram.")
parser.add_argument("-k", "--kIndices",type=int, nargs=2,\
  help="Starting and final index (k_start, k_final) of the considered data. Required.")
parser.add_argument("-s", "--save", type=str, default=None, \
  help="Identifier (name) for the saved figures. Default=None")
parser.add_argument("-of", "--outputToFile", type=str, default=None, \
  help="Name of the file to output analysis results. Default=None")
parser.add_argument("--gridOn", action="store_true", default=False,\
  help="Grid on the figure.")
parser.add_argument("--limsOn", action="store_true", default=False,\
  help="User defined limits.")
args = parser.parse_args()    
#==========================================================# 
# Rename ...
filename  = args.filename
fileout   = args.fileout
varname   = args.varname
kIds      = args.kIndices
mode      = args.mode
power     = args.power
complexOn = args.complex 
ofile     = args.outputToFile
saveFig   = args.save
gridOn   = args.gridOn
limsOn   = args.limsOn
#==========================================================# 
'''
Establish two boolean variables which indicate whether the created variable is an
independent or dependent variable in function createNetcdfVariable().
'''
parameter = True;  variable  = False


# Read data from NETCDF file
cl = 1
ncDict = read3dDataFromNetCDF( filename , [varname] , cl )
v = ncDict.pop(varname) 

# Spatial coords and time
x = ncDict.pop('x')
y = ncDict.pop('y')
z = ncDict.pop('z')
time = ncDict.pop('time')

# Plot coord. information. This aids the user in the beginning.
infoStr = '''
  Coord. range:
  min(x)={0} ... max(x)={1}, nx = {2}
  min(y)={3} ... max(y)={4}, ny = {5}
  min(z)={6} ... max(z)={7}, nz = {8}
'''.format(\
    np.min(x), np.max(x), len(x),\
    np.min(y), np.max(y), len(y),\
    np.min(z), np.max(z), len(z) )
#print(infoStr)

# Now check whether the given indices make sense 
# Here we set i = 0 and j = 0.
ijk1 = sensibleIds( np.array([0,0,kIds[0]]), x, y, z )
ijk2 = sensibleIds( np.array([0,0,kIds[1]]), x, y, z )
print(' Check (1): i, j, k = {}'.format(ijk1))
print(' Check (2): i, j, k = {}'.format(ijk2))

Nt = 4. 
dtime = np.mean( time[1:] - time[:-1] )
ltime = (time[-1] - time[0])/Nt
f1    = 1./ltime
f2    = 1./(2.*dtime)
nsteps= int( (f2-f1)/f1 )
freq  = np.linspace(f1,f2,nsteps,endpoint=False)


print(' f1 = {}, f2 = {}, nsteps = {} '.format(f1,f2,nsteps))
print(' freq = {} '.format( freq ))


sig  = v[:,kIds[0], ijk1[1], ijk1[0]] - np.mean( v[:,kIds[0], ijk1[1], ijk1[0]] )
fig1 = plt.figure()
plotStr = [" Signal "," time "," v(time) "]
fig1 = addToPlot(fig1, time, \
  ndimage.gaussian_filter(sig, sigma=60.),'Filtered', plotStr )
fig1 = addToPlot(fig1, time,sig,'Signal', plotStr )



WD1 = wtDataset(sig,  time, f=freq )

if( complexOn ):  tt = "complex"
else:             tt = "real"


if('hist' in mode):
  Qstr = '''    Enter the index of desired {} values
    given the range values {}...{} with indecies {}...{}.'''
  if( 'spectro' in mode ):
    histmode = 'frequency'
    print(Qstr.format(histmode,freq[0],freq[-1], 0, len(freq)-1))
  else:
    histmode = 'scale'
    print(Qstr.format(histmode,WD1.scales[0],WD1.scales[-1], 0, len(WD1.scales)-1))
    
  idx = input('    Use comma (,) to separate index values: ')
  if( not isinstance( idx, tuple ) ):
    if( isinstance( idx, int ) ):
      idx = (idx,)
    else:
      sys.exit(' Error: wrong entry for index values ')
      
  for i in idx:
    pk, bins, patches = WD1.MorletHistogram(i, mode=histmode, nbins=32, plotHistogram=True )
    pk=map(lambda x: x*(bins[1]-bins[0]), pk)
    
    

else:
  wDict  = dict()
  exDict = WD1.getBounds()
  
  if( mode == 'spectro' ):
    
    if( power ): 
      wDict['R'] = WD1.PowerMorletSpectrogram()
      wDict['title'] = ' Power Morlet Spectrogram '
      bounds = exDict['freq']
    else:
      wDict['R'] = WD1.SigMorletSpectrogram(ttype=tt)
      wDict['title'] = ' {} Morlet Spectrogram '.format(tt.upper())
      bounds = exDict['freq']
      
    wDict['ylabel'] = 'frequency (Hz)'
    
      
  elif( mode == 'scalo' ):
    if( power ): 
      wDict['R'] = WD1.PowerMorletScalogram()
      wDict['title'] = ' Power Morlet Scalogram '
      bounds = exDict['scales']
    else:
      wDict['R'] = WD1.SigMorletScalogram(ttype=tt)
      wDict['title'] = ' {} Morlet Scalogram '.format(tt.upper())
      bounds = exDict['scales']
      
    wDict['ylabel'] = 'scales (s)'

  wDict['extent'] = bounds
  wDict['xlabel'] = 'time (s)'
  wDict['gridOn'] = gridOn; wDict['limsOn'] = limsOn
  wDict['cmap'] = 'bwr'
  
  fig2 = plt.figure()
  fig2 = addImagePlotDict(fig2, wDict )
  
  '''
  # Set up the axes with gridspec
  figx = plt.figure(figsize=(6, 6))
  grid = plt.GridSpec(4, 4, hspace=0.2, wspace=0.2)
  main_ax = figx.add_subplot(grid[:-1, 1:])
  y_hist = figx.add_subplot(grid[:-1, 0], xticklabels=[], sharey=main_ax)
  x_hist = figx.add_subplot(grid[-1, 1:], yticklabels=[], sharex=main_ax)

  # scatter points on the main axes
  main_ax.imshow(np.real(wDict['R']), aspect='auto')

  # histogram on the attached axes
  x_hist.hist(x, 40, histtype='stepfilled',
            orientation='vertical', color='gray')
  x_hist.invert_yaxis()

  y_hist.hist(y, 40, histtype='stepfilled',
            orientation='horizontal', color='gray')
  y_hist.invert_xaxis()
  '''


plt.show()

'''

# = = = = = = = = = = = = =
# Mean and variance
vmean = np.mean(v, axis=(0))

# Extract fluctuating part and normalize by variance
# Reuse the v variable
v -= vmean

# Assemble the quadrant analysis dict 
wlDict = dict()
wlDict['wavelet'] = wavelet
wlDict['nlevel']  = nlevel
# = = = = = = = = = = = = = = =
kList = np.arange(kIds[0],kIds[1]+1)
for kt in kList:
  vt = v[:,kt, ijk1[1], ijk1[0]]
  
  values, labels = discreteWaveletAnalysis( vt, wlDict )


  fig = plt.figure()
  fig.subplots_adjust(hspace=0.2, bottom=.03, left=.07, right=.97, top=.92)
  plt.subplot(2, 1, 1)
  plt.title(" signal")
  plt.plot(time, vt, 'b')
  plt.xlim(time[0], time[-1])

  ax = plt.subplot(2, 1, 2)
  plt.title("Wavelet packet coefficients at level {}".format(nlevel))
  plt.imshow(values, interpolation='nearest', aspect="auto",\
    origin="lower", extent=[0, 1, 0, len(values)])
  plt.yticks(np.arange(0.5, len(labels) + 0.5), labels)
  #plt.setp(ax.get_xticklabels(), visible=False)

  #plt.figure(2)
  #plt.specgram(data, NFFT=64, noverlap=32, cmap=cmap)
  #plt.imshow(values, origin='upper', extent=[-1,1,-1,1],
  # interpolation='nearest')
  
  if( saveFig ):
    plt.savefig( saveFig+'_{}m.jpg'.format(int(z[kt])), format='jpg', dpi=300)
    
  fig = None

'''


'''
# = = output file = = = = =
# Create a NETCDF output dataset (dso) for writing out the data.
dso = netcdfOutputDataset( fileout )
xv = createNetcdfVariable( dso, xa  , 'x'   , len(xa)   , 'm', 'f4', ('x',)   , parameter )
yv = createNetcdfVariable( dso, ya  , 'y'   , len(ya)   , 'm', 'f4', ('y',)   , parameter )
zv = createNetcdfVariable( dso, ka  , 'z'   , len(ka)   , 'm', 'f4', ('z',)   , parameter )
Qv = createNetcdfVariable( dso, Qa, 'Q', dims[0], 'm-2', 'f4',('z','y','x',) , variable )

# - - - - Done , finalize the output - - - - - - - - - -
netcdfWriteAndClose( dso )
'''
