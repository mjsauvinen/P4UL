#!/usr/bin/env python
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as scs
from plotTools import addToPlot

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def readColsFromFile( filename , cols ):
  
  csvOn = 'csv' in filename 
  
  try:    
    if( csvOn ): dat = np.loadtxt( filename , usecols=cols, delimiter = ',', skiprows=1 )
    else       : dat = np.loadtxt( filename , usecols=cols )
  except: 
    sys.exit(' Error! Cannot read cols={} from {}. Exiting ...'.format(cols, filename))
  
  if( len(cols) == 1 ):
    t = None; v = dat
  else:
    t = dat[:,0]; v = dat[:,1]
  
  print(' v.shape = {}'.format(v.shape))
  
  return t, v

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def timeSeriesFromFiles( fileNos, fileList, cols ):
  tx = None; vx = None
  print(' Reading data from files ... ')
  
  for fn in fileNos:
    if( tx == None and vx == None ):
      tx, vx = readColsFromFile( fileList[fn] , cols )
    else:
      ttmp, vtmp  = readColsFromFile( fileList[fn] , cols )
      if( ttmp != None ): time = np.concatenate( (tx , ttmp) )
      vx = np.concatenate( (vx , vtmp) )
    
    ttmp = None; vtmp = None
  
  print('... done!')

  return tx, vx

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def samplingFrequency(time, sfreqUser ):
  if( time is not None ):
    dt   = (time[1]-time[0])
    sfreq = 1./dt
  else:
    sfreq = sfreqUser
    
  return sfreq

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def timeInterval( time, v, sfreq ):
  if( time != None ):
    deltaT = (time[-1]-time[0])
  else:
    deltaT = float( len(v) )/ sfreq
    
  return deltaT

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def applyTapering( v , dtime , sfreq ):
  
  # m = 2^k    
  k = int( np.log( dtime * sfreq )/ np.log(2.) )
  m = 2**k
  if( m > len(v) ):
    print(' m = {} vs. len(v) = {} '.format(m,len(v)))
    sys.exit(' Error! Exiting.')
  
  idx = np.arange(m)
  hw  = np.hamming(m) # Hamming window, dealing with edges in fft (or BellTaper, chebwin...)
  
  # Kaimal & Finnigan 1994 p.267, used for compensating for windowing, multiply the power by this
  Hamming_normfact = 2.52
  # Linear detrending
  vd = scs.detrend(v[idx]); v = None
  
  # Tapering with a chosen window
  vw = hw*vd; vd = None  # the time series may have aquired a non-zero mean after windowing   
  
  return scs.detrend(vw)

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def evalSpectra( v, sfreq, normalizeOn=False ):
  
  # Power is the square of the amplitudes
  P = np.abs( np.fft.fft(v) )**2

  ''' 
  Select first half and let 1st point (DC value) be replaced with 2nd point
  and specify the frequency distribution that goes from min(freqs)=2*df, max(freqs)=sfreq/2.
  That is, from the largest wave lengths to the smallest ...
  NOTE: Smallest frequency is 2*[sampling frequency].
  '''
  nv  = len(v)
  ids = np.arange(1,(nv/2)+1); ids[0]=2 
  df = sfreq/float(nv)  # Frequency interval
  freqs = df * ids # Actual frequencies, min(freqs)=2*df, max(freqs)=sfreq/2
  
  E = 2.*P[ids]/nv**2    # Power spectral energy (e.g. Stull p.313) [units of variance]
  S = E/df
  P = P[ids]
  '''
  Power spectral density [variance per frequency interval]
  Integration over all frequencies gives variance. 
  '''
  if( normalizeOn ):
    # Variance (of windowed v) -> the power does not need to be multiplied by Hamming_normfact
    var_v = np.std(v)**2
    S *= freqs    # freqs/df * E  
    S /= var_v    #  Normalized power spectra [dimensionless]
  
  
  # Return Power, Power spectral energy, and power spectral density
  return P, E, S, freqs    

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def frequencyBins( Q , freqs, Nbins ):
  
  minLogFz = np.log(np.min(freqs))
  maxLogFz = np.log(np.max(freqs))
  
  df  = (maxLogFz-minLogFz)/float(Nbins)
  logBins  = np.linspace( minLogFz, maxLogFz, Nbins, endpoint=True )
  freqBins = np.exp(logBins); logBins = None
  
  fbin = np.zeros(Nbins); fbin[:] = None   # None -> Nan
  Qbin = np.zeros(Nbins); Qbin[:] = None
  for i in xrange(Nbins-1):
    ieff = (freqs>freqBins[i]) * (freqs<=freqBins[i+1])
    if( any(ieff) ):
      #print('{}: N of ieff = {} '.format(i,np.count_nonzero(ieff)))
      Qbin[i] = np.nanmean( Q[ ieff ] )
    fbin[i] = ( freqBins[i]+freqBins[i+1] )/2.

  return Qbin, fbin

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def spectraAnalysis(fig, v , time, vName, Nbins, mode, normalize=False ):
      
  nterms  = np.shape(v) # Number of variables in the file, number of entries in the data. 
  print(' Number of terms in the {} data series, N = {}'.format(vName, nterms))

  # Determine the sampling frequency.
  samplingFreq   = samplingFrequency( time, None )
  print(' sampling frequency = {}'.format(samplingFreq))

  deltaT = (time[-1]-time[0])
  vw     = applyTapering( v , deltaT , samplingFreq )

  # Evaluate, Power (P), power spectral energy (E), and power spectral density (S).
  P, E, S, freqs = evalSpectra( vw, samplingFreq, normalize )

  
  if(mode == 'S'): 
    Sbin, fbin = frequencyBins( S , freqs, Nbins )
    Vbin = Sbin
  if(mode == 'E'): 
    Ebin, fbin = frequencyBins( E , freqs, Nbins )
    Vbin = Ebin
  if(mode == 'P'): 
    Pbin, fbin = frequencyBins( P , freqs, Nbins )
    Vbin = Pbin

  Refbin = 1.E-1 * fbin[Nbins/2:]**(-5./3.)


  if( mode == 'S' ):
    if( normalize ):
      labelStr = "{}, normalized power spectral density".format(vName)
      plotStr  = [ "Normalized power spectral density","Frequency [Hz]","f*S/$\sigma^2$ "] 
    else:
      labelStr = "{}, power spectral density: $\Phi(f)$".format(vName)
      plotStr  = ["Power spectral density" ,"Frequency [Hz]","$\Phi(f)$"] 
  elif( mode == 'E' ):
    labelStr = "{}, energy spectrum: E(f)".format(vName)
    plotStr  = ["Energy spectrum","Frequency [Hz]","E(f)"]
  else:
    labelStr = "{}, power spectrum: P(f)".format(vName)
    plotStr  = ["Power spectrum","Frequency [Hz]","P(f)"]

    
  if( fig is None ):
    fig = plt.figure(num=1, figsize=(12.,10.))
    fig = addToPlot(fig,fbin[Nbins/2:],np.nanmean(Sbin)*Refbin,"Model -5/3 curve",plotStr,logOn=True)
  fig = addToPlot(fig, fbin, Vbin, labelStr, plotStr, logOn=True)
  
  return fig

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*



# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*