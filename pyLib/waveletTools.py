from scipy import signal
import numpy as np
import matplotlib.pyplot as plt

class wtDataset:
  """
  The element of the class constructor are 
  -data: time series to analyse
  -t: time array associated to the time series
  -s (optional): : scales to use for the wavelet transform and the scalogram; if None (by default),  a scales arrays is constructed  starting from the given frequency range.
  -f (optionnal): frequencies to use for the wavelet transform and the scalogram; if None (by default), a frequencies arrays is constructed starting from the given scales range.
  -omega0 (optional): the central pulsation of the Morlet/RealMorlet wavelet.

  """
  def __init__(self, data,  t, s=None, f=None, omega0=6.):
    self.data = data
    self.scales = s
    self.t=t
    self.omega0=omega0
    self.freq=f
  
    assert( (f is not None) or (s is not None) ), "Neither scales nor frequencies provided for the wavelet transform"

    if(f is None):
      f1 = self.omega0/(2*np.pi*max(self.scales))
      f2 = self.omega0/(2*np.pi*min(self.scales))
      ns = len(self.scales)
      self.freq=np.linspace(f1,f2,ns)
    if(s is None):
      s1 = self.omega0/(2*np.pi*max(self.freq))
      s2 = self.omega0/(2*np.pi*min(self.freq))
      ns = len(self.freq)
      self.scales=np.linspace(s1,s2,ns)

#==========================================================#

  def getBounds( self ):
    extDict = dict()
    extDict['scales'] = list( [self.t[0], self.t[-1], self.scales[-1], self.scales[0] ] )
    extDict['freq']   = list( [self.t[0], self.t[-1], self.freq[-1],   self.freq[0]   ] )
    
    return extDict

#==========================================================#
  def SigMorletScalogram(self,ttype="complex", plotOn=True):
    """
    Morlet wavelet transform of the signal with respect to scale parameter,
    Arguments:
    ttype (optional): it's the string "complex" by default; it can be set to "real"
    plot (optional): boolean variable which plots the scalogram when set to True; it is True by default
    """
    assert (ttype=="real" or ttype=="complex"), "Transform type must be in string format: \"real\" or \"complex\" "
    
    dt=np.mean( self.t[1:]-self.t[:-1] )

    output = np.zeros( [len(self.scales), len(self.data)] , dtype=complex)
    
    for ind, width in enumerate(self.scales):
      wavelet_data = Morlet(min(width*10./dt, len(self.data)), width,self.omega0, dt)
      output[ind, :] = 1./np.sqrt(width)*dt*\
        signal.convolve(self.data, np.real(wavelet_data), mode='same') + \
          1./np.sqrt(width)*dt*signal.convolve(self.data, np.imag(wavelet_data), mode='same')*1j

    if(ttype=="real"): output=np.real(output)
    
    if(plotOn):
      fig3=plt.figure()
      
      if( ttype =="real" ):
        ax3a=fig3.add_subplot(1,1,1)
        ax3a.set_xlabel('time (s)'); ax3a.set_ylabel('scale (s)')
        fig3.suptitle("Real Morlet wavelet signal scalogram",fontsize=16)
        cxa=ax3a.imshow(output, extent=[self.t[0], self.t[-1], self.scales[-1], self.scales[0] ], cmap='PRGn',\
          aspect='auto', vmax=abs(output).max(), vmin=-abs(output).max())
        plt.colorbar(cxa)
        
      else: # "complex"
        fig3.suptitle("Morlet Wavelet signal scalogram",fontsize=16)
        ax3a=fig3.add_subplot(2,1,1)
        ax3a.set_xlabel('time (s)'); ax3a.set_ylabel('scale (s)')
        ax3a.set_title("Real Part")
        cxa=ax3a.imshow(np.real(output), extent=[self.t[0], self.t[-1], self.scales[-1], self.scales[0] ],\
          cmap='PRGn', aspect='auto', vmax=abs(output).max(), vmin=-abs(output).max())
        plt.colorbar(cxa)
        
        ax3b=fig3.add_subplot(2,1,2)
        ax3b.set_title("Imaginary Part")
        ax3b.set_xlabel('time (s)')
        ax3b.set_ylabel('scale (s)')
        cxb=ax3b.imshow(np.imag(output), extent=[self.t[0], self.t[-1], self.scales[-1], self.scales[0] ],\
          cmap='PRGn', aspect='auto', vmax=abs(output).max(), vmin=-abs(output).max())
        plt.colorbar(cxb)
      
    return output

#==========================================================#

  def SigMorletSpectrogram(self,ttype="complex",plotOn=True):
    """
    RealMorlet wavelet transform of the signal  with respect to frequency, taking into account the set 
    central frequency.
    Arguments:
    ttype (optional): it's the string "complex" by default; it can be set to "real"
    plot (optional): boolean variable which plots the spectrogram when True; it is set to True by default
    N.B.: omega0 must be greater or equal to 5!
    """
    assert (ttype=="real" or ttype=="complex"), "Transform type can be only in a string format: \"real\" or \"complex\" "

    if( self.omega0<5 ): 
      raise ValueError, "invalid omega0 value"
    
    dt=np.mean(self.t[1:]-self.t[:-1])

    output = np.zeros([len(self.freq), len(self.data)], dtype=complex )
    
    for ind, width in enumerate(1./self.freq*self.omega0/(2*np.pi)):
        wavelet_data = Morlet(min(width*10./dt,len(self.data)), width,self.omega0,dt)
        output[ind, :] = \
          1./np.sqrt(width)*dt * signal.convolve(self.data, np.real(wavelet_data), mode='same') + \
          1./np.sqrt(width)*dt * signal.convolve(self.data, np.imag(wavelet_data), mode='same')*1j
   
    if(ttype=="real"): output=np.real(output)
    
    if( plotOn ):
      fig3=plt.figure()
      
      if(ttype=="real"):
        ax3a=fig3.add_subplot(1,1,1)
        fig3.suptitle("RealMorlet wavelet signal spectrogram",fontsize=16) 
        cax=ax3a.imshow(output, extent=[self.t[0], self.t[-1], self.freq[-1], self.freq[0] ],\
          cmap='PRGn', aspect='auto', vmax=abs(output).max(), vmin=-abs(output).max())
        plt.colorbar(cax)
        ax3a.set_xlabel('time (s)')
        ax3a.set_ylabel('frequency(Hz)')
 
      else: # "complex"
        ax3a=fig3.add_subplot(2,1,1)
        fig3.suptitle("Morlet Wavelet signal spectrogram",fontsize=16)
        ax3a.set_title("Real Part")
        cax=ax3a.imshow(np.real(output), extent=[self.t[0], self.t[-1], self.freq[-1], self.freq[0] ],\
          cmap='PRGn', aspect='auto', vmax=abs(output).max(), vmin=-abs(output).max())
        plt.colorbar(cax)
        
        ax3b=fig3.add_subplot(2,1,2)
        ax3a.set_xlabel('time (s)')
        ax3a.set_ylabel('frequency (Hz)')
        ax3b.set_title("Imaginary Part")
        cax=ax3b.imshow(np.imag(output), extent=[self.t[0], self.t[-1], self.freq[-1], self.freq[0] ],\
          cmap='PRGn', aspect='auto', vmax=abs(output).max(), vmin=-abs(output).max())
        plt.colorbar(cax)
        ax3b.set_xlabel('time (s)')
        ax3b.set_ylabel('frequency (Hz)')

    return output

#==========================================================#
  def PowerMorletScalogram(self,ttype="complex", plotOn=True):
    """
    square modulus of the wavelet transform of the signal: it's the power of every component
    Arguments:
    ttype (optional): it's the string "complex" by default; it can be set to "real"
    """
    assert (ttype=="real" or ttype=="complex"), "Transform type can be only in a string format: \"real\" or \"complex\" "
    
    output=abs(self.SigMorletScalogram(ttype=ttype,plotOn=False))**2
    
    
    if( plotOn ):
      fig3=plt.figure()
      ax3a=fig3.add_subplot(1,1,1)
      
      if(ttype=="real"):
        fig3.suptitle("Real Morlet power scalogram",fontsize=16)
      else:
        fig3.suptitle("Morlet power scalogram",fontsize=16) 
      cax=ax3a.imshow(output, extent=[self.t[0], self.t[-1], self.scales[-1], self.scales[0] ],\
        cmap='PRGn', aspect='auto', vmax=abs(output).max(), vmin=-abs(output).max())
      plt.colorbar(cax)
      ax3a.set_xlabel('time (s)')
      ax3a.set_ylabel('scale (s)')

    return output

#==========================================================#
  def PowerMorletSpectrogram(self,ttype="complex", plotOn=True):
    """
    square modulus of the wavelet transform of the signal: it's the power of every component
    Arguments:
    ttype (optional): it's the string "complex" by default; it can be set to "real"
    N.B.: omega0 must be greater or equal to 5!
    """
    assert (ttype=="real" or ttype=="complex"), "Transform type can be only in a string format: \"real\" or \"complex\" "



    if( self.omega0<5 ): 
      raise ValueError, "invalid omega0 value"  

    output=abs(self.SigMorletSpectrogram(ttype=ttype,plot=False))**2
    
    if( plotOn ):
      fig3=plt.figure()
      ax3a=fig3.add_subplot(1,1,1)

      if(ttype=="real"):
        fig3.suptitle("Real Morlet power spectrogram",fontsize=16) 
      else:
        fig3.suptitle("Morlet power spectrogram",fontsize=16)
      
      cax=ax3a.imshow(output, extent=[self.t[0], self.t[-1], self.freq[-1], self.freq[0] ],\
        cmap='PRGn', aspect='auto', vmax=abs(output).max(), vmin=-abs(output).max())
      plt.colorbar(cax)
      ax3a.set_xlabel('time (s)')
      ax3a.set_ylabel('frequency (Hz)')
  
    return output 

#==========================================================#
  def MorletHisto(self,ix, mode="frequency",nbins=30):
    """
    Plots the histogram of the input signal at a specific frequency or scale (via specified index, ix), once mode (i.e., the second argument) is equal to either \"frequency\" or \"scale\"
    Arguments:
    ix
    mode (optional), which is "frequency" by default
    nbin(optional,set to 30 by default): it gives the number of bins in the histogram to plot
    NB It always calculated with respect to the REAL part of the Morlet transform
    """
    assert (mode=="frequency" or mode=="scale"), "Histogram specification can be only in a string format: \"frequency\" or \"scale\" "
 
    fig=plt.figure()
    tStr = "Histogram of fluctuations at {} = {}  ( corresponding to {} = {} ) "
    if(mode=="frequency"):
      titleStr = tStr.format(mode,self.freq[ix], "scale", self.omega0/(2.*np.pi*self.freq[ix])) 
      plt.title(titleStr)
      CWT=self.SigMorletSpectrogram(ttype="real",plotOn=False)
    if(mode=="scale"):
      titleStr = tStr.format(mode,self.scales[ix], "freq", self.omega0/(2.*np.pi*self.scales[ix]))
      plt.title(titleStr)
      CWT=self.SigMorletScalogram(ttype="real",plotOn=False)
      
    plt.hist(CWT[ ix ,:],nbins)
      
#==========================================================#
# Module Functions
#==========================================================#

def Morlet(npoints, width, omega0, dt, ttype="complex"):
  """
  The  central frequency is: Omega0/(2*pi*s)
  The band spans from Omega0/s-3/s and Omega0/s+3/s
  Sampling pulsation must therefore be >2*Omega0/s+3/s
  s here is the scale parameter
  """
  vec= dt*(np.arange(0, npoints) - (npoints - 1.0)/2)
  c=1./(np.sqrt(1.+np.exp(-omega0**2)-2*np.exp(-3./4.*omega0**2)))
  k=np.exp(-omega0**2/2.)
  total = c/np.pi**(1./4.) * \
    np.exp(-vec**2/(2*width**2)) * \
      (np.cos(omega0*vec/width)+(np.sin(omega0*vec/width))*1j-k)

  if(ttype=="real"):
    return np.real(total)
  else:
    return total

#==========================================================#

#==========================================================#