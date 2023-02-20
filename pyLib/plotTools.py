#!/usr/bin/env python3
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from utilities import dataFromDict
from  matplotlib.ticker import FormatStrFormatter

# 08.04.2016:  Mona added an option for colorbar bounds to addImagePlot

plt.rc('xtick', labelsize=24); #plt.rc('ytick.major', size=10)
plt.rc('ytick', labelsize=24); #plt.rc('ytick.minor', size=6)
plt.rcParams["font.family"] = "serif"
plt.rcParams["legend.fontsize"] = "large"
#plt.rcParams["font.serif"] = "Utopia"
#plt.rcParams["font.family"] = "monospace"
#plt.rcParams["font.monospace"] = "Courier"
#plt.rcParams["legend.labelspacing"] = 1.


iCg = 0  # Global integer for color
iMg = 0  # Global integer for markers
iLg = 0  # Global integer for line styles
gxI = -1     # Global x-index for csv animations
gyLst = []   # Global y-value list for csv animations

# The available color maps:
cmaps = { 1:'viridis', 2:'inferno', 3:'plasma',  4:'magma',     5:'Blues', 
          6:'BuGn',    7:'BuPu',    8:'GnBu',    9:'Greens',   10:'Greys', 
         11:'Oranges', 12:'OrRd',  13:'PuBu',   14:'PuBuGn',   15:'PuRd', 
         16:'Purples', 17:'RdPu',  18:'afmhot', 19:'autumn', 
         20:'bone',    22:'cool',  23:'copper', 24:'gist_heat', 
         25:'gray',    26:'hot',   27:'pink',   28:'spring',    29:'summer', 
         30:'winter',  31:'Reds',  32:'YlGn',   33:'YlGnBu',    34:'YlOrBr', 
         35:'YlOrRd',  36:'BrBG',  37:'bwr',    38:'coolwarm',  39:'PiYG', 
         40:'PRGn',    41:'PuOr',      42:'RdBu',    43:'RdGy', 44:'RdYlBu', 
         45:'RdYlGn',  46:'Spectral',  47:'seismic', 48:'Accent', 49:'Dark2', 
         50:'Paired',  51:'Pastel1',   52:'Pastel2', 53:'Set1',   54:'Set2', 
         55:'Set3',    56:'gist_earth',57:'terrain', 58:'ocean',  59:'gist_stern',
         60:'brg',     61:'CMRmap',    62:'cubehelix', 63:'gnuplot', 64:'gnuplot2', 
         65:'gist_ncar',66:'nipy_spectral', 67:'jet',  68:'rainbow', 69:'gist_rainbow', 
         70:'hsv',      71:'flag',          72:'prism', 73:'cividis', 74:'Wistia',
         75:'turbo'}


# =*=*=*=* FUNCTION DEFINITIONS *=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def addFigAxes( fig ):
  
  if( len(fig.get_axes()) == 0 ):
    ax = fig.add_axes( [0.115, 0.09 , 0.85 , 0.81] ) #[left, up, width, height]
  else:
    ax = fig.get_axes()[0]

  return ax

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def printDict( D , ncols=3 ):
  i = 0; pStr = str()
  for k, v in D.items():
    i+=1
    # use at least 13 chars to make columns line up
    pStr += ' {}: {:13s} \t'.format(k,v) 
    if( i%ncols == 0 ):
      print(pStr); pStr = str()
  
  # print whatever is left at the end    
  print(pStr+'\n'); pStr = None; i = None  
  
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def setColormap( img ):
  global cmaps
  # Select the desired colormap
  try:
    printDict( cmaps, 3 )
    icmap = int(input(' Enter integer key for the colormap = '))
    try:    nc = int(input(' Number of discrete colors in colormap = '))
    except: nc = None
    cm = plt.get_cmap( cmaps[icmap], nc )
    img.set_cmap(cm)
  except:
    print(' Using default colormap.')
    pass

  return img

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
def setColorbarLims( img, lMax=None, lMin=None ):
  # Specify the bounds in the colorbar
  if( (lMax is None) or (lMin is None) ):
    try:
      mm = input(' Enter limits for colorbar: <min> <max> =')
      lMin,lMax = list( map(float, mm.split()) )
      img.set_clim([lMin,lMax])
    except:
      pass
  else:
    try:
      lMin = float(lMin); lMax = float(lMax)
      img.set_clim([lMin,lMax])
    except:
      pass
    
  return img

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def userColormapSettings( fig, im, Rmax=None, Rmin=None ):
  uticks =None # User-defined ticks. <None> leads to default setting.
  eformat=None
  
  im = setColorbarLims( im )
  im = setColormap( im )
    
  try:
    uticks=list( map(float, input(' Enter ticks separated by comma (empty=default):').split(',')) )
  except:
    uticks=None
  
  if(Rmax is not None):
    if(Rmax<1.e-3): 
      eformat='%.2e'
  
  cb = fig.colorbar(im, ticks=uticks, format=eformat)
  
  return cb

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def random_marker():
  markerList = ['x','s','p','h','d','*','o','+']
  nm = len(markerList)
  im = np.random.random_integers(nm) - 1
  mrk = markerList[im]
  
  return mrk

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
def marker_stack():
  global iMg
  markerList = ['+','s','D','o','h','p','*','x']
  mrk = markerList[ iMg ]
  iMg = min( ( iMg + 1 ), ( len(markerList)-1 ) )

  return mrk

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def linestyle_stack(lm=1, il=None):
  global iLg
  # '-'  : solid line style,    '--': dashed line style
  # '-.' : dash-dot line style, ':' : dotted line style
  
  if( lm == 1 ):
    lstyleList = ['-','--','-.',':']
  elif( lm == 2 ):
    lstyleList = ['-','--'] # ['x','+'] # ['-','--'] #
  elif( lm == 3 ):
    lstyleList = ['o','x','D','s'] # ['x','+']
  else:
    lstyleList = ['-','--','-.',':']

  nlinestyles = len(lstyleList)
  
  if( il is not None and np.isscalar(il) ):
    iLg = min( int(il) , (nlinestyles-1) )
  
  lstyle = lstyleList[iLg] 
  
  iLg += 1
  if( iLg > (nlinestyles-1) ):
    iLg = 0
    
  return lstyle

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
def color_stack(lm=1, ic=None):
  global iCg
  '''
  Brown       '#A52A2A',
  DeepPink    '#FF1493',
  BlueViolet  '#8A2BE2',
  DarkCyan    '#008B8B',
  DarkOrange  '#FF8C00',
  DarkMagenta '#8B008B',
  GoldenRod   '#DAA520',
  SeaGreen    '#2E8B57',
  OrangeRed   '#FF4500',
  SlateBlue   '#6A5ACD'
  '''
  if( lm == 1 ):
    colorList = ['b','r','c','k','#FF8C00','g','#8B008B',\
      '#FF1493','#8A2BE2','#008B8B','m',\
        '#2E8B57','#FF4500','#6A5ACD',\
          '#A52A2A','#DAA520']
  else:
    colorList = ['b','b','r','r','c','c','k','k','#FF8C00','#FF8C00','g','g','#8B008B','#8B008B',\
      '#FF1493','#FF1493','#8A2BE2','#8A2BE2','#008B8B','#008B8B','m','m',\
        '#2E8B57','#2E8B57','#FF4500','#FF4500','#6A5ACD','#6A5ACD',\
          '#A52A2A','#A52A2A','#DAA520','#DAA520']
  
  ncolors = len(colorList)
  
  if( ic is not None and np.isscalar(ic) ):
    iCg = min( int(ic) , ( ncolors-1 ) ) 
  clr = colorList[iCg]
  iCg += 1
  if( iCg > (ncolors-1) ): 
    iCg = 0
  
  return clr

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def plotBar(fig, xb, yb, labelStr, plotStr=["","",""], wb=0.6, errb=0):
  ax = addFigAxes( fig )
  bars=ax.bar(xb,yb,width=wb, label=labelStr, yerr=errb, ecolor='r')
  ax.set_title( plotStr[0], fontsize=22)
  ax.set_xlabel(plotStr[1], fontsize=22)
  ax.set_ylabel(plotStr[2], fontsize=22); ax.grid(True)
  
  return fig
  
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def addImagePlot( fig, R, titleStr, gridOn=False, limsOn=False, plotNests=False):
  global cmaps
  ax = addFigAxes( fig )
  im = ax.imshow(np.real(R), aspect='auto')       

  while (plotNests):
    try:
      nestParams=list( map(float, input(' Please enter nest location (top left x, top left y, width, height).\n' 
                                          ' Leave empty to continue plotting.\n').split(',')) )
      annotation=str(input(' Please enter annotation for nest.\n'))
    except:
      break
    try: 
      nesti = patches.Rectangle((nestParams[0],nestParams[1]),nestParams[2],nestParams[3], linewidth=1, edgecolor='r', 
                                  facecolor='none')
      ax.add_patch(nesti)
      ax.annotate(annotation,(nestParams[0],nestParams[1]),textcoords='offset pixels',xytext=(4,-18),color='r',size='medium')
    except:
      print(' Nest drawing failed.')

  ax.set_title(titleStr)
  ax.grid(gridOn)
  
  if(limsOn):
    cbar = userColormapSettings( fig, im, np.nanmax(R) )
  else:
    minval  = np.nanmin(R); maxval = np.nanmax(R)
    minSign = np.sign( minval )
    maxSign = np.sign( maxval )
    vmin = min( np.abs(minval), np.abs(maxval) )
    vmax = max( np.abs(minval), np.abs(maxval) )
    
    if( vmax/(vmin+1.E-5) < 1.5 ):
      vmax *= maxSign; vmin = minSign * vmax 
    else:
      vmax *= maxSign; vmin *= minSign
    
    im   = setColorbarLims( im, vmax, vmin )
    cbar = fig.colorbar(im)
  
  return fig

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def addImagePlotDict(fig, RDict ):
  global cmaps 
  
  R    = dataFromDict('R',      RDict, allowNone=False)
  ex   = dataFromDict('extent', RDict, allowNone=True)
  ttl  = dataFromDict('title',  RDict, allowNone=True)
  xlbl = dataFromDict('xlabel', RDict, allowNone=True)
  ylbl = dataFromDict('ylabel', RDict, allowNone=True)
  gOn  = dataFromDict('gridOn', RDict, allowNone=False)
  lOn  = dataFromDict('limsOn', RDict, allowNone=False)
  cm   = dataFromDict('cmap',   RDict, allowNone=True)
  orig = dataFromDict('origin', RDict, allowNone=True)
  
  ax = addFigAxes( fig )
  im = ax.imshow(np.real(R), origin=orig, extent=ex, aspect='auto', cmap=cm)
  
  ax.set_title(ttl); ax.set_xlabel(xlbl); ax.set_ylabel(ylbl)
  ax.grid(gOn)
  
  if(lOn):
    cbar = userColormapSettings( fig, im, np.nanmax(R), np.nanmin(R) )
  else:
    cbar = fig.colorbar(im)
  
  return fig

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
def addToPlot(fig, x,y,labelStr, plotStr=["","",""], logOn=False):
  '''
  Add variables x,y to a given plot.
  Test whether y has multiple columns --> Require different treatment.
  '''
  ax = addFigAxes( fig )
  
  d = np.size(np.shape(y)) # Test if y has multiple columns

  for i in range(d):
    if(d==1):
      yt = y
    else:
      yt = y[:,i]; labelStr+='['+str(i)+']'
    if(logOn):
      lines=ax.semilogy(x,yt,linestyle_stack(1, None), linewidth=2.0, label=labelStr)
      #lines=ax.loglog(x,yt,'-', linewidth=1.3, label=labelStr)
    else:
      lines=ax.plot(x,yt,linestyle_stack(1, None), linewidth=2.0, label=labelStr)
  ax.set_title( plotStr[0], fontsize=22)
  ax.set_xlabel(plotStr[1], fontsize=28)
  ax.set_ylabel(plotStr[2], fontsize=28); ax.grid(True)  
  return fig
  
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def plotXX( fig, pDict, ax=None ):
  global iCg, iMg, iLg
  fileStr = dataFromDict('filename', pDict, allowNone=False)
  logOn   = dataFromDict('logOn',    pDict, allowNone=False)
  llogOn  = dataFromDict('llogOn',   pDict, allowNone=False)
  Cx      = dataFromDict('Cx',       pDict, allowNone=False)
  Cy      = dataFromDict('Cy',       pDict, allowNone=False)
  revAxes = dataFromDict('revAxes',  pDict, allowNone=False)
  linemode = dataFromDict('lm',      pDict, allowNone=False)
  linewidth= dataFromDict('lw',      pDict, allowNone=False)
  ylims   = dataFromDict('ylims',    pDict, allowNone=True)
  xlims   = dataFromDict('xlims',    pDict, allowNone=True)
  reset   = dataFromDict('reset',    pDict, allowNone=True)
  
  try:    x = np.loadtxt(fileStr)
  except: x = np.loadtxt(fileStr,delimiter=',')
  
  if( ax is None ):
    ax = addFigAxes( fig )

  # Reset global integer for color, marker and linestyle.
  if( reset ):
    iCg = 0; iMg = 0; iLg = 0 


  labelStr = labelString( fileStr )
  #lStr = fileStr.rsplit(".", 1)[0]  # Remove the ".dat" 
  #rStr = lStr.rsplit("_")[-1]
  #tStr = lStr.split("/", 2)
  #if( tStr[0] is "." ):
  #  lStr = tStr[1]
  #else:
  #  lStr = tStr[0]
  
  #labelStr = lStr+"_"+rStr
  

  # Print each column separately
  amax = 0.
  Ny = (x.shape[1]-1)
  
  for i in range(Ny):
    if( Ny == 1 ):
      labelXX = labelStr
    else:
      labelXX = labelStr+'['+str(i)+']' 
    
    if( revAxes ):
      yp = Cy*x[:,0];  xp = Cx*x[:,i+1]; dp = xp
    else:
      xp = Cx*x[:,0];  yp = Cy*x[:,i+1]; dp = yp
    
    if( logOn or llogOn ):
      if( revAxes ): 
        xp = np.abs( xp )
        plotf = ax.semilogx 
      else:
        yp = np.abs( yp )
        plotf = ax.semilogy
        
      if( llogOn ):
        plotf = ax.loglog
        
    else:
      plotf = ax.plot
    
    
    lines = plotf( xp, yp, \
      linestyle_stack(lm=linemode), linewidth=linewidth, \
        label=labelXX, color=color_stack(lm=linemode))
    
    lmax = np.abs(np.max(dp))  # Local maximum
    if( lmax > amax ): amax = lmax
    
  #if( amax <5.e-4 and revAxes): 
  #  if( revAxes ): ax.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
  #  else:          ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
  
  ax.set_ybound(lower=ylims[0], upper=ylims[1] )
  ax.set_xbound(lower=xlims[0], upper=xlims[1] )
  
  ax.set_xlabel(" ")
  ax.set_ylabel(" ")
  return fig

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def ciDataFromFile( filename ):
  
  try:    x = np.loadtxt(filename)
  except: x = np.loadtxt(filename,delimiter=',')
  
  nrows, ncols = x.shape
  #print(' nrows, ncols = {}, {}'.format(nrows,ncols))
  if( ncols > 3 ):
    # Copy values and clear memory
    d = x[:,0]; v = x[:,1]; vl = x[:,2]; vu = x[:,3]
  elif( ncols == 2 ):
    d = x[:,0]; v = x[:,1]; vl = x[:,1]; vu = x[:,1]
  else:
    msg = '''
    Error! ncols has a strange value {}. 
    The data must be in [x, v, v_lower, v_upper, (possibly something else)] format.
    Or alternatively [x,v] format in which case no confidence intervals will be present.
    Exiting...'''.format( ncols )
    sys.exit(msg)
  
  # clear memory
  x = None 
  
  return d, v, vl, vu

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
def ciScaleVals( d, v, vl, vu, Cx, Cy, revAxes ):
  if( revAxes ):
    xx  = Cx*v
    vl *= Cx; vu *= Cx
    d  *= Cy
    yy  = d
    
  else:
    yy  = Cy*v
    vl *= Cy; vu *= Cy
    d  *= Cx
    xx  = d
    
  return d, xx, yy, vl, vu

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def ciDiffVals( x1, y1, v1l, v1u, x2, y2, v2l, v2u, revAxes ):
  
  # Because it is possible that all v2u > v1u (or v2u < v1u) we have to prepare for that.
  id2 = (v2u>v1u)
  id1 = ~id2
  
  id2On=False; id1On=False
  if( np.count_nonzero(id2) > 0 ):
    id2On = True
    v1mu = np.abs( np.mean(v1u[id2]) )
    
  if( np.count_nonzero(id1) > 0 ): 
    id1On = True
    v1ml = np.abs( np.mean(v1l[id1]) )
  
  if( revAxes ):
    #if( id2On ): x1[id2] =np.maximum( ((v2l[id2]-v1u[id2])/v1mu)*100., 0.) # If diff is pos, there's overlap
    #if( id1On ): x1[id1] =np.minimum( ((v2u[id1]-v1l[id1])/v1ml)*100., 0.) # If diff is neg, -- " --
    if( id2On ): x1[id2] =np.maximum( (v2l[id2]-v1u[id2]), 0.) # If diff is pos, there's overlap
    if( id1On ): x1[id1] =np.minimum( (v2u[id1]-v1l[id1]), 0.) # If diff is neg, -- " --
    y1 = 0.5*( y1 + y2 )
    dm = np.mean( np.abs(x1) )
  else:
    #if( id2On ): y1[id2] =np.maximum( ((v2l[id2]-v1u[id2])/v1mu)*100., 0.) # If diff is pos, there's overlap
    #if( id1On ): y1[id1] =np.minimum( ((v2u[id1]-v1l[id1])/v1ml)*100., 0.) # If diff is neg, -- " --  
    if( id2On ): y1[id2] =np.maximum( (v2l[id2]-v1u[id2]), 0.) # If diff is pos, there's overlap
    if( id1On ): y1[id1] =np.minimum( (v2u[id1]-v1l[id1]), 0.) # If diff is neg, -- " --  
    x1 = 0.5*( x1 + x2 )
    dm = np.mean( np.abs(y1) )
    
  return x1, y1, dm

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def labelString( fname ):
  
  ls = fname
  if( "." in ls ): 
    ls = ls.rsplit(".", 1)[0]
  
  if( "/" in ls ):
    sL = ls.split('/')
    
    if( len(sL) > 1 ):
      lL = list(map( len, sL ))
      if(  (lL[0] > 1) and ("." not in sL[0]) ):
        ls = sL[0]
      elif((lL[1] > 1) and ("." not in sL[1]) ):
        ls = sL[1]
      else:
        ls = sL[-1]  

  return ls

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def plotCiXY( fig, pDict ):
  fn      = dataFromDict('filename', pDict, allowNone=False)
  Cx      = dataFromDict('Cx',       pDict, allowNone=True)
  Cy      = dataFromDict('Cy',       pDict, allowNone=True)
  linemode= dataFromDict('lm',       pDict, allowNone=False)
  logOn   = dataFromDict('logOn',    pDict, allowNone=True)
  revAxes = dataFromDict('revAxes',  pDict, allowNone=True)
  ylims   = dataFromDict('ylims',    pDict, allowNone=True)
  xlims   = dataFromDict('xlims',    pDict, allowNone=True)
  
  labelStr = labelString( fn )
  
  if( Cx is None ): Cx = 1.
  if( Cy is None ): Cy = 1.
  
  d, v, v_l, v_u = ciDataFromFile( fn )
  
  ax = addFigAxes( fig )
  
  d, xp, yp, v_l, v_u = ciScaleVals( d, v, v_l, v_u, Cx, Cy, revAxes )
  
  if( revAxes ): xlb = 'V(d)'; ylb = 'd'
  else:          ylb = 'V(d)'; xlb = 'd'
  
  
  if( logOn ):
    if( revAxes ):
      plotf  = ax.semilogx
      fillbf = ax.fill_betweenx
    else:
      plotf = ax.semilogy
      fillbf= ax.fill_between
  else:
    plotf = ax.plot
    if( revAxes ):
      fillbf = ax.fill_betweenx
    else:
      fillbf = ax.fill_between
  
  lines = plotf( xp, yp, linestyle_stack(lm=linemode), lw=3., \
    label=labelStr, color=color_stack(lm=linemode))
  linef = fillbf( d, v_u, v_l, facecolor='grey', edgecolor='black', alpha=0.25)
  ax.set_ybound(lower=ylims[0], upper=ylims[1] )
  ax.set_xbound(lower=xlims[0], upper=xlims[1] )
  ax.set_xlabel(xlb)
  ax.set_ylabel(ylb)

  return fig

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def plotCiDiffXY( fig, pDict ):  
  f1      = dataFromDict('fileref',  pDict, allowNone=False)
  fn      = dataFromDict('filename', pDict, allowNone=False)
  Cx      = dataFromDict('Cx',       pDict, allowNone=True)
  Cy      = dataFromDict('Cy',       pDict, allowNone=True)
  linemode= dataFromDict('lm',       pDict, allowNone=False)
  logOn   = dataFromDict('logOn',    pDict, allowNone=True)
  revAxes = dataFromDict('revAxes',  pDict, allowNone=True)
  ylims   = dataFromDict('ylims',    pDict, allowNone=True)
  xlims   = dataFromDict('xlims',    pDict, allowNone=True)
  
  labelStr = labelString( fn )
  
  if( Cx is None ): Cx = 1.
  if( Cy is None ): Cy = 1.
  
  d1, v1, v1_l, v1_u = ciDataFromFile( f1 )
  d2, v2, v2_l, v2_u = ciDataFromFile( fn )
  
  if( d2[-1] != d1[-1] ):
    if( d2[-1] > d1[-1] ):  # Quick and dirty handling for cases when d2[-1] > d1[-1]
      idx = ( d2 <= d1[-1] ) # Take the terms where values match 
      d2 = d2[idx]; v2 = v2[idx]; v2_l = v2_l[idx]; v2_u = v2_u[idx] # Shorten
      
      # Compute the ratio to match the resolutions (roughly)
      r = np.round( (d2[1]-d2[0])/(d1[1]-d1[0]) ).astype(int)
      
      # Use the matching indecies only
      idm = ( np.mod((np.arange(len(d1))+1) , r) == 0 )
      d1 = d1[idm]; v1 = v1[idm]; v1_l = v1_l[idm]; v1_u = v1_u[idm]


  Lm  = min( len(v2), len(v1) )
  d2 = d2[:Lm]; v2 = v2[:Lm]; v2_l = v2_l[:Lm]; v2_u = v2_u[:Lm]
  d1 = d1[:Lm]; v1 = v1[:Lm]; v1_l = v1_l[:Lm]; v1_u = v1_u[:Lm]
  
  d1, x1, y1, v1_l, v1_u = ciScaleVals( d1, v1, v1_l, v1_u, Cx, Cy, revAxes )
  d2, x2, y2, v2_l, v2_u = ciScaleVals( d2, v2, v2_l, v2_u, Cx, Cy, revAxes )
  
  xp, yp, dm = ciDiffVals( x1, y1, v1_l, v1_u, x2, y2, v2_l, v2_u, revAxes )
  
  if( revAxes ): xlb = 'D(d)'; ylb = 'd'
  else:          ylb = 'D(d)'; xlb = 'd'
  
  ax = addFigAxes( fig )
  
  if( logOn ):
    if( revAxes ):
      plotf  = ax.semilogx
      fillbf = ax.fill_betweenx
    else:
      plotf = ax.semilogy
      fillbf= ax.fill_between
  else:
    plotf = ax.plot
    if( revAxes ):
      fillbf = ax.fill_betweenx
    else:
      fillbf = ax.fill_between
  
  lines = plotf( xp, yp, linestyle_stack(lm=linemode), lw=3., \
    label=labelStr+r': $\left< | \Delta | \right>$={:.2g}'.format(dm) , color=color_stack(lm=linemode))
  #label=r': $\left< | \Delta | \right>$={:.2f}'.format(dm) , color=color_stack(lm=linemode))
  
  #linef = fillbf( d, v_u, v_l, facecolor='gray', alpha=0.25)
  
  ax.set_ybound(lower=ylims[0], upper=ylims[1] )
  ax.set_xbound(lower=xlims[0], upper=xlims[1] )
  ax.set_xlabel(xlb)
  ax.set_ylabel(ylb)

  return fig

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
def plotDY( fig, fileStr, dim=3,  revAxes=False ):
  dim = min( dim, 3 ); dim=max(dim , 1)
  x = np.loadtxt(fileStr)
  r = np.zeros( len(x[:,0]), float )
  for i in range(dim):
    x0 = np.min( x[:,i] )
    r += (x[:,i]-x0)**2 
    
  d = np.sqrt(r)
  
  ax = addFigAxes( fig )
  

  labelStr = labelString( fileStr )

  # Print each column separately
  for i in range((x.shape[1]-dim)):
    if( revAxes ):
      lines=ax.plot(x[:,i+dim],d[:],marker=marker_stack(),
                    color=color_stack(), fillstyle='none', ls='None' , label=labelStr+'['+str(i)+']' )
    else:
      lines=ax.plot(d[:],x[:,i+dim],marker=marker_stack(), mew=1.7,
                    color=color_stack(), fillstyle='none', ls='None', label=labelStr+'['+str(i)+']')

  if( revAxes ):
    ax.set_ylabel(" D(X,Y,Z) "); ax.set_xlabel(" F(D) ")
  else:
    ax.set_xlabel(" D(X,Y,Z) "); ax.set_ylabel(" F(D) ")

  return fig

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def plotYX( fig, fileStr, logOn ):
  x = np.loadtxt(fileStr)
  y = x[:,1]
  ax  = addFigAxes( fig )

  # Print each column separately
  for i in range((x.shape[1]-3)):
    if( logOn ):
      lines=ax.semilogy(np.abs(x[:,i+3]), y[:] , linewidth=1.1 , label=fileStr+'_'+str(i))
    else:
      lines=ax.plot(x[:,i+3], y[:], linewidth=1.1, label=fileStr+'_'+str(i) )

  ax.set_xlabel(" F(Y) ")
  ax.set_ylabel(" Y ")
  return fig


# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def fullPlotXY(fig,fileStr,figStr,xlabelStr,ylabelStr,lwidth=1.2,fsize=16,logOn=False):
  x = np.loadtxt(fileStr)
  y = x[:,1]
  ax  = addFigAxes( fig )

  # Print each column separately
  for i in range((x.shape[1]-3)):
    if( logOn ):
      lines=ax.semilogy(np.abs(x[:,i+3]), y[:] , linewidth=lw , label=figStr+'_'+str(i))
    else:
      lines=ax.plot(x[:,i+3], y[:], linewidth=lwidth, label=figStr+'_'+str(i) )

  ax.set_xlabel(xlabelStr, fontsize=fsize)
  ax.set_ylabel(ylabelStr, fontsize=fsize)
  return fig

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def plotCSV( fig, fileStr, revAxes=False, magY=False, globalValues=False ):
  global gxI
  global gyLst
  fl = open( fileStr, 'r' )
  line = fl.readline() # Read first line which contains all variable names as str.
  fl.close()
  varList = line.split(',')
  for i in range(len(varList)):
    varList[i]=varList[i].strip("\"")

  x = np.loadtxt(fileStr, delimiter=',', skiprows=1)
  
  if( not globalValues or (globalValues and gxI == -1) ):
    n = 0
    for v in varList:
      print("  => ["+str(n)+"]: "+ v)
      n+=1

    try:  
      xI = int(input(" X [index]  = "))
    except:
      sys.exit(' No selection. Exiting program. ')
      
    e = input(" Y [List] = ")
    
    if( e == ''): 
      select=input(" Select All? [1-9]=> Yes, [Empty]=> No: ")
      if( select == ''): 
        sys.exit(' Exiting program.')
      else:
        yLst = list(range(len(fileList)))
        
    else:
      try: yLst = list( map( int, e.split(',') ) )
      except: sys.exit(' Bad entry. Exiting program.')
    
    if( globalValues and gxI == -1 ):
      gxI   = xI   # Store the global values
      gyLst = yLst
      
  else: # (globalValues and gxI /= -1)
    #print ' Utilizing global values '
    xI   = gxI     # Extract the global values
    yLst = gyLst 


  labelStr = fileStr.split(".")[0]

  ax = addFigAxes( fig )
  
  if( not magY ):
    yLbl = ""   # Start with empty label
    for yJ in yLst:
      yLbl = yLbl+varList[yJ]+"; "  # Compile label
      
      if( revAxes ):
        lines=ax.plot(x[:,yJ],x[:,xI],'-', markersize=6, linewidth=1.5,  label=labelStr+": "+varList[yJ])
      else:
        lines=ax.plot(x[:,xI],x[:,yJ],'o-', markersize=6, linewidth=1.5, label=labelStr+": "+varList[yJ])
        #writeXY( x[:,xI],x[:,yJ], 'out.dat' )

      
  else:
    yt = np.zeros(len(x[:,0]))
    yLbl = " Mag(y[:]) "   # Set fixed label 
    for yJ in yLst:
      yt += x[:,yJ]**2
    if( revAxes ):
      lines=ax.plot(np.sqrt(yt),x[:,xI],'-', markersize=6, linewidth=1.5,  label=labelStr)
    else:
      lines=ax.plot(x[:,xI],np.sqrt(yt),'o-', markersize=6, linewidth=1.5, label=labelStr)


  if( revAxes ):
    ax.set_ylabel(varList[xI]); ax.set_xlabel(yLbl)   
  else:
    ax.set_xlabel(varList[xI]); ax.set_ylabel(yLbl)

  return fig


# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def extractFromCSV( csvFile, varNames ):
  fl = open( csvFile, 'r' )
  line = fl.readline() # Read first line which contains all variable names as str.
  fl.close()
  varList = line.split(',')
  for i in range(len(varList)):
    varList[i]=varList[i].strip("\"")
    varList[i]=varList[i].strip("\""+"\n")  # This is in case the line contain '\n'
    
  Ix = []
  for varStr in varNames:
    try: Ix.append( varList.index(varStr) )#; print "Index List= {}".format(Ix)
    except: None
  
  if (len(Ix) == 0):
    print("None of the variables in {0} were found in {1}".format(varNames,varList))
    print("Exiting program. ")
    sys.exit(1)
  
  x = np.loadtxt(csvFile, delimiter=',', skiprows=1)
  data = []
  for jcol in Ix:
    data.append( np.array(x[:,jcol]) )

  return np.array(data)

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def extractFromRAW( rawFile, varNames ):
  fl = open( rawFile, 'r' )
  # Read (first or second) line which contains all var names as str.
  while 1:
    line = fl.readline()
    if('#' and 'x' in line):
      break
  fl.close()
  varList = line.split(); varList.remove('#')
  #print varList
  
  Ix = []
  for varStr in varNames:
    try: Ix.append( varList.index(varStr) )#; print "Index List= {}".format(Ix)
    except: None
  
  #print Ix
  if (len(Ix) == 0):
    print("None of the variables in {0} were found in {1}".format(varNames,varList))
    print("Exiting program. ")
    sys.exit(1)
  
  x = np.loadtxt(rawFile)
  data = []
  for jcol in Ix:
    data.append(x[:,jcol])

  return data


# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def addQuiver( X, Y, Ux, Uy , fc,  labelStr, titleStr=" " ):
  plt.figure()
  Q = plt.quiver(X[::fc, ::fc],Y[::fc, ::fc],Ux[::fc, ::fc],Uy[::fc, ::fc],\
  pivot='tail', color='b', units='xy', scale=1.5 )
  #qk = plt.quiverkey(Q, 0.9, 1.05, 1, labelStr, labelpos='E',fontproperties={'weight': 'bold'})

  plt.title(titleStr)
  return Q

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*  

def addContourf( X, Y, Q, CfDict=None ):
  Xdims = np.array(X.shape)
  figDims = 12.*(Xdims[::-1].astype(float)/np.max(Xdims))
  figDims[0] = max( figDims[0] , 9 )
  figDims[1] = max( figDims[1] , 8 )
  #figDims  = (11,11)
  #figDims  = (9,11)
  fig = plt.figure(figsize=figDims)
  #fig, ax = plt.subplots()
  ax = addFigAxes( fig )
  
  # Default values 
  labelStr = ' Q(X,Y) '
  titleStr = ' Title:  Q(X,Y) '
  cm   = None
  vx   = None
  vn   = None
  levels = None
  N = 12
  
  if( CfDict is not None ):
    titleStr = dataFromDict('title',  CfDict, allowNone=False)
    labelStr = dataFromDict('label',  CfDict, allowNone=False)
    cm       = dataFromDict('cmap',   CfDict, allowNone=True )
    N        = dataFromDict('N',      CfDict, allowNone=True )
    vn       = dataFromDict('vmin',   CfDict, allowNone=True )
    vx       = dataFromDict('vmax',   CfDict, allowNone=True )
    levels   = dataFromDict('levels', CfDict, allowNone=True )
    if( N is None ): N = 12
  
  #print(' vmax = {}, vmin = {} '.format(vx,vn))
  #levels = [-1e-6, -1e-7, 0, 1e-7, 1e-6]
  #CO = plt.contourf(X,Y,Q, levels )
  
  if( levels is not None ): CO = ax.contourf(X,Y,Q, levels, cmap=cm, vmin=vn, vmax=vx )
  else:                     CO = ax.contourf(X,Y,Q, N     , cmap=cm, vmin=vn, vmax=vx )
  
  ax.set_title( titleStr )
  
  cbar = fig.colorbar(CO)
  if( vx is not None ): cbar.vmax = vx
  if( vn is not None ): cbar.vmin = vn
  cbar.ax.set_ylabel(labelStr, fontsize=20, fontstyle='normal', fontweight='book', fontname='serif')
  return CO


# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*  
def addScatterPlot(fig, X, Y, C, fc=4 ):
  ax = addFigAxes( fig )
  dims = np.array(np.shape(X))//fc  # NOTE: integer division necessary
  N = np.prod(dims)
  ax.scatter(X[::fc,::fc].reshape(N), Y[::fc,::fc].reshape(N), s=10, \
             c=C[::fc,::fc].reshape(N), marker=',', cmap=plt.cm.rainbow)
  return fig

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*  

def arrow2DPlot( fig, fileStr , scale=1.0, ic=0, fillOn=True ):
  d = np.loadtxt(fileStr)
  labelStr = fileStr.split(".")[0]
  
  try:
    x = d[:,0]; y =d[:,1]; dx = d[:,2]; dy =d[:,3]
  except:  
    print(' The file must contain (at least) 4 columns: x, y, dx, dy ')
    sys.exit(1)
    
  ax  = addFigAxes( fig )
  
  lx = max(scale, 0.825 )*0.0008
  lx = min( lx, 0.0016 )
  for i in range( len(x) ):
    ax.arrow( x[i], y[i], scale*dx[i], scale*dy[i], color=color_stack(ic) , width=lx, \
             head_width=5.85*lx, head_length=2.85*lx, overhang=0.25, fill=fillOn )

  return fig

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def writeXY( x , y , fileName ):
  f = open( fileName ,'w')             #'w' = for writing
  for i in range(len(x)):
    f.write("%13.7e \t %13.7e \n" %(x[i], y[i]) )
  print('Writing file '+fileName)
  f.close()

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def wavePlot( fig, fileStr, logOn ):
  x = np.loadtxt(fileStr)
  ax  = addFigAxes( fig )

  labelStr = fileStr.split(".")[0]

  # Print each column separately
  Ny = (x.shape[1]-1)
  for i in range(Ny):
    if( Ny == 1 ):
      labelXX = labelStr
    else:
      labelXX = labelStr+'['+str(i)+']' 
    
    if( logOn ):
      #lines=ax.loglog(x[:,0],np.abs(x[:,i+1]),'o-', linewidth=1.3 , label=labelXX)
      lines=ax.semilogy(x[:,0],np.abs(x[:,i+1]),'-', linewidth=1.1 , label=labelXX)
    else:
      lines=ax.plot(x[:,0],x[:,i+1],'o', linewidth=1.1, label=labelXX)

  ax.set_xlabel(" X ")
  ax.set_ylabel(" Y ")
  return fig 

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def userLabels( pFig ):
  #print('A) pFig.get_axes()[] ')
  ax=pFig.get_axes()[0]  # Get a handle on the first axes
  #print('B) pFig.get_axes()[] ')
  #pl.rc('text', usetex=True )
  #pl.rc('xtick', labelsize=24)
  #pl.rc('ytick', labelsize=24)

  titleStr = strEntry( " Plot Title = " , " " )
  yLbl     = strEntry( " Y Label = "    , " Y " )
  xLbl     = strEntry( " X Label = "    , " X " )
 
  """
  fontname: [ FONTNAME | 'serif' | 'sans-serif' | 'cursive' | 'fantasy' | 'monospace' ]
  fontsize: [ size in points ]
  fontweight: [ a numeric value in range 0-1000 | 'ultralight' | 'light' | 'normal' | 'regular' | 'book' | 'medium' | 'roman' | 'semibold' |
  'demibold' | 'demi' | 'bold' | 'heavy' | 'extra bold' | 'black' ]
  fontstyle: [ 'normal' | 'italic' | 'oblique']
  """
 
  ax.set_title(titleStr, fontsize=20, fontstyle='normal', fontweight='demibold', fontname='serif')
  ax.set_ylabel(yLbl, fontsize=20, fontstyle='normal', fontweight='book', fontname='serif')
  ax.set_xlabel(xLbl, fontsize=20, fontstyle='normal', fontweight='book', fontname='serif')

  return pFig

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def strEntry( questionStr , defaultStr ):
  try:
    oStr = input(str(questionStr))
  except:
    oStr = str(defaultStr)

  return oStr

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def numEntry( questionStr , defaultValue ):
  try:
    value = input(str(questionStr))
  except:
    value = float(defaultValue)

  return value

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def maxValues( fileStr ):
  x = np.loadtxt(fileStr)
  mv = []
  for i in range(x.shape[1]):
    mv.append(np.max(x[:,i]))

  return mv

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# ADDED MY MONA KURPPA, 2016:
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
def addToPlot_marker(fig, x,y,labelStr, plotStr=["","",""], logOn=False, marker='-'):
  '''
  Add variables x,y to a given plot.
  Test whether y has multiple columns --> Require different treatment.
  
  e.g. marker = '-' or '--' or 'v-'
  '''
  ax = addFigAxes( fig )
  d = np.size(np.shape(y)) # Test if y has multiple columns

  for i in range(d):
    if(d==1):
      yt = y
    else:
      yt = y[:,i]; labelStr+='['+str(i)+']'
    if(logOn):
      lines=ax.loglog(x,yt,marker, linewidth=1.3, label=labelStr)
    else:
      lines=ax.plot(x,yt,marker, linewidth=1.6, label=labelStr)
  ax.set_title( plotStr[0], fontsize=22)
  ax.set_xlabel(plotStr[1], fontsize=22)
  ax.set_ylabel(plotStr[2], fontsize=22); ax.grid(True)  
  return fig

