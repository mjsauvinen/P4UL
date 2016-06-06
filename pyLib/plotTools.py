#!/usr/bin/env python
import sys
import glob
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt

# 08.04.2016:  Mona added an option for colorbar bounds to addImagePlot

iCg = 0  # Global integer for color
iMg = 0  # Global integer for markers
gxI = -1     # Global x-index for csv animations
gyLst = []   # Global y-value list for csv animations


# =*=*=*=* FUNCTION DEFINITIONS *=*=*=*=*=*=*=*=*=*=*=*
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
def color_stack(ic=None):
  global iCg
  colorList = ['b', 'r' ,'g', 'c', 'm', 'k']
  if( ic is not None and np.isscalar(ic) ):
    iCg = min( int(ic) , ( len(colorList)-1 ) ) 
  clr = colorList[iCg]
  iCg = min( ( iCg + 1 ), ( len(colorList)-1 ) )
  
  return clr

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def plotBar(fig, xb, yb, labelStr, plotStr=["","",""], wb=0.6, errb=0):
  ax = fig.add_axes( [0.115, 0.075 , 0.85 , 0.81] ) #[left, up, width, height]
  bars=ax.bar(xb,yb,width=wb, label=labelStr, yerr=errb, ecolor='r')
  ax.set_title( plotStr[0], fontsize=22)
  ax.set_xlabel(plotStr[1], fontsize=22)
  ax.set_ylabel(plotStr[2], fontsize=22); ax.grid(True)
  
  return fig
  
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
def addImagePlot( fig, X, labelStr, gridOn=False, limsOn=False ):
  ax = fig.add_axes( [0.1, 0.075 , 0.875 , 0.81] ) #[left, up, width, height]
  im = ax.imshow(X)
  im.set_cmap('rainbow')
  ax.set_title(labelStr)
  ax.grid(gridOn)
  
  if(limsOn): # bounds for a colorbar
    lMin, lMax = raw_input('Enter limits for the colorbar: <min> <max> = ').split()
    im.set_clim([lMin,lMax])
  cbar = fig.colorbar(im)
  
  return fig

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
def addToPlot(fig, x,y,labelStr, plotStr=["","",""], logOn=False):
  '''
  Add variables x,y to a given plot.
  Test whether y has multiple columns --> Require different treatment.
  '''
  ax = fig.add_axes( [0.115, 0.075 , 0.85 , 0.81] ) #[left, up, width, height]
  d = np.size(np.shape(y)) # Test if y has multiple columns

  for i in xrange(d):
    if(d==1):
      yt = y
    else:
      yt = y[:,i]; labelStr+='['+str(i)+']'
    if(logOn):
      lines=ax.loglog(x,yt,'-', linewidth=1.3, label=labelStr)
    else:
      lines=ax.plot(x,yt,'-', linewidth=1.6, label=labelStr)
  ax.set_title( plotStr[0], fontsize=22)
  ax.set_xlabel(plotStr[1], fontsize=22)
  ax.set_ylabel(plotStr[2], fontsize=22); ax.grid(True)  
  return fig
  
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def plotXX( fig, fileStr, logOn, Cx=1., Cy=1. ):
  try:    x = np.loadtxt(fileStr)
  except: x = np.loadtxt(fileStr,delimiter=',')
  ax  = fig.add_axes( [0.15, 0.075 , 0.8 , 0.81] ) #[left, up, width, height], fig.add_subplot(111)

  labelStr = fileStr.split(".")[0]

  # Print each column separately
  Ny = (x.shape[1]-1)
  for i in xrange(Ny):
    if( Ny == 1 ):
      labelXX = labelStr
    else:
      labelXX = labelStr+'['+str(i)+']' 
    
    if( logOn ):
      #lines=ax.loglog(x[:,0],np.abs(x[:,i+1]),'o-', linewidth=1.3 , label=labelXX)
      lines=ax.semilogy(Cx*x[:,0], Cy*np.abs(x[:,i+1]),'-', linewidth=1.1 , label=labelXX)
    else:
      lines=ax.plot(Cx*x[:,0], Cy*x[:,i+1],'-', linewidth=1.1, label=labelXX)

  ax.set_xlabel(" X ")
  ax.set_ylabel(" Y ")
  return fig

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def plotDY( fig, fileStr, dim=3,  revAxes=False ):
  dim = min( dim, 3 ); dim=max(dim , 1)
  x = np.loadtxt(fileStr)
  r = np.zeros( len(x[:,0]), float )
  for i in xrange(dim):
    x0 = np.min( x[:,i] )
    r += (x[:,i]-x0)**2 
    
  d = np.sqrt(r)
  ax  = fig.add_axes( [0.115, 0.075 , 0.85 , 0.81] ) #[left, up, width, height], fig.add_subplot(111)

  labelStr = fileStr.split(".")[0]

  # Print each column separately
  for i in xrange((x.shape[1]-dim)):
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

def plotCXY( fig, fileStr, Cx, Cy ):
  x = np.loadtxt(fileStr)
  ax  = fig.add_axes( [0.115, 0.075 , 0.85 , 0.81] ) #[left, up, width, height], fig.add_subplot(111)

  labelStr = fileStr.split(".")[0]

  # Print each column separately
  for i in xrange((x.shape[1]-1)):
    lines=ax.plot(Cx*x[:,0],Cy*x[:,i+1],'o', markersize=6, linewidth=1.1, label=labelStr+'['+str(i)+']')
    
  return fig

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def plotYX( fig, fileStr, logOn ):
  x = np.loadtxt(fileStr)
  y = x[:,1]
  ax  = fig.add_axes( [0.115, 0.075 , 0.85 , 0.81] ) #[left, up, width, height], fig.add_subplot(111)

  # Print each column separately
  for i in xrange((x.shape[1]-3)):
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
  ax  = fig.add_axes( [0.115, 0.075 , 0.85 , 0.81] ) #[left, up, width, height], fig.add_subplot(111)

  # Print each column separately
  for i in xrange((x.shape[1]-3)):
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
  for i in xrange(len(varList)):
    varList[i]=varList[i].strip("\"")

  x = np.loadtxt(fileStr, delimiter=',', skiprows=1)
  
  if( not globalValues or (globalValues and gxI == -1) ):
    n = 0
    for v in varList:
      print "  => ["+str(n)+"]: "+ v
      n+=1

    try:  
      xI = input(" X [index]  = ")
    except:
      print ' No selection. Exiting program. '
      sys.exit(1)
      
    yLst = []
    try:
      yLst.extend(input(" Y [List] = "))
    except:
      try:
        select=input(" Select All? [1-9]=> Yes, [Empty]=> No: ")
      except:
        print ' Exiting program. '
        sys.exit(1)

    if( len(yLst) == 0 ):
      yLst.extend(range(len(varList)))
      
      
    if( globalValues and gxI == -1 ):
      gxI   = xI   # Store the global values
      gyLst = yLst
      
  else: # (globalValues and gxI /= -1)
    #print ' Utilizing global values '
    xI   = gxI     # Extract the global values
    yLst = gyLst 


  labelStr = fileStr.split(".")[0]

  ax = fig.add_axes( [0.115, 0.075 , 0.85 , 0.81] ) #[left, up, width, height]
  
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
  for i in xrange(len(varList)):
    varList[i]=varList[i].strip("\"")
    varList[i]=varList[i].strip("\""+"\n")  # This is in case the line contain '\n'
    
  Ix = []
  for varStr in varNames:
    try: Ix.append( varList.index(varStr) )#; print "Index List= {}".format(Ix)
    except: None
  
  if (len(Ix) == 0):
    print "None of the variables in {0} were found in {1}".format(varNames,varList)
    print "Exiting program. "
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
    print "None of the variables in {0} were found in {1}".format(varNames,varList)
    print "Exiting program. "
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

def addContourf( X, Y, Q, labelStr, titleStr=" " ):
  Xdims = np.array(X.shape)
  figDims = 12.*(Xdims[::-1].astype(float)/np.max(Xdims))
  plt.figure(figsize=figDims)
  #levels = [-1e-6, -1e-7, 0, 1e-7, 1e-6]
  #CO = plt.contourf(X,Y,Q, levels )
  CO = plt.contourf(X,Y,Q, 20 )
  plt.title( titleStr )
  cbar = plt.colorbar(CO)
  cbar.ax.set_ylabel(labelStr)
  return CO


# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*  
def addScatterPlot(fig, X, Y, C, fc=4 ):
  ax = fig.add_axes( [0.115, 0.075 , 0.85 , 0.81] ) #[left, up, width, height]
  dims = np.array(np.shape(X))/fc
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
    print ' The file must contain (at least) 4 columns: x, y, dx, dy '
    sys.exit(1)
    
  ax  = fig.add_axes( [0.075, 0.075 , 0.85 , 0.85] ) #[left, up, width, height], fig.add_subplot(111)
  
  lx = max(scale, 0.825 )*0.0008
  lx = min( lx, 0.0016 )
  for i in xrange( len(x) ):
    ax.arrow( x[i], y[i], scale*dx[i], scale*dy[i], color=color_stack(ic) , width=lx, \
             head_width=5.85*lx, head_length=2.85*lx, overhang=0.25, fill=fillOn )

  return fig

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def writeXY( x , y , fileName ):
  f = open( fileName ,'w')             #'w' = for writing
  for i in xrange(len(x)):
    f.write("%13.7e \t %13.7e \n" %(x[i], y[i]) )
  print 'Writing file '+fileName
  f.close()

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def wavePlot( fig, fileStr, logOn ):
  x = np.loadtxt(fileStr)
  ax  = fig.add_axes( [0.15, 0.075 , 0.8 , 0.81] ) #[left, up, width, height], fig.add_subplot(111)

  labelStr = fileStr.split(".")[0]

  # Print each column separately
  Ny = (x.shape[1]-1)
  for i in xrange(Ny):
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
  ax=pFig.get_axes()[0]  # Get a handle on the first axes
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
  ax.set_ylabel(yLbl, fontsize=16, fontstyle='normal', fontweight='book', fontname='serif')
  ax.set_xlabel(xLbl, fontsize=16, fontstyle='normal', fontweight='book', fontname='serif')

  return pFig

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def strEntry( questionStr , defaultStr ):
  try:
    oStr = raw_input(str(questionStr))
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
  for i in xrange(x.shape[1]):
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
  ax = fig.add_axes( [0.115, 0.075 , 0.85 , 0.81] ) #[left, up, width, height]
  d = np.size(np.shape(y)) # Test if y has multiple columns

  for i in xrange(d):
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

