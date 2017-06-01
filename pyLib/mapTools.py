import operator
import numpy as np
import sys
''' 
Description:


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def arrangeTileGrid( dictList, fileTypes ):
  coordList = []
  ascii = fileTypes[0]; npz = fileTypes[1]
  
  XO_TL = np.zeros(2)  # Initialize the Top Left Origin.
  
  for d in dictList:
    # The last two indecies are for row / col addresses.
    if( ascii ):
      coordList.append( [d['id'], d['xllcorner'], d['yllcorner'], 0, 0] )
    else: # .npz
      coordList.append( [d['id'], d['xtlcorner'], d['ytlcorner'], 0, 0] )

  # Sort the list according to y-values
  coordListSorted = sorted( coordList, key=operator.itemgetter(2) ) 
  #print ' y-sorted : {} '.format( coordListSorted )


  # Determine the Top Left Origin (y-value).
  ltmp = coordListSorted[-1]  # Take the last entry.
  dtmp = dictList[ltmp[0]]    # Fetch the desired dict. ltmp[0] == id.
  if( ascii ):
    XO_TL[1]= dtmp['yllcorner']+dtmp['nrows']*dtmp['cellsize']
  else:
    XO_TL[1]= dtmp['ytlcorner']

  irow = 0; maxVal = 0.
  for t in reversed(coordListSorted):
    if( t[2] >= maxVal ):
      pass
    else:
      irow+=1  # Change row
    t[3] = irow; maxVal = t[2]
  
  imax = irow+1  # Store the number of rows.

  # Sort the list according to x-values
  coordListSorted = sorted( coordList, key=operator.itemgetter(1) ) 
  #print ' x-sorted : {} '.format( coordListSorted )

  # Determine the Top Left Origin (x-value).
  ltmp = coordListSorted[0]
  dtmp = dictList[ltmp[0]]
  if( ascii ):
    XO_TL[0]= dtmp['xllcorner']
  else:
    XO_TL[0]= dtmp['xtlcorner']

  jcol = 0; minVal = 1.e12
  for t in coordListSorted:
    if( t[1] <= minVal ):
      pass
    else:
      jcol+=1  # Change column
    t[4] = jcol; minVal = t[1]
  
  jmax = jcol+1   # Store the number of columns
  
  ijList = []
  for t in coordListSorted:
    ijList.append( [ t[0], t[3], t[4] ] )# id, irow, jcol
  
  return ijList, XO_TL, imax, jmax

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def compileTileGrid( dictList, ijList, Mrows, Mcols, fileTypes ):
  M = []  # An empty array to start with.
  ascii = fileTypes[0]; npz = fileTypes[1]
  
  for i in xrange(Mrows):
    for j in xrange(Mcols):
      for idTile, irow, jcol in ijList:
      
        if(irow == i and jcol == j ):
          d = dictList[idTile]
          if( ascii ):
            r = readAsciiGrid( d['name']+'.asc' )
          elif( npz ):
            r, a, b, c = readNumpyZTile(d['name']+'.npz')
            a = b = c = None   # Throw the rest away.
          M.append(r); r = None
  
  print(' M.shape = {}'.format(np.shape(M)))

  T = None
  for i in xrange(Mrows):
    c1 = i*Mcols; c2 = (i+1)*Mcols
    print 'c1={}, c2={}'.format(c1,c2)
    if( T == None ):
      T = np.hstack(M[c1:c2])
    else:
      T = np.vstack( (T,np.hstack(M[c1:c2])) )
  
    print(' np.shape(T) = {}'.format(np.shape(T)))

  M = None
  dims = np.shape(T)
  
  return T, np.array(dims)

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def readAsciiGridHeader( filename, idx=0 ):
  fl = open( filename , 'r')
  name = filename.strip('.asc') # Extract the tile name.
  hdict = {'id':idx,'name': name, 'ncols':None,'nrows':None,\
            'xllcorner':None,'yllcorner':None,'cellsize':None,\
            'NODATA_value':None}

  for i in xrange(len(hdict.keys())):
    try:
      s = fl.readline().split()
      hdict[s[0]] = float( s[1] )
    except: 
      print('Unexpected ascii grid header format. Exiting.')
      sys.exit(1)
  
  idx += 1
  fl.close()
  return hdict, idx 
  
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*  

def readNumpyZGridData( filename, idx=0 ):
  Rx, Rxdims, RxOrig, dPx = readNumpyZTile(filename, dataOnly=True)
  name = filename.strip('.npz') # Extract the tile name.
  hdict = {'id':idx,'name': name, 'ncols':Rxdims[1],'nrows':Rxdims[0],\
           'xtlcorner':RxOrig[1],'ytlcorner':RxOrig[0],\
           'cellsize':int(dPx[0]),'NODATA_value':None}
  idx += 1
  return hdict, idx
  
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*  

def resolutionFromDicts( dictList ):
  d1 = dictList[0]
  dPxRef = d1['cellsize']
  
  for d in dictList:
    dPx = d['cellsize']
    if( dPx != dPxRef ):
      print 'The tile resolutions do not match. Exiting.'
      sys.exit(1)
      
  return np.array([dPxRef,dPxRef])

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*  

def readAsciiGrid( filename ):
  try:
    rx = np.loadtxt( filename, skiprows=6 ) # Note! skiprows=6.
    print(' File {} read successfully.'.format(filename))
  except:
    print(' Could not read ascii grid file: {}. Exiting.'.format(filename))
    sys.exit(1)

  return rx

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def saveTileAsNumpyZ( filename, Rx, Rxdims, RxOrig, dPx):
  '''
  The saved .npz file contains R,dims,dpx and XOrig
  '''
  try:
    np.savez_compressed(filename, R=Rx, dims=Rxdims, dpx=dPx, XOrig=RxOrig)
    print(' {}.npz saved successfully!'.format(filename))
  except:
    print(' Error in saving {}.npz in saveTileAsNumpyZ().'.format(filename))
    
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def readNumpyZTile( filename, dataOnly=False ):
  print(' Read filename {} '.format(filename))
  dat = np.load(filename)
  if(dataOnly):
    Rx = []
  else:
    Rx = dat['R']
    
  Rxdims = dat['dims']; RxOrig=dat['XOrig']; dPx=dat['dpx']
  dat.close()
  
  return Rx, Rxdims, RxOrig, dPx

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def readNumpyZTileForMesh( filename ):
  Rx, Rxdims, RxOrig, dPx = readNumpyZTile( filename )
  
  # N,E - coords, start from top left.
  rowCoords = np.arange(RxOrig[0],(RxOrig[0]-Rxdims[0]*dPx[0]),-dPx[0]) # N
  colCoords = np.arange(RxOrig[1],(RxOrig[1]+Rxdims[1]*dPx[1]), dPx[1]) # E
  
  #Nx, Ex = np.meshgrid(ni,ej)
  
  return Rx, rowCoords, colCoords, Rxdims, RxOrig, dPx 
  

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def rotateGridAroundPivot( X, Y, xp, yp, theta, deg=True ):
  
  if( deg ):
    theta = theta * (np.pi/180.)
  
  CtX = np.array([ np.cos(theta), -np.sin(theta)  ])
  CtY = np.array([ np.sin(theta) , np.cos(theta) ])
  #print ' CtX = {} , CtY = {} '.format(CtX, CtY) 
  
  Mdims = np.shape(X)
  XR = np.zeros( Mdims, float )
  YR = np.zeros( Mdims, float )

  
  for i in xrange( Mdims[0] ):
    XR[i,:] = xp + (X[i,:]-xp)*CtX[0] + (Y[i,:]-yp)*CtX[1] # E-X 
    YR[i,:] = yp + (X[i,:]-xp)*CtY[0] + (Y[i,:]-yp)*CtY[1] # N-Y 

  return XR, YR


# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def entry2Int( ax ):
  try:
    ax = np.mean(np.abs(ax))
  except:
    pass

  return int(ax)
  
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
def applyMargins( Rx, Mw, Mr, Mh ):
  
  Rxdims = np.shape(Rx)
  if( Mw.count(None) == 0 ):
    print(' Zero (or non-zero) margins: L={}, R={}, B={}, T={}'.format(Mw[0],Mw[1],Mw[2],Mw[3]))
    L1= 0; L2 = max( int(Mw[0]*Rxdims[1]-1), 0 )  # These can never be -1.
    R1= int((1-Mw[1])*Rxdims[1]); R2 = Rxdims[1]
    B1= int((1-Mw[2])*Rxdims[0]); B2 = Rxdims[0] 
    T1= 0; T2 = max( int(Mw[3]*Rxdims[0]-1) , 0 )  # These can never be -1.

    Rx[:,L1:L2] = Mh[0]; Rx[:,R1:R2] = Mh[1]
    Rx[T1:T2,:] = Mh[3]; Rx[B1:B2,:] = Mh[2]
  else:
    L1= L2 = 0.
    R1= R2 = Rxdims[1]-1
    B1= B2 = Rxdims[0]-1 
    T1= T2 = 0

  if( Mr.count(None) == 0 ):
    print(' Ramp margins: L={}, R={}, B={}, T={}'.format(Mr[0],Mr[1],Mr[2],Mr[3]))
    dL  = int(Mr[0]*Rxdims[1]); dR = int(Mr[1]*Rxdims[1])
    dB  = int(Mr[2]*Rxdims[0]); dT = int(Mr[3]*Rxdims[0])
    L11 = L2    ; L22 = L2+dL
    R11 = R1-dR ; R22 = R1
    B11 = B1-dB ; B22 = B1
    T11 = T2    ; T22 = T2+dT
  
    #Rc = Rx.copy()
    Rx = applyRamp( Rx, L11, L22, 1, 0 )
    Rx = applyRamp( Rx, R11, R22, 1, 1 )
    Rx = applyRamp( Rx, B11, B22, 0, 1 )
    Rx = applyRamp( Rx, T11, T22, 0, 0 )

    #Rx -= Rc
    
  return Rx


# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def applyRamp( Rz, L1, L2, LeftRight, End ):
  dL = (L2-L1)
  w = np.arange( L1, L2 ).astype(float)
  w -= np.min(w); w /= np.max(w)
  w *= np.pi    ; w -= (np.pi/2.)
  w = np.sin(w)/2. + 0.5 
  if( End ):
    w = (1.-w)
  #print ' w = {}, len(w) = {}, len(dL) = {}'.format(w,len(w),dL)
  if( LeftRight ):
    for i in xrange(dL):
      Rz[:,L1+i] *= w[i]
  else: # TopBottom
    for i in xrange(dL):
      Rz[L1+i,:] *= w[i]

  return Rz
      
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def filterAndScale(Rxo, Rx, filterInfo, sx=1.0, ix=None, jx=None):
  
  # Check if the indecies are explicitly given.
  inxOn = True
  if( ix==None or jx==None ):
    inxOn = False
    
  if( filterInfo.count(None) == 0):
    if( 'user' in filterInfo[0] ):
      nI = int(filterInfo[1])
      for i in xrange(nI):
        ftmp = raw_input(' Enter <method>, <num> = ').split(',')
        if( i == 0 and inxOn ):  Rxf = applyFilter(Rx[ix,jx], ftmp)
        else:                    Rxf = applyFilter(Rx, ftmp)       
        Rx = Rxf.copy()
      Rx = None
    else:
      if( inxOn ): Rxf = applyFilter(Rx[ix,jx], filterInfo)
      else:        Rxf = applyFilter(Rx, filterInfo)
      Rx = None
    
    Rxo += sx*Rxf

  else:
    if( inxOn ):
      Rxo += sx*Rx[ix,jx]
    else:
      Rxo += sx*Rx

  return Rxo 

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def applyFilter(Rx, filterInfo ):
  import scipy.ndimage as sn # contains the filters
  
  if( 'gauss' in filterInfo[0] ):
    try: 
      Nf = float(filterInfo[1])
    except: 
      print(' Failed to obtain <sigma> for the Gaussian filter. Exiting.')
      sys.exit(1)
  else:
    try: 
      Nf = int(filterInfo[1])
    except: 
      print(' Failed to obtain <size> for the filters. Exiting.')
      sys.exit(1)
  
  if( 'median' in filterInfo[0] ):
    print(' Median {0}x{0} filter applied. '.format(Nf))
    Rf = sn.median_filter(Rx, size=Nf) 
  elif( 'perc' in filterInfo[0] ):
    print(' Percentile 60 {0}x{0} filter applied. '.format(Nf))
    Rf = sn.percentile_filter(Rx, 60, size=Nf)
  elif( 'rank' in filterInfo[0] ):
    print(' Rank 5 {0}x{0} filter applied. '.format(Nf))
    Rf = sn.rank_filter(Rx, 5, size=Nf)
  elif( 'gauss' in filterInfo[0] ):
    print(' Gaussian sigma={} filter applied. '.format(Nf))
    Rf = sn.gaussian_filter(Rx, sigma=Nf)
  elif( 'local' in filterInfo[0] ):
    print(' Local mean {0}x{0} filter applied. '.format(Nf))
    Rf = sn.uniform_filter(Rx, size=Nf)
  elif( 'max' in filterInfo[0] ):
    print('Max {0}x{0} filter applied. '.format(Nf))
    Rf = sn.maximum_filter(Rx, size=Nf)
  else:
    print(' No filter applied. ')
    Rf = Rx
    
  return Rf

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


def openTifAsNumpy(tifFile):
  from PIL import Image
  
  im = Image.open(tifFile)
  #im.show()
  a = np.array(im); a_dims = a.shape
  
  return a, a_dims
  
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def numpyArray2Tif( arr ):
  from PIL import Image
  
  return Image.fromarray( arr )
  
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def farFieldIds( xc, exclx ):
  # Exclude the given percentile of the nearest upstream field.

  exclx = max( 1. , exclx )
  exclx = min( 99., exclx )
  clw   = exclx / 100.
  xth   = (clw)*np.max(xc) + (1.-clw)*np.min(xc)
  idx   = (xc < xth )
  
  return idx

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def farFieldMean( dat, Xcoord, excludeNearest ):
  if( len(dat) != len(Xcoord) ):
    sys.exit(' Error! The data and coord arrays are different length. Exiting ...')
  
  idx = farFieldIds( Xcoord, excludeNearest )
  
  return np.mean( dat[idx] )

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

