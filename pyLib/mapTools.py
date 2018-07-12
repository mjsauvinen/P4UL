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
    XO_TL[0]= dtmp['yllcorner']+dtmp['nrows']*dtmp['cellsize']
  else:
    XO_TL[0]= dtmp['ytlcorner']

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
    XO_TL[1]= dtmp['xllcorner']
  else:
    XO_TL[1]= dtmp['xtlcorner']

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
            Rdict = readNumpyZTile(d['name']+'.npz')
            r=Rdict['R']
            Rdict = None   # Throw the rest away.
          M.append(r); r = None

  print(' M.shape = {}'.format(np.shape(M)))

  T = None
  for i in xrange(Mrows):
    c1 = i*Mcols; c2 = (i+1)*Mcols
    print 'c1={}, c2={}'.format(c1,c2)
    if( T is None ):
      T = np.hstack(M[c1:c2])
    else:
      T = np.vstack( (T,np.hstack(M[c1:c2])) )

    print(' np.shape(T) = {}'.format(np.shape(T)))

  M = None
  Rdict = {'R' : T}
  return Rdict

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

def asciiTileToNumpyZ(filename):

  Rdict, idx = readAsciiGridHeader( filename )
  R = readAsciiGrid( filename )

  Rdict['id'] = idx
  Rdict['ytlcorner'] = Rdict['yllcorner'] + Rdict['cellsize']* Rdict['nrows']
  Rdict['xtlcorner'] = Rdict['xllcorner']

  # These are ofter used later.
  Rdict['GlobOrig'] = np.array([ Rdict['ytlcorner'], Rdict['xtlcorner']]) # [N,E]
  Rdict['dPx'] = np.array([ Rdict['cellsize'], Rdict['cellsize'] ])
  Rdict['R'] = R

  saveTileAsNumpyZ( filename.strip('.asc'), Rdict )

  return Rdict

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def readNumpyZGridData( filename, idx=0 ):
  Rdict = readNumpyZTile(filename, dataOnly=True)
  Rxdims=np.array(np.shape(Rdict['R']))
  Rdict['R'] = []
  RxOrig=Rdict['GlobOrig']
  dPx=Rdict['dPx']
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

def saveTileAsNumpyZ( filename, Rdict):
  '''
  The saved npz file doesn't contain Rdict, but separate numpy arrays matching key names.
  Therefore np.load(filename) is equal to the saved Rdict.
  '''
  try:
    np.savez_compressed(filename, **Rdict)
    print(' {} saved successfully!'.format(filename))
  except:
    print(' Error in saving {}.npz in saveTileAsNumpyZ().'.format(filename))

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def readNumpyZTile( filename, dataOnly=False, verbose=True):
  if (verbose):
    print(' Read filename {} '.format(filename))
  # dat must be closed to avoid leaking file descriptors.
  try:
    dat = np.load(filename)
    Rdict = dict(dat)
    dat.close()
  except IOError as e:
    print('Error reading file {0}: {1}'.format(filename, e.strerror))
    sys.exit(e.errno)

  #if(dataOnly):
    #Rdict['R'] = []

  # Backwards compatibility for variable name change.
  if ('XOrig' in Rdict and not('GlobOrig' in Rdict)):
    Rdict['GlobOrig']=Rdict['XOrig'];
  # For some reason dPx arrays were saved as 'dpx' in the past hardcoded versions of saveTileAsNumpyZ.
  if ('dpx' in Rdict and not('dPx' in Rdict)):
    Rdict['dPx']=Rdict['dpx']
  return Rdict

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def readNumpyZTileForMesh( filename ):
  Rdict = readNumpyZTile( filename )
  Rx = Rdict['R']
  Rxdims = np.array(np.shape(Rx))
  RxOrig = Rdict['GlobOrig']
  dPx = Rdict['dPx']
  try:
    gridRot = Rdict['gridRot']
  except:
    gridRot = Rdict['gridRot'] = 0

  # N,E - coords, start from top left.
  # Sometimes the dPx[0] is correctly negative, but sometimes not.
  # Here, we need to make sure it's own sign is irrelevant
  dPN = np.abs(dPx[0]); dPE = np.abs(dPx[1])
  Rdict['rowCoords'] = np.arange(RxOrig[0],(RxOrig[0]-Rxdims[0]*dPN),-dPN) # N
  Rdict['colCoords'] = np.arange(RxOrig[1],(RxOrig[1]+Rxdims[1]*dPE), dPE) # E

  #Nx, Ex = np.meshgrid(ni,ej)

  return Rdict


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

def rotatePoint(pivot, point, angle):
  # Simple 2D rotation matrix
  rotatedPoint = np.zeros(2)
  rotatedPoint[1] = pivot[1] + np.cos(angle) * (point[1] - pivot[1]) - np.sin(angle) * (point[0] - pivot[0])
  rotatedPoint[0] = pivot[0] + np.sin(angle) * (point[1] - pivot[1]) + np.cos(angle) * (point[0] - pivot[0])
  return rotatedPoint

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def entry2Int( ax ):
  try:
    ax = np.mean(np.abs(ax))
  except:
    pass

  return int(ax)


# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def marginIds( Rxdims, Mw ):

  Li = np.zeros(2, int); Ri = Li.copy(); Bi = Li.copy(); Ti = Li.copy()

  Li[0]= 0
  Li[1]= max( int( np.ceil(Mw[0]*Rxdims[1]-1) ), 1 )  # These can never be -1.

  Ri[0]= min( int((1.-Mw[1])*Rxdims[1]+1), Rxdims[1]-1 )
  Ri[1]= Rxdims[1]

  Bi[0]= min( int((1.-Mw[2])*Rxdims[0]+1), Rxdims[0]-1 )
  Bi[1]= Rxdims[0]

  Ti[0]= 0
  Ti[1]= max( int( np.ceil(Mw[3]*Rxdims[0]-1) ), 1 )  # These can never be -1.

  return Li, Ri, Bi, Ti


# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
def applyMargins( Rx, Mw, Mr, Mh ):

  Rxdims = np.shape(Rx)
  if( Mw.count(None) == 0 ):
    print(' Zero (or non-zero) margins: L={}, R={}, B={}, T={}'.format(Mw[0],Mw[1],Mw[2],Mw[3]))
    L12, R12, B12, T12 = marginIds( Rxdims, Mw )
    L1 = L12[0]; L2 = L12[1]
    R1 = R12[0]; R2 = R12[1]
    B1 = B12[0]; B2 = B12[1]
    T1 = T12[0]; T2 = T12[1]
    #print('Margin\nL:{},{},R:{},{},T:{},{},B:{},{}'.format(L1,L2,R1,R2,T1,T2,B1,B2))

    if( not all( L12 == 0 ) ): Rx[:,L1:L2] = Mh[0]
    if( not all( R12 == 0 ) ): Rx[:,R1:R2] = Mh[1]
    if( not all( T12 == 0 ) ): Rx[T1:T2,:] = Mh[3]
    if( not all( B12 == 0 ) ): Rx[B1:B2,:] = Mh[2]


  else:
    L1=0; L2=1
    R1=Rxdims[1]-1; R2=Rxdims[1]
    B1=Rxdims[0]-1; B2=Rxdims[0]
    T1=0; T2=1

  if( Mr.count(None) == 0 ):
    print(' Ramp margins: L={}, R={}, B={}, T={}'.format(Mr[0],Mr[1],Mr[2],Mr[3]))
    dL  = int(Mr[0]*Rxdims[1]); dR = int(Mr[1]*Rxdims[1])
    dB  = int(Mr[2]*Rxdims[0]); dT = int(Mr[3]*Rxdims[0])

    L11 = max(L2-1,0) ; L22 = L2+dL
    R11 = R1-dR       ; R22 = min(R1+1, Rxdims[1])
    B11 = B1-dB       ; B22 = min(B1+1, Rxdims[0])
    T11 = max(T2-1,0) ; T22 = T2+dT
    #print('Ramp\nL:{},{},R:{},{},T:{},{},B:{},{}'.format(L11,L22,R11,R22,T11,T22,B11,B22))


    if( dL != 0 ):
      if( (Mw[0] is None) or (Mw[0] ==0.) ):
        Rx = applyRamp( Rx, L11, L22, 1, 0, Mh )
      else:
        Rx = applyRamp( Rx, L11, L22, 1, 0 )

    if( dR != 0 ):
      if( (Mw[1] is None) or (Mw[1] ==0.) ):
        Rx = applyRamp( Rx, R11, R22, 1, 1, Mh )
      else:
        Rx = applyRamp( Rx, R11, R22, 1, 1 )

    if( dB != 0 ):
      if( (Mw[2] is None) or (Mw[2] ==0.) ):
        Rx = applyRamp( Rx, B11, B22, 0, 1, Mh )
      else:
        Rx = applyRamp( Rx, B11, B22, 0, 1 )

    if( dT != 0 ):
      if( (Mw[3] is None) or (Mw[3] ==0.) ):
        Rx = applyRamp( Rx, T11, T22, 0, 0, Mh )
      else:
        Rx = applyRamp( Rx, T11, T22, 0, 0 )

  return Rx


# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def applyRamp( Rz, L1, L2, LeftRight, End, Mh=None ):
  dL = (L2-L1)
  w = np.arange( L1, L2 ).astype(float)
  w -= np.min(w); w /= np.max(w)
  w *= np.pi    ; w -= (np.pi/2.)
  w = np.sin(w)/2. + 0.5

  if  ( LeftRight and not End ):      # Left
    if( Mh is None ): Rm = Rz[:,L1]
    else:             Rm = Mh[0]
    #
  elif( LeftRight and End ):          # Right
    if( Mh is None ): Rm = Rz[:,L2]
    else:             Rm = Mh[1]
  elif( not LeftRight and End ):      # Bottom
    if( Mh is None ): Rm = Rz[L2,:]
    else:             Rm = Mh[2]
  else:                               # Top
    if( Mh is None ): Rm = Rz[L1,:]
    else:             Rm = Mh[3]



  if( End ):
    w = (1.-w)
  #print ' w = {}, len(w) = {}, len(dL) = {}'.format(w,len(w),dL)
  if( LeftRight ):
    for i in xrange(dL):
      Rz[:,L1+i] = w[i]*Rz[:,L1+i] + (1.-w[i])*Rm
  else: # TopBottom
    for i in xrange(dL):
      Rz[L1+i,:] = w[i]*Rz[L1+i,:] + (1.-w[i])*Rm

  return Rz

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def filterAndScale(Rxo, Rx, filterInfo, sx=1.0, ix=None, jx=None):

  # Check if the indecies are explicitly given.
  inxOn = True
  if( ix is None or jx is None ):
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

def labelRaster(R, maskId=None):
  import scipy.ndimage.measurements as snms
  Rm = np.zeros( R.shape, type(R[0,0]) )
  if( maskId is not None ):
    mIds = list()
    if(   isinstance( maskId, list) ):
      mIds.extend(maskId)
    elif( isinstance( maskId, int ) ):
      mIds.append(maskId)
    else:
      sys.exit(' Error in labelRaster: maskId is not a list or int. It is {}'.format(type(maskId)))

    idx = np.zeros( R.shape , bool )
    for im in mIds:
      idx = np.maximum( idx , (R == im ) )

    Rm[idx] = R[idx]   # Set desired mask values


  Rl, shapeCount = snms.label(Rm) # this might be slow for unfiltered data
  Rm = None
  print(' Found {} shapes from the data.'.format(shapeCount))

  return Rl, shapeCount


# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


def openTifAsNumpy(tifFile):
  from PIL import Image

  im = Image.open(tifFile)
  #im.show()
  a = np.array(im)

  return a

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

def canopyBetaFunction(height,dpz,alpha,beta,lai):
  '''
  Calculate lead area index using beta probability density function
  (Markkanen et al., 2003, BLM 106, 437-459).
  '''
  from scipy.stats import beta as betadist
  z_col=np.arange(0.,height+dpz[2],dpz[2])  # BUG HERE: np.arange(0.,height/dpz[2],dpz[2])
  z_col=np.divide(z_col,height) # z/H
  #print(' z_col (2) = {}'.format(z_col))

  lad_d = betadist.pdf(z_col,alpha,beta)/height
  #print(' lad_d = {}'.format(lad_d))
  lad=np.multiply(lad_d,lai)
  return lad

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def totalArea( Rdims, dx ):
  #Calculate total area of the domain
  Npx  = np.prod( Rdims ) # Number of pixels
  At = Npx*np.abs(np.prod(dx))
  return At

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def frontalAreas( Ri ):
  # Calculate frontal areas of the domain
  Ae = 0.
  for i in xrange( Ri.shape[0] ):
    ce = (Ri[i,1:] > Ri[i,:-1]).astype(float)
    he = (Ri[i,1:]-Ri[i,:-1]); he[(he<4.)] = 0. # Height, clip out non-buildings
    Ae += np.sum( ce * he )

  An = 0.
  for j in xrange( Ri.shape[1] ):
    cn = (Ri[1:,j] > Ri[:-1,j]).astype(float)
    hn = (Ri[1:,j]-Ri[:-1,j]); hn[(hn<4.)] = 0.
    An += np.sum( cn* hn )

  return Ae, An

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def maskMeanValues(Rm, Ri, mlist):
  dims = np.shape(mlist)

  m_mean = np.zeros(dims)
  m_var  = np.zeros(dims)
  m_std  = np.zeros(dims)
  j = 0
  for im in mlist:
    idm = (Rm == im)
    m_mean[j] = np.mean( Ri[idm] )
    m_var[j]  = np.var( Ri[idm] )
    m_std[j]  = np.std( Ri[idm] )
    print(' Mask {} mean, var, std = {:.2f}, {:.2f}, {:.2f} '.format(im, m_mean[j], m_var[j], m_std[j]))
    j += 1

  return m_mean, m_var, m_std

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def planAreaFractions( Ri, mlist ):
  Npx = np.prod( np.array(Ri.shape) )
  r = np.zeros( np.shape(mlist) )
  j = 0
  for im in mlist:
    r[j] = np.count_nonzero( Ri == im )/float( Npx )
    print(' Mask {} plan area fraction = {:.2f} '.format(im, r[j]))
    j += 1

  return r

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
