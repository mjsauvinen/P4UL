from osgeo import gdal
from gdalconst import *
import operator
import scipy.ndimage as sn 
import numpy as np
import sys
''' 
Description:


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''

gdal.UseExceptions()

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
'''
Global dict for determining steps
Note: These are steps in [row/col] sense.
At level 4 when the alphabets A-H are used
the E-direction is split into 4 subtiles instead of 
the usual 2. This requires some care with the step
lengths.
'''

EStepUtmDict = {1:0,2:0,3:1,4:1,\
      'A':0,'B':0,'C':1,'D':1,\
      'E':2,'F':2,'G':3,'H':3}
NStepUtmDict = {1:-1,3:-1,2:0,4:0,\
      'A':-1,'B':0,'C':-1,'D':0,\
      'E':-1,'F':0,'G':-1,'H':0} 

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def UtmTileDims( level=1 ):
  # Tile dimensions [N,E] (or [row,col]) in meters.
  LN = 96; LE = 192
  if( level >= 4 ): LE/=2
  return np.array([LN, LE])*1000

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def openGeoTiff( filename ):
  dataset = gdal.Open( filename, GA_ReadOnly)
  if (dataset is None):
    sys.exit(' Unable to open file: {}. Exiting.'.format(filename))
    
  return dataset
  
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def getGeoTransform(ds):
  '''
  Returns the top left origin and the pixel resolution.
  The pixel resolution is in meters and includes a sign.
  The N/row-direction usually has negative dPxl value.
  '''
  try:
    gt = ds.GetGeoTransform()
    XO_TL = np.array([gt[3], gt[0]],int) # [N,E] top left origin
    dPxl  = np.array([gt[5], gt[1]],int) # [N,E]/[row/col] pixel size.
  except:
    print('Could not obtain geo transformation info.')
    XO_TL = None; dPxl = None 
    
  return XO_TL, dPxl

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


def numberOfRasterBands( ds , printOn=False):
  nb = ds.RasterCount 
  if( printOn ):
    print(" RASTER BAND COUNT: {} ".format(nb))

  return nb

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def selectBand( nBands , promtUser=True, defInt=1):
  defInt = int(defInt)
  if( nBands > 1 and promptUser ):
    ib = int(input('Select Band (1-{}): ib = '.format(nBands)))
    ib = min(abs(ib), nBands)
  else:
    ib = defInt

  return ib

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def getRasterBand( dataset, iband ):
  ibandMax = dataset.RasterCount
  if( iband > ibandMax ):
    sys.exit('Error in extractRasterBand: iband > ibandMax')
  
  try:
    rBand = dataset.GetRasterBand(iband)
  except:
    print('Error in extracting band {}'.format(iband))
    print('{}'.format(e))
    sys.exit(1)

  return rBand  

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def printRasterBandStatistics(rBand):
  try:
    Rmin,Rmax,Rmean,Rstd = rBand.GetStatistics(True,True)
    Rscale = rBand.GetScale()
    Runit =  rBand.GetUnitType()
    print('''
    [ STATS ]:  
              Minimum = {0:6.3f}
              Maximum = {1:6.3f}
              Mean    = {2:6.3f}
              Std     = {3:6.3f}
              Scale   = {4}
              Unit    = {5}
    '''.format(Rmin,Rmax,Rmean,Rstd,Rscale,Runit))
  except Exception as e:
    print('{}'.format(e))

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def UtmReference( LocCode ):
  RefCode  = 'L4'
  RefChars  = list(RefCode)
  LocChars  = list(LocCode.upper()) # Make sure it's upper case.
  
  RefCoord = np.array([6666, 308])*1000 # L4 Bottom left: etrs-tm35fin [N,E]
  dNE = UtmTileDims()
  

  # Have to perform the char operation by an own function.
  # ndN must be appended with +1 because we want the Top left origin.
  ndN = subtractUtmChars(LocChars[0], RefChars[0])+1 # Letters go N-S (y)
  ndE = int(LocChars[1]) - int(RefChars[1]) # Number go E-W (x)
  nNE  = np.array([ndN,ndE])
  
  return RefCoord+nNE*dNE

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def extractSubTile( rBand, tileCode, XOrg, dPx):
  
  Xtmp = XOrg.copy()    # Make a proper copy of the XOrg.
  nPx  = None # np.array([ rBand.YSize , rBand.XSize ], int)/dPx 
  nPxOffset = None
  
  if( tileCode is not None ):
    tileChars = list(tileCode)
    code = tileChars[0]; code+=tileChars[1]
    pStr = '{}: XOrig = {}, nPx = {} [N,E]'
    print(pStr.format(code,Xtmp,None))
  
    if( len(tileCode) == 2 ): # The user could ask for the entire tile.
      nPx = UtmTileDims()/np.abs(dPx)
    else:  
      for level in range(1,len(tileCode)-1):
        code+=tileChars[level+1]
        Xtmp, nPx = newTileCoords( Xtmp, tileChars, level, dPx )
        print(pStr.format(code,Xtmp,nPx))

    nPxOffset = np.abs( (Xtmp-XOrg)/dPx )
    print(' Number of Offset Pixels = {}'.format(nPxOffset))
  
  Rb = readAsNumpyArray( rBand, nPxOffset, nPx)
  Rdict = {'R' : Rb, 'GlobOrig' : Xtmp}
  
  return Rdict
  
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


def newTileCoords( XOld, tileChars, level, dPx ):
  # First test whether it's a number or an alphabet.
  m = tileChars[1+level] # Min 2, skip 0,1.
  if (level == 4 ):
    if(m not in NStepUtmDict.keys()):
      print(' Letter {} is not acceptable. Only A-H allowed.'.format(m))
      sys.exit(1)
  else:
    m = int(m)
    if(m not in NStepUtmDict.keys()):
      print(' Integer {} is not acceptable. Only 1-4 allowed.'.format(m))
      sys.exit(1)
    
  nSteps = np.array([ NStepUtmDict[m],EStepUtmDict[m] ], int)
  Xdims = UtmTileDims(level)/2**(level)
  dXO   = (nSteps*Xdims).astype(int) 
  XOnew = XOld + dXO
  nPx   = Xdims/np.abs(dPx)     # Number of pixels [N,E] or [row/col]
  
  return XOnew, nPx

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def extractSubTileFromNumpyArray(R, XO, LocCode, level, dPx ):
  
  dPx = np.abs(dPx) # We want positive resolution values only.
  LocChars = list(LocCode)
  if( level > len(LocCode)-2 ):
    return R

  # First test whether it's a number or a char.
  m = LocChars[1+level] # Min 2, skip 0,1.
  if (level == 4 ):
    #this is going to be for char.
    m = subtractUtmChars(m,'A')+1
  else:
    m = int(m)
  
  nSteps = np.array([ NStepUtmDict[m],EStepUtmDict[m] ], int)
  Xdims = UtmTileDims()/2**(level)
  XOnew = XO + nSteps*Xdims
    
  # Compute the number of points extracted from the raster data.
  # dPx is the resolution of raster data.
  nPx = Xdims/dPx   
  iPx1 =  np.abs(nSteps)   *Xdims/dPx  # Starting indecies.
  iPx2 = (np.abs(nSteps)+1)*Xdims/dPx  
  print(' Start, Stop : {}, {} '.format(iPx1, iPx2))
    
    
  return R[iPx1[0]:iPx2[0], iPx1[1]:iPx2[1]], XOnew 

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def readAsNumpyArray( rBand, nPxoff=None, nPx=None ):
  ''' 
  The dXoff and nPx come in as [N,E] arrays.
  Now we need to specify x,y values.
  '''
  if( nPxoff==None and nPx==None ):
    try:
      dat = rBand.ReadAsArray()
    except:
      sys.exit('Error in (2) readAsNumpyArray. Exiting.')
  else:
    try:
      #  ReadAsArray(xoff,yoff, xsize, ysize)
      # Make sure the type is int and not numpy.int64 !
      i1=int(nPxoff[1]); i2=int(nPxoff[0])
      j1=int(nPx[1]);    j2=int(nPx[0])
      dat = rBand.ReadAsArray(i1,i2,j1,j2)
    except:
      sys.exit('Error in (1) readAsNumpyArray. Exiting.')
    
  return dat

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def subtractUtmChars(c2,c1):
  ''' 
  In the Utm tile system 'O' is skipped.
  This routine removes the headache of dealing with 
  the related pains.
  '''
  c2 = c2.upper(); c1.upper()
  v2 = ord(c2); v1 = ord(c1); vO = ord('O')
  diff = v2 - v1
  
  if( v2 > vO and v1 < vO ):
    return diff-1
  else:
    return diff

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

