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
#==========================================================#

def dataFromDict( keyStr, dataDict, allowNone=True ):
  data = None
  if( keyStr in dataDict.keys() ):
    data = dataDict[keyStr]
  elif( not allowNone ):
    sys.exit(' Error in dataFromDict: {} not found. Exiting ...'.format(keyStr))
  else:
    pass
  
  return data 

#==========================================================#
def sensibleIds( ijk, x, y, z ):
  '''
  Check whether the chosen i,j,k indices make sense for given x, y, z coords.
  '''
  ijk[0] = np.minimum( ijk[0] , len(x)-1 ); ijk[0] = np.maximum( ijk[0], 0 )
  ijk[1] = np.minimum( ijk[1] , len(y)-1 ); ijk[1] = np.maximum( ijk[1], 0 )
  ijk[2] = np.minimum( ijk[2] , len(z)-1 ); ijk[2] = np.maximum( ijk[2], 0 )

  return ijk

#==========================================================#
def groundOffset( vx ):
  '''
  Find ground offset (in z-direction) for given velocity array vx(t, z, y, x) 
  '''
  k_offset = 0
  while 1:
    idNz = (vx[:,k_offset,1,1] > 0.)
    if( any( idNz ) ):
      break
    else:
      k_offset += 1
      
  return k_offset

#==========================================================#

def quadrantAnalysis( v1, v2, qDict ):
  debug = True
  
  # Extract data from dict. Using dict makes future modifications easy.
  ijk1     = dataFromDict('ijk1',     qDict, allowNone=False )
  ijk2     = dataFromDict('ijk2',     qDict, False )
  nkpoints = dataFromDict('nkpoints', qDict, True )
  npx      = dataFromDict('npixels',  qDict, False )
  axisLim  = dataFromDict('axisLim',  qDict, False )
  holewidth= dataFromDict('holewidth',qDict, False )
  weighted = dataFromDict('weighted',  qDict, True )
  
  # Create arrays for running loops over the selected coordinates.
  iList = np.arange(ijk1[0],ijk2[0]+1)
  jList = np.arange(ijk1[1],ijk2[1]+1)
  kList = np.arange(ijk1[2],ijk2[2]+1)

  '''
  In quadrant analysis context, using a stride in z-direction is usually not desired.
  By default npoints is None.'''
  if( nkpoints is None): stride = 1
  else:                  stride = max( ((kList[-1]-kList[0])/nkpoints)+1 , 2 )


  # Compute the covariance term (for example, u'w')
  v = v1*v2
  
  if( debug ):
    print('min(v1)={}, max(v1)={}, min(v2)={}, max(v2)={}'\
      .format(np.abs(np.min(v1)), np.max(v1), np.abs(np.min(v2)), np.max(v2)))
  
  maxLim = np.abs(axisLim)
  minLim = -1.*maxLim

  # Determine if some grid points are under the ground level.
  k_off = max( groundOffset( v1 ), groundOffset( v2 ) )
  if( k_off > 0 ):
    print(' {}: ground offset (k_off) = {}'.format(filename, k_off))

  x = np.linspace(minLim,maxLim,npx+1)
  y = np.linspace(minLim,maxLim,npx+1)
  dx = (maxLim-minLim)/(npx)
  X,Y = np.meshgrid(x,y)
  Qi  = np.zeros( np.shape(X), float )


  nTot = 0
  nQTot= 0
  nQ1  = 0  # q1: u'(+), w'(+), OutwardInteraction
  nQ2  = 0  # q2: u'(-), w'(+), Ejection
  nQ3  = 0  # q3: u'(-), w'(-), Inward Interaction
  nQ4  = 0  # q4: u'(+), w'(-), Sweep
  
  SQ1 = 0.; SQ2 = 0.; SQ3 = 0.; SQ4 = 0. 
  STot= 0.
  
  cf = 1.
  for i in iList:
    for j in jList:
      for k in kList[::stride]:
        vt  = v[:,k+k_off,j,i]
        
        vt_mean = np.mean( np.abs(vt) )

        v1t = v1[:,k+k_off,j,i] 
        v2t = v2[:,k+k_off,j,i] 
        
        for l in xrange( len(vt) ):
          STot += vt[l]; nTot += 1
          if( np.abs(vt[l]) > (holewidth*vt_mean) ): 
            n = np.minimum( int((v1t[l] - minLim)/dx) , npx )
            n = np.maximum( n , 0 )
            m = np.minimum( int((v2t[l] - minLim)/dx) , npx )
            m = np.maximum( m, 0 )
            Qi[m,n] += 1.; nQTot += 1
          
            if( v1t[l] > 0. and v2t[l] > 0. ):
              nQ1 += 1; SQ1 += vt[l]
            elif( v1t[l] < 0. and v2t[l] > 0. ):
              nQ2 += 1; SQ2 += vt[l]  # Ejection
            elif( v1t[l] > 0. and v2t[l] < 0. ):
              nQ4 += 1; SQ4 += vt[l]  # Sweep
            else:
              nQ3 += 1; SQ3 += vt[l]

  v = None; v1 = None; v2 = None
  Qi /= (np.float(nQTot)*dx**2)  # Obtain the PDF
  
  if( weighted ):
    Qi *= np.abs(X*Y)
  
  STot /= np.float(nTot)
  SQ1 /= np.float(nQTot); SQ2 /= np.float(nQTot)
  SQ3 /= np.float(nQTot); SQ4 /= np.float(nQTot)
  
  
  # Assemble the result dict 
  rDict = dict()
  rDict['STot'] = STot
  rDict['nQ1']  = nQ1; rDict['SQ1'] = SQ1
  rDict['nQ2']  = nQ2; rDict['SQ2'] = SQ2   # Ejection
  rDict['nQ3']  = nQ3; rDict['SQ3'] = SQ3
  rDict['nQ4']  = nQ4; rDict['SQ4'] = SQ4   # Sweep
  #rDict['klims']= np.array([ kList[0], kList[-1] ])
  
  return Qi, X, Y, rDict
