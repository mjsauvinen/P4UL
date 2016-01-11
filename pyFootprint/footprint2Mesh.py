#!/usr/bin/env python
from utilities import filesFromList, vtkWriteDataStructured2d, vtkWriteHeaderAndGridStructured2d
from utilities import vtkWritePointDataHeader, vtkWritePointDataStructured2D
from plotTools import addContourf, extractFromCSV
from footprintTools import *
from mapTools import readNumpyZTile
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
''' 
Description:


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''
# = # = # = # Function definitions # = # = # = # = # = # = # 

def centralValue( a , tol, nitv=40):
  res = 1.
  ia  = np.ones( np.shape(a) , bool )
  ac_old = 0.
  
  icheck = 0
  while (res > tol and icheck < int(0.001*len(a))):
    
    for i in xrange(nitv):
      imax = np.argmax( a[ia] ); ia[imax] = False
      imin = np.argmin( a[ia] ); ia[imin] = False
      icheck += 2
    
    ac = 0.5*(np.max(a[ia]) + np.min(a[ia]))
    res= abs( ac - ac_old )/max( abs(ac), abs(ac_old) )
    
    ac_old = ac
    #print('res={:e}, i={}'.format(res,icheck))
    
  return ac

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #

def polarPercentileAvg(dat, vleft, vright ):
  return 0.5*(np.percentile(dat,vleft) + np.percentile(dat,vright))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #

def farFieldIds( xc, exclx ):
  # Exclude the given percentile of the nearest upstream field.

  exclx = max( 1. , exclx )
  exclx = min( 99., exclx )
  clw   = exclx / 100.
  xth   = (clw)*np.max(xc) + (1.-clw)*np.min(xc)
  idx   = (xc < xth )
  
  return idx

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #

def farFieldMean( dat, Xcoord, excludeNearest ):
  if( len(dat) != len(Xcoord) ):
    sys.exit(' Error! The data and coord arrays are different length. Exiting ...')
  
  idx = farFieldIds( Xcoord, excludeNearest )
  
  return np.mean( dat[idx] )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #

def meanFromExternal( dat, xe, ye, ze, xl, yl, zl):
  # Define the vertical mean flow from the external data.
  # <#>e: external, <#>l: local
  
  zbmin = max( np.min(zl) , np.min(ze) )
  zbmax = min( np.max(zl) , np.max(ze) )
  if( zbmin > zbmax ): zbmin = zbmax
  
  ybmin = max( np.min(yl) , np.min(ye) )
  ybmax = min( np.max(yl) , np.max(ye) )
  if( ybmin > ybmax ): ybmin = ybmax
  
  xbmin = max( np.min(xl) , np.min(xe) )
  xbmax = min( np.max(xl) , np.max(xe) )
  if( xbmin > xbmax ): xbmin = xbmax
  
  # print(' zbmin = {}, zbmax = {}\n zm = {}'.format(zbmin, zbmax, zm))
  ke = (ze>=zbmin) * (ze<=zbmax)
  je = (ye>=ybmin) * (ye<=ybmax)
  ie = (xe>=xbmin) * (xe<=xbmax)
  
  kz, jy, ix = np.meshgrid(ke, je, ie, indexing='ij', sparse=True)
  idw = (kz * jy * ix).astype(bool)
  
  id0 = ~(dat == 0)   # Indecies where zero occurs.
  idw *= id0  # Take away the zero values.
  #print(' ke = {}'.format(ke) ); print(' je = {}'.format(je) ); print(' ie = {}'.format(ie) )
  return np.mean(dat[idw]), ke, je, ie

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #
# = # = # = # End Function definitions # = # = # = # = # = #

#========================================================== #
parser = argparse.ArgumentParser(prog='footprint2Mesh.py')
parser.add_argument("fileKey", help="Search string for collecting (.npz) files.",\
  nargs='?', default="npz")
parser.add_argument("-fo", "--fileout", type=str, default='FP',\
  help="Brief prefix for the footprint output file. (npz format)")
parser.add_argument("-ft", "--filetopo", type=str,\
  help="File containing the topography data. (npz format)")
parser.add_argument("-fm", "--filemean", type=str,\
  help="File containing the mean velocity data. (npz format)", default=None)
parser.add_argument("-N","--NxG", type=int,nargs=2,\
  help="Number of points [Nx, Ny] in the 2D Palm grid.")
parser.add_argument("-dx","--dxG", type=float,nargs=2,\
  help="Resolution [dx, dy] of the 2D Palm grid.")
#parser.add_argument("-fm", "--filemean", type=str,\
#  help="Name of the mean velocity .csv file.", default=None)
parser.add_argument("-v", "--vtk", help="Write VTK-files.",\
  action="store_true", default=False)
parser.add_argument("-i", "--ijk", help="Files contain ijk info.",\
  action="store_true", default=False) 
parser.add_argument("-cw", "--coefwm", type=float, default=1.,\
  help="Coefficient for scaling <w> for mean correction.")
parser.add_argument("-p", "--printOn", help="Print the extracted tile.",\
  action="store_true", default=False) 
parser.add_argument("-pp", "--printOnly", help="Only print the extracted tile. Don't save.",\
  action="store_true", default=False) 
args = parser.parse_args() 
#writeLog( parser, args )
#========================================================== #

# Rename ... that's all.
fileKey = args.fileKey
fileout  = args.fileout
filetopo = args.filetopo
filemean = args.filemean

NxG = args.NxG
dxG = args.dxG
cw_init  = args.coefwm

ijkOn     = args.ijk
vtkOn     = args.vtk
printOn   = args.printOn
printOnly = args.printOnly


# For writing the header once.
writeHeader = True

# Gather raw footprint data files: 
fileNos, fileList = filesFromList( fileKey+"*" )


if( filemean ):
  dat = np.load(filemean)
  wm  = dat['w']; xm = dat['x']; ym = dat['y']; zm = dat['z']
  dat = None
else:
  sys.exit(' Error. File for the mean values was not provided.')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = =   #
# xO := origin coords. # xt := target coords. # ut := target speed

for fn in fileNos:

  print(' Processing file: {}'.format(fileList[fn]))

  xO, yO, zO,\
    xt, yt, zt,\
      ut, vt, wt = readNumpyZFootprintRaw( fileList[fn] )

  # = = = Positive/Negative contributions = = = = = = = = = = #
  ''' 
  Now we make a dangerous assumption that the mean value of vertical
  velocity w_mean can be obtained from the particle data. This is strictly 
  not true and its effect must be carefully examined.
  '''
  
  wtm_1 = np.mean(wt)
  wtm_2, km, jm, im = meanFromExternal( wm, xm, ym, zm, xt, yt, zt )
  
  #wtm_3 = centralValue( wt, 1.e-6, 50 ) # Approximation for mean value.
  #wtm_4 = polarPercentileAvg( wt, 5, 95 )
  #wtm_5 = farFieldMean( wt, xO, 25. )   # Exclude the first 25%.
  '''
  - - - - - - - - - - - - -
  Location: x={0}, y={1}, z={2}\n
  Directly from particle data:  w_mean = {3:5.2f}
  Central (skimmed) value:      w_mean = {4:5.2f}
  Polar percentile (5,95):      w_mean = {5:5.2f}
  Far field mean (25% exluded): w_mean = {6:5.2f}
  Mean from external:           w_mean = {7:5.2f} 
  
  Selected:  w_mean = {8:6.3f}

  - - - - - - - - - - - - -
  '''
  
  # Store the mean values just for printing.
  xim = int(np.mean(xt)); yim = int(np.mean(yt)); zim = int(np.mean(zt))
  
  idx   = farFieldIds( xO, 10. )  # Consider the first 20% of the x-range.
  cw    = cw_init
  count = 0; ro    = None
  dr    = -1000.
  dc    = 0.005
  count_max = 100
  r_lim = 100
  infoStr = '''
  #- # - # - # - # - # - # - # - # - #
  Location: x={0}, y={1}, z={2}
  Directly from particle data:  wtm_1 = {3:5.2f}
  Mean from external:           wtm_2 = {4:5.2f}
  - - - - - - - - - - - - -
  '''.format(xim,yim,zim, wtm_1, wtm_2)
  print(infoStr)
 
  if( np.sum(km)>1 and np.sum(jm)>1 and np.sum(im)>1 ):
    wtm_ref =  wtm_2
  else:
    wtm_ref = wtm_1
  
  while( 1 ):

    if( count == 50 ):
      dc *= 2.

    wt_mean = cw * wtm_ref
    
    ipos  = ( (wt-wt_mean) > 0.)   # Boolean array for positive values. 
    ineg  = ~ipos                  # Boolean array for negative values. 
  
    # Form a loop that aims to equalize sum(ipos) and sum(ineg) when x < x_lim. 

    # Function evaluation
    r = abs( np.sum(ipos[idx]) - np.sum(ineg[idx]) )
    
    if( ro ):
      dr = (r - ro )
    
    if( dr > 0 ):
      dc *= -1.         # Bad direction, must change.
      cw += 1.5*dc      # Do not return to cwo
    else:
      cw += dc
  
    ro  = r
    count += 1
    
    if( (r < r_lim) or (count > count_max)):
      break

  
  infoItr = '''
  Iteration = {0}
  w_mean = {1:6.3f}\t vs. w_mean_orig = {2:6.3f}
  cw = {3}
  r  = {4}\t dr = {5}
  - - - - - - - - - - - - -
  '''.format(count, wt_mean, wtm_ref, cw, r, dr)
  print(infoItr)
  

  
  xt = None; yt = None; zt = None
  ut = None; vt = None; wt = None

  # = = = = 2d footprint evaluation. = = = = = = = = = = = =  # 

  # Determine grid coordinates for the footprint domain:
  xD, yD = coordsFootprintGrid( NxG, dxG, xO, yO )
  
  Nt = len(xO)
  print(' Total number of hits: Nt = {} '.format( Nt ))
  print(' Number of w>0 / w<0 hits: {} / {} '.format( len(xO[ipos]), len(xO[ineg] )))

  print(' Processing positive flux contributions ...')
  FMpos, XM, YM, ZMpos = fp2mshIJ( xO[ipos], yO[ipos], zO[ipos], xD, yD, dxG[0], dxG[1] )
  #print(' mean( FMpos ) = {}'.format(np.mean(FMpos)))

  print(' Processing negative flux contributions ...')
  FMneg, XM, YM, ZMneg = fp2mshIJ( xO[ineg], yO[ineg], zO[ineg], xD, yD, dxG[0], dxG[1] )
  #print(' mean( FMneg ) = {}'.format(np.mean(FMneg)))
  print(' ... done!')

  # Clear memory
  xO = None; yO = None; zO = None

  # Gather all the recorded z-coordinates.
  ZM = np.maximum( ZMpos , ZMneg )
  ZMpos = None; ZMneg = None

  print(' Gathering and Normalizing the footprint array ...')
  Cnorm = (Nt*dxG[0]*dxG[1])   # Coefficient for normalization.  
  FM = (FMpos - FMneg)/Cnorm;  FMpos = None;  FMneg = None
  print(' ... done!')
  
  
  # = = = = Output procedure = = = = = = = = = = = = = = =  #  
  fileId, varId = idAppendices(fileList[fn], ijkOn )
  writeNumpyZFootprint(fileout+fileId, FM, XM, YM, ZM, Cnorm )
  
  if( vtkOn ):
    if( writeHeader ):
      R, Rdims, ROrig, dPx = readNumpyZTile( filetopo )
      if( all(Rdims != np.shape(XM)) ):
        print(' Error! Mismatch Topo_dims={} vs. fp_dims={}'.format(Rdims,np.shape(XM)))
        sys.exit(1)
      
      f_vtk = vtkWriteHeaderAndGridStructured2d( XM, YM, R[::-1,:], fileout, 'Footprints')
      f_vtk = vtkWritePointDataHeader( f_vtk, FM, len(fileNos) )
      writeHeader = False; R=None 
   
    f_vtk = vtkWritePointDataStructured2D( f_vtk, FM   , XM, 'fp_'+varId )
    
    '''
    if( printOn or printOnly ):
      Cfp = addContourf( XM, YM, FM   , 'F(x,y)', fileout )
      Cfm = addContourf( XM, YM, FMneg, 'F_neg(x,y)', fileout+'_neg' )
      Cz = addContourf( XM, YM, ZM, 'Z(x,y)', ' Topography Height (m) ' )
      plt.show()
    '''
  FM = ZM = XM = YM = None

# Close the file at the end.
if( vtkOn ): f_vtk.close()
