#!/usr/bin/env python
from utilities import filesFromList, vtkWriteDataStructured2d, vtkWriteHeaderAndGridStructured2d
from utilities import vtkWritePointDataHeader, vtkWritePointDataStructured2D
from utilities import writeLog
from plotTools import addContourf, extractFromCSV
from footprintTools import *
from mapTools import readNumpyZTile, farFieldIds, farFieldMean
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

def externalGridIds( ce, cl ):
  cbmin = max( np.min(cl) , np.min(ce) ) 
  cbmax = min( np.max(cl) , np.max(ce) )
  #if( cbmin > cbmax ): cbmin = cbmax + 1e-5
  
  jg = (ce>=cbmin); jl = (ce<=cbmax)
  
  # Check if 'greater than' and 'smaller than' are complements 
  jcompl = all(~jl==jg) and not all(jl==True) and not all(jl==False)
  
  # Check resolutions
  dce = ce[1]-ce[0]  # delta of external grid 
  dcl = cl[1]-cl[0]  # delta of local box
  finer_dcl = dcl < dce
  
  if( finer_dcl and jcompl ):
    idN = np.arange(len(jl))
    id2 = np.array([ idN[jl][-1], idN[jg][0] ])
    je  = np.zeros( len(jl) , bool )
    je[id2] = True
  else:
    je = jg * jl

  return je

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #

def meanFromExternal( dat, xe, ye, ze, xl, yl, zl):
  # Define the vertical mean flow from the external data.
  # <#>e: external, <#>l: local

  ke = externalGridIds( ze, zl )
  je = externalGridIds( ye, yl )
  ie = externalGridIds( xe, xl )
  
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
parser.add_argument("-a", "--allfiles", help="Select all files automatically.",\
  action="store_true", default=False) 
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
parser.add_argument("-b","--hybrid", help="Hybrid approach with far field correction.",\
  action="store_true", default=False)
parser.add_argument("--vtk", help="Write VTK-files.",\
  action="store_true", default=False)
parser.add_argument("-i", "--ijk", help="Files contain ijk info.",\
  action="store_true", default=False) 
parser.add_argument("-cw", "--coefwm", type=float, default=1.,\
  help="Coefficient for scaling <w> for mean correction.")
help_px ='''Percentage of first x-koords where the footprint is set to zero (fp=0).
If not specified, (Default=None) the farfield correction is not performed.''' 
parser.add_argument("-px","--pxzero", type=float, default=None, help=help_px)
parser.add_argument("-p", "--printOn", help="Print the extracted tile.",\
  action="store_true", default=False) 
parser.add_argument("-pp", "--printOnly", help="Only print the extracted tile. Don't save.",\
  action="store_true", default=False) 
parser.add_argument("-v", "--verbose", help="Print all information on screen.",\
  action="store_true", default=False) 
args = parser.parse_args() 
writeLog( parser, args )
#========================================================== #

# Rename ... that's all.
fileKey = args.fileKey
fileout  = args.fileout
filetopo = args.filetopo
filemean = args.filemean

NxG = args.NxG
dxG = args.dxG
cw_init  = args.coefwm
pxz      = args.pxzero

allFiles  = args.allfiles
hybridOn  = args.hybrid
ijkOn     = args.ijk
vtkOn     = args.vtk
printOn   = args.printOn
printOnly = args.printOnly
verbose   = args.verbose


# For writing the header once.
writeHeader = True

# Gather raw footprint data files: 
fileNos, fileList = filesFromList( fileKey+"*", allFiles )


if( filemean ):
  dat = np.load(filemean)
  wm  = dat['w']; xm = dat['x']; ym = dat['y']; zm = dat['z']
  dat = None
else:
  sys.exit(' Error. File for the mean values was not provided.')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = =   #
# xO := origin coords. # xt := target coords. # ut := target speed

for fn in fileNos:

  if( verbose ): print(' Processing file: {}'.format(fileList[fn]))

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
  
  # Choose the reference vertical velocity.
  if( np.sum(km)>0 and np.sum(jm)>0 and np.sum(im)>0 ):
    wtm_ref     =  wtm_2
    meanFromExt = True       # mean from external 
  else:
    wtm_ref     = wtm_1
    meanFromExt = False
  
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
  
  # Boolean to determine whether far field correction is performed.
  farFieldCorrOn = (pxz is not None)
  if( hybridOn ):
    farFieldCorrOn = farFieldCorrOn and (not meanFromExt)
  
  
  if( farFieldCorrOn ):
    print(' Far Field Correction! ')
  
    idx   = farFieldIds( xO, pxz )  # Consider the first 15% (=default) of the x-range.
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
    if( verbose ): print(infoStr)
  
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

    # end while
  
    infoItr = '''
    Iteration = {0}
    w_mean = {1:6.3f}\t vs. w_mean_orig = {2:6.3f}
    cw = {3}
    r  = {4}\t dr = {5}
    - - - - - - - - - - - - -
    '''.format(count, wt_mean, wtm_ref, cw, r, dr)
    if( verbose ): print(infoItr)
  
  elif( meanFromExt ): # no farfield correction
    ipos  = ( (wt-wtm_ref) > 0.)   # Boolean array for positive values.
    ineg  = ~ipos                  # Boolean array for negative values.
  else:
    continue
  

  # Clear memory
  xt = None; yt = None; zt = None
  ut = None; vt = None; wt = None

  # = = = = 2d footprint evaluation. = = = = = = = = = = = =  # 

  # Determine grid coordinates for the footprint domain:
  xD, yD = coordsFootprintGrid( NxG, dxG, xO, yO, verbose )
  
  Nt = len(xO); Ntp = len(xO[ipos]); Ntn = len(xO[ineg])
  print(' > Nt, Nt(pos), Nt(neg), <xO>+, <xO>- xt, yt, zt = {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'\
    .format( Nt, Ntp, Ntn, np.mean(xO[ipos]), np.mean(xO[ineg]), xim, yim, zim ))
  
  
  #print(' Number of w>0 / w<0 hits: {} / {} '.format( Ntp, Ntn))

  if( verbose ): print(' Processing positive flux contributions ...')
  FMpos, XM, YM, ZMpos = fp2mshIJ( xO[ipos], yO[ipos], zO[ipos], xD, yD, dxG[0], dxG[1] )
  #print(' mean( FMpos ) = {}'.format(np.mean(FMpos)))

  if( verbose ): print(' Processing negative flux contributions ...')
  FMneg, XM, YM, ZMneg = fp2mshIJ( xO[ineg], yO[ineg], zO[ineg], xD, yD, dxG[0], dxG[1] )
  #print(' mean( FMneg ) = {}'.format(np.mean(FMneg)))
  if( verbose ): print(' ... done!')

  # Clear memory
  xO = None; yO = None; zO = None

  # Gather all the recorded z-coordinates.
  ZM = np.maximum( ZMpos , ZMneg )
  ZMpos = None; ZMneg = None

  if( verbose ): print(' Gathering and Normalizing the footprint array ...')
  Cnorm = (Nt*dxG[0]*dxG[1])   # Coefficient for normalization.  
  FM = (FMpos - FMneg)/Cnorm;  FMpos = None;  FMneg = None
  if( verbose ): print(' ... done!')
  
  
  # = = = = Output procedure = = = = = = = = = = = = = = =  #  
  fileId, varId = idAppendices(fileList[fn], ijkOn )
  writeNumpyZFootprint(fileout+fileId, FM, XM, YM, ZM, Cnorm )
  
  if( vtkOn ):
    if( writeHeader ):
      Rdict = readNumpyZTile( filetopo )
      R = Rdict['R']
      Rdims = np.array(np.shape(R))
      ROrig = Rdict['LocalOrig']
      dPx = Rdict['dPx']
      Rdict = None
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
