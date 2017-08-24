#!/usr/bin/env python
from utilities import filesFromList, vtkWriteDataStructured2d, vtkWriteHeaderAndGridStructured2d
from utilities import vtkWritePointDataHeader, vtkWritePointDataStructured2D
from utilities import writeLog
from plotTools import addContourf, addToPlot
from footprintTools import *
from mapTools import readNumpyZTile, filterAndScale
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
''' 
Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''
# = # = # = # Function definitions # = # = # = # = # = # = # 
# = # = # = # = # = # = # = # = # = # = # = # = # = # = #
# = # = # = # End Function definitions # = # = # = # = # = #

#========================================================== #
parser = argparse.ArgumentParser(prog='footprintGather.py')
parser.add_argument("fileKey", help="Search string for collecting (.npz) files.",\
  nargs='?', default="npz")
parser.add_argument("-a", "--allfiles", help="Select all files automatically.",\
  action="store_true", default=False) 
parser.add_argument("-fo", "--fileout", type=str, default='fp_gather',\
  help="Footprint output file. (npz format)")
parser.add_argument("-ft", "--filetopo", type=str,\
  help="File containing the topography data. (npz format)", default='')
helpFlt = ''' Filter type and its associated number. Available filters:
 median, percentile, rank, gaussian, local, max. 
 Entering \"user num\" allows the user to specify <num> different filters consecutively.
 Example entry: median 5'''
parser.add_argument("-fl","--filter",type=str,nargs=2,default=[None,None], help=helpFlt)
parser.add_argument("-n1", "--norm2one", help="Normalize by making global sum = 1.",\
  action="store_true", default=False)
parser.add_argument("-v","--vtk", help="Write the results in VTK format with topography.",\
  action="store_true", default=False) 
parser.add_argument("-p", "--printOn", help="Print the contour of the footprint.",\
  action="store_true", default=False) 
parser.add_argument("-pp", "--printOnly", help="Only print the contour. Don't save.",\
  action="store_true", default=False) 
args = parser.parse_args() 
writeLog( parser, args, args.printOnly )
#========================================================== #

# Rename ... that's all.
fileKey   = args.fileKey
fileout   = args.fileout
filetopo  = args.filetopo
flt       = args.filter
allFiles  = args.allfiles
norm2one  = args.norm2one
printOn   = args.printOn or args.printOnly
printOnly = args.printOnly
vtkOn     = args.vtk

if( vtkOn and (filetopo == '')):
  sys.exit(' Error! VTK results require -ft/--filetopo. Exiting ...')


# Gather footprint data files: 
fileNos, fileList = filesFromList( "*"+fileKey+"*", allFiles )

# = = = = = = = = = = = = = = = = = = = = = = = = = = = =   #
# xO := origin coords. # xt := target coords. # ut := target speed

Ft = None; Ct = None; Zt = None

for fn in fileNos:

  print(' Processing file: {}'.format(fileList[fn]))

  Fi, X, Y, Z, Ci = readNumpyZFootprint( fileList[fn] )
  Fi *= Ci  # Return the footprint into unscaled state.
  
  if( Ft is None ):
    Ft = Fi.copy()
    Ct = Ci
    Zt = Z.copy(); Xt = X.copy(); Yt = Y.copy()
  else:
    Ft += Fi  # Accumulate the footprint data.
    Ct += Ci  # Accumulate the coefficient for normalization.
    Zt = np.maximum( Zt, Z )

# = = = = = = = = = = = = = = = = = = = = = = = = = = = =   #
# Resolution:
dPx = np.array([ (Xt[0,1]-Xt[0,0]) , (Yt[1,0]-Yt[0,0]) ])

# = = = = = = = = = = = = = = = = = = = = = = = = = = = =   #
# Compute the final footprint: 
Fi = X = Y = Z = Ci = None
Ft /= Ct

# = = = = = = = = = = = = = = = = = = = = = = = = = = = =   #
# Apply filter if desired.
if( flt.count(None) == 0):
  Fft = np.zeros( np.shape(Ft) , float)
  Fft = filterAndScale(Fft, Ft, flt )
  Ft  = Fft.copy(); Fft = None

# = = = = = = = = = = = = = = = = = = = = = = = = = = = =   #
# Kormann et. Meixner Analytical footprint.
# Tower: u = 4.86,       sigma_v = 0.75, 
# LES  : u = 4.5 - 5.1,  sigma_v = 0.72-0.74
# Upwind LES mean over 1200 m: 
#        u = 6.1,        sigma_v = 0.95
L =10000.; z_m = (60.-14.9); z_0 = 1.4; sigma_v = 0.75; u=4.86
x_off = 2.*228.; y_off = 2.*508.
F_km  = kormann_and_meixner_fpr(z_0, z_m, u, sigma_v, L, Xt, Yt, x_off, y_off)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = =   #
# Make the analytical and LES footprints comparable
# Force the global sum = 1.
if( norm2one ):
  print(' Normalizing the footprints such that SUM(Fp) = 1 ...')
  Cn = 1./np.sum( Ft  * np.prod(dPx));  Ft  *= Cn
  Ca = 1./np.sum( F_km* np.prod(dPx));  F_km*= Ca
  print('... done! C_les = {} and C_ana = {}'.format(Cn, Ca))

# = = = = = = = = = = = = = = = = = = = = = = = = = = = =   #
# Extract indecies for partial (%) footprints

id50 = percentileFootprintIds( Ft, 50. )
id75 = percentileFootprintIds( Ft, 75. )
id90 = percentileFootprintIds( Ft, 90. )

id90_km = percentileFootprintIds( F_km, 90. ) # 90% 
id75_km = percentileFootprintIds( F_km, 75. ) # 75% 

# = = = = = = = = = = = = = = = = = = = = = = = = = = = =   #
# Output to npz and vtk formats.

if( not printOnly ):
  
  # Compute the cross wind mean of the footprint.
  Ftm   = writeCrossWindSum( Ft, Xt, fileout )
  #Ftm50 = writeCrossWindSum( Ft, Xt, fileout+'-50' , id50 )
  #Ftm75 = writeCrossWindSum( Ft, Xt, fileout+'-75' , id75 )
  #Ftm90 = writeCrossWindSum( Ft, Xt, fileout+'-90' , id90 )

  Fm_km   = writeCrossWindSum( F_km, Xt, fileout+'_km' )
  #Fm90_km = writeCrossWindSum( F_km, Xt, fileout+'-90_km', id90_km )
  
  
  # Write the footprint in npz format.
  IDict = {}
  IDict[50] = id50[::-1,:]; IDict[75] = id75[::-1,:]; IDict[90] = id90[::-1,:]
  writeNumpyZFootprint(fileout, Ft[::-1,:], Xt, Yt, Zt, Ct, IDict )
  
  # Write also the Kormann-Meixner footprint
  IDict = {}
  IDict[75] = id75_km[::-1,:]; IDict[90] = id90_km[::-1,:]
  writeNumpyZFootprint(fileout+'_KormannMeixner', F_km[::-1,:], Xt, Yt, Zt, Ct, IDict )
  
  
  if( vtkOn ):
    
    # Footprint to VTK-format together with the complete topography. 
    Ftmp = np.zeros( np.shape(Ft), float )
    R, Rdims, ROrig, dPx = readNumpyZTile( filetopo )
    if( all(Rdims != np.shape(Xt)) ):
      sys.exit(' Error! Mismatch Topo_dims={} vs. Fp_dims={}'.format(Rdims,np.shape(Xt)))
  

    f_vtk = vtkWriteHeaderAndGridStructured2d( Xt, Yt, R[::-1,:], fileout, 'Footprint'); R=None 
    f_vtk = vtkWritePointDataHeader( f_vtk, Ft, 5 )
  
    # ======= Write 100% Ft ================
    f_vtk = vtkWritePointDataStructured2D( f_vtk, Ft , Xt, 'fp' )
  
    # ======= Write 75% Ft ================
    Ftmp[:,:] = 0.; Ftmp += Ft*id75
    f_vtk = vtkWritePointDataStructured2D( f_vtk, Ftmp , Xt, 'fp75' )
  
    # ======= Write 90% Ft ================
    Ftmp[:,:] = 0.; Ftmp += Ft*id90
    f_vtk = vtkWritePointDataStructured2D( f_vtk, Ftmp , Xt, 'fp90' )
  
    # ======= Write 100% F_km ================
    f_vtk = vtkWritePointDataStructured2D( f_vtk, F_km, Xt, 'fp_km' )
  
    # ======= Write 00% F_km ================
    Ftmp[:,:] = 0.; Ftmp += F_km*id90_km
    f_vtk = vtkWritePointDataStructured2D( f_vtk, Ftmp , Xt, 'fp90_km' )

    # Close the file at the end.
    f_vtk.close(); Ftmp = None

if( printOn ):
  Cfp = addContourf( Xt, Yt, Ft  , 'F(x,y)'        , fileout )
  Cfa = addContourf( Xt, Yt, F_km, 'F_km(x,y), Ana', fileout+'_km' )
  
  Fym = writeCrossWindSum( Ft, Xt, None, None )
  pfig = plt.figure(num=3, figsize=(12.,9.))
  varLabel = '$fp_y(x) = \sum_y fp(x,y)$'
  axLabels = ['Cross Wind Integrated Footprint', 'x', 'sum_y fp(x,y) ']
  pfig = addToPlot(pfig, Xt[0,:], Fym, varLabel, axLabels, False )


Ft = Zt = Xt = Yt = None
F_km = None

plt.show()
