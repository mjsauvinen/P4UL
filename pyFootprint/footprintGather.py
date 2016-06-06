#!/usr/bin/env python
from utilities import filesFromList, vtkWriteDataStructured2d, vtkWriteHeaderAndGridStructured2d
from utilities import vtkWritePointDataHeader, vtkWritePointDataStructured2D
from plotTools import addContourf, addToPlot
from footprintTools import *
from mapTools import readNumpyZTile, filterAndScale
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
# = # = # = # = # = # = # = # = # = # = # = # = # = # = #
# = # = # = # End Function definitions # = # = # = # = # = #

#========================================================== #
parser = argparse.ArgumentParser(prog='footprintGather.py')
parser.add_argument("fileKey", help="Search string for collecting (.npz) files.",\
  nargs='?', default="npz")
parser.add_argument("-fo", "--fileout", type=str, default='fp_gather',\
  help="Footprint output file. (npz format)")
parser.add_argument("-ft", "--filetopo", type=str,\
  help="File containing the topography data. (npz format)")
helpFlt = ''' Filter type and its associated number. Available filters:
 median, percentile, rank, gaussian, local, max. 
 Entering \"user num\" allows the user to specify <num> different filters consecutively.
 Example entry: median 5'''
parser.add_argument("-fl","--filter",type=str,nargs=2,default=[None,None], help=helpFlt)
parser.add_argument("-p", "--printOn", help="Print the contour of the footprint.",\
  action="store_true", default=False) 
parser.add_argument("-pp", "--printOnly", help="Only print the contour. Don't save.",\
  action="store_true", default=False) 
args = parser.parse_args() 
#writeLog( parser, args )
#========================================================== #

# Rename ... that's all.
fileKey  = args.fileKey
fileout  = args.fileout
filetopo = args.filetopo
flt      = args.filter
printOn   = args.printOn or args.printOnly
printOnly = args.printOnly

# Gather footprint data files: 
fileNos, fileList = filesFromList( "*"+fileKey+"*" )

# = = = = = = = = = = = = = = = = = = = = = = = = = = = =   #
# xO := origin coords. # xt := target coords. # ut := target speed

Ft = None; Ct = None; Zt = None

for fn in fileNos:

  print(' Processing file: {}'.format(fileList[fn]))

  Fi, X, Y, Z, Ci = readNumpyZFootprint( fileList[fn] )
  Fi *= Ci  # Return the footprint into unscaled state.
  
  if( Ft == None ):
    Ft = Fi.copy()
    Ct = Ci
    Zt = Z.copy(); Xt = X.copy(); Yt = Y.copy()
  else:
    Ft += Fi  # Accumulate the footprint data.
    Ct += Ci  # Accumulate the coefficient for normalization.
    Zt = np.maximum( Zt, Z )
    

# Compute the final footprint: 
Fi = X = Y = Z = Ci = None
Ft /= Ct

id50 = percentileFootprintIds( Ft, 50. )
id75 = percentileFootprintIds( Ft, 75. )
id90 = percentileFootprintIds( Ft, 90. )

# = = = = = = = = = = = = = = = = = = = = = = = = = = = =   #
# Apply filter if desired.
if( flt.count(None) == 0):
  Fft = np.zeros( np.shape(Ft) , float)
  Fft = filterAndScale(Fft, Ft, flt )
  Ft  = Fft.copy(); Fft = None

# = = = = = = = = = = = = = = = = = = = = = = = = = = = =   #
# Kormann et. Meixner Analytical footprint.
L = 2000.; z_m = (74.-20.); z_0 = 1.; sigma_v = 0.7; u=10.
x_off =  2.*228.; y_off = 2.*508.
F_km  = kormann_and_meixner_fpr(z_0, z_m, u, sigma_v, L, Xt, Yt, x_off, y_off)


# Compute the cross wind mean of the footprint.
Ftm   = writeCrossWindSum( Ft, Xt, fileout )
Ftm50 = writeCrossWindSum( Ft, Xt, fileout+'-50' , id50 )
Ftm75 = writeCrossWindSum( Ft, Xt, fileout+'-75' , id75 )
Ftm90 = writeCrossWindSum( Ft, Xt, fileout+'-90' , id90 )

Fm_km = writeCrossWindSum( F_km, Xt, fileout+'-ana' ) 



# = = = = = = = = = = = = = = = = = = = = = = = = = = = =   #
# Output to npz and vtk formats.
Ftmp = np.zeros( np.shape(Ft), float )

if( not printOnly ):
  IDict = {}
  IDict[50] = id50; IDict[75] = id75; IDict[90] = id90
  writeNumpyZFootprint(fileout, Ft, Xt, Yt, Zt, Ct, IDict )
  writeNumpyZFootprint('Fp_KormanMeixner', F_km, Xt, Yt, Zt, Ct )

  # Footprint to VTK-format together with the complete topography. 
  R, Rdims, ROrig, dPx = readNumpyZTile( filetopo )
  if( all(Rdims != np.shape(Xt)) ):
    sys.exit(' Error! Mismatch Topo_dims={} vs. fp_dims={}'.format(Rdims,np.shape(XM)))
  
  # Dnorm = np.percentile(Ft[Ft!=0], 80)
  Dnorm = np.max(Ft)
  Ft /= Dnorm
  
  f_vtk = vtkWriteHeaderAndGridStructured2d( Xt, Yt, R[::-1,:], fileout, 'Footprint'); R=None 
  f_vtk = vtkWritePointDataHeader( f_vtk, Ft, 5 )
  f_vtk = vtkWritePointDataStructured2D( f_vtk, Ft , Xt, 'fp' )
  
  Ftmp += Ft*id50
  f_vtk = vtkWritePointDataStructured2D( f_vtk, Ftmp , Xt, 'fp50' )
  
  Ftmp = 0.; Ftmp += Ft*id75
  f_vtk = vtkWritePointDataStructured2D( f_vtk, Ftmp , Xt, 'fp75' )
  
  Ftmp = 0.; Ftmp += Ft*id90
  f_vtk = vtkWritePointDataStructured2D( f_vtk, Ftmp , Xt, 'fp90' )
  
  # Dnorm = np.percentile(F_km[F_km!=0], 80)
  Dnorm = np.max(F_km)
  F_km /= Dnorm
  
  f_vtk = vtkWritePointDataStructured2D( f_vtk, F_km, Xt, 'fp_km' )

  # Close the file at the end.
  f_vtk.close(); Ftmp = None

if( printOn ):
  Cfp = addContourf( Xt, Yt, Ft   , 'F(x,y)', fileout )
  Cfa = addContourf( Xt, Yt, F_km   , 'F_km(x,y), Ana', fileout+'-ana' )
  
  pfig = plt.figure(num=3, figsize=(12.,9.))
  varLabel = 'fp_y(x) = sum_y fp(x,y)'
  axLabels = ['Cross Wind Integrated Footprint', 'x', 'sum_y fp(x,y) ']
  pfig = addToPlot(pfig, Xt[0,:], Ftm, varLabel, axLabels, False )


Ft = Zt = Xt = Yt = None
F_km = None

plt.show()
