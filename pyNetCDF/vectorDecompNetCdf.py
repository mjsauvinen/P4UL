#!/usr/bin/env python3
from netcdfTools import *
import sys
import argparse
import numpy as np
from utilities import filesFromList, writeLog
''' 
Description: Vector decomposition U = <ubar> + utilde(x) + up(x,t)


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''
#==========================================================#

def decomp3( q ):
  idn = (q[-1,:,:,:]==0.) # take the zeros from last time step
  nt = np.shape(q)[0]
  for i in range(nt):
    q[i,idn] = np.nan 
  idn = None
  
  qtilde  = np.nanmean( q , axis=0 ) # mean at the moment
  qp      = q - qtilde
  qda     = np.nanmean( qtilde )     # double average
  qtilde -= qda
  
  return np.array([qda]), qtilde, qp

#==========================================================#

def normReynodsStressTensor(up,vp,wp):
  
  print(' Computing Reynolds stresses and their norm ...')
  R11 = np.nanmean(up*up, axis=0)
  R22 = np.nanmean(vp*vp, axis=0)
  R33 = np.nanmean(wp*wp, axis=0)
  
  tr = (1./3.)*( R11 + R22 + R33 )
  norm  = R11**2; R11 -= tr; devnorm  = R11**2; R11 = None
  norm += R22**2; R22 -= tr; devnorm += R22**2; R22 = None
  norm += R33**2; R33 -= tr; devnorm += R33**2; R33 = None
  
  R12p = 2.*np.nanmean(up*vp, axis=0)**2
  norm += R12p; devnorm += R12p; R12p = None
  
  R13p = 2.*np.nanmean(up*wp, axis=0)**2
  norm += R13p; devnorm += R13p; R13p = None
  
  R23p = 2.*np.nanmean(vp*wp, axis=0)**2
  norm += R23p; devnorm += R23p; R23p = None
  
  norm    **=(0.5)
  devnorm **=(0.5)
  print(' ... done!')
  
  return norm, devnorm, (3./2.)*tr

#==========================================================#

def cleanValues( D, Dstr, rV=0. ):
  # Set all Nans and bad values to zeros
  if( len(np.shape(D)) == 4 ):
    
    Nt = np.shape(D)[0]
    idn = np.isnan( D[-1,:,:,:] )
    nn = np.count_nonzero( idn )
    if( nn > 0 ):
      print(' {}: Number of nan values = {}'.format(Dstr,nn))
      for i in range(Nt): 
        D[i,idn] = rV
      
    idn = ( np.abs(D[-1,:,:,:]) > 1.E6 )
    nn = np.count_nonzero( idn )
    if( nn > 0 ):
      print(' {}: Number of ill defined values = {}'.format(Dstr,nn))
      for i in range(Nt):
        D[i,idn] = rV
    idn = None
    
  else:
    
    idn = np.isnan( D )
    nn = np.count_nonzero( idn )
    if( nn > 0 ):
      print(' {}: Number of nan values = {}'.format(Dstr,nn))
      D[idn] = rV
    
    idn = ( np.abs(D) > 1.E6 )
    nn = np.count_nonzero( idn )
    if( nn > 0 ):
      print(' {}: Number of ill defined values = {}'.format(Dstr,nn))
      D[idn] = rV
    
    idn = None
  
  return D

#==========================================================#
parser = argparse.ArgumentParser(prog='vectorDecompNetCdf.py')
parser.add_argument('-f', '--filename', type=str, required=True,\
  help="Name of input NETCDF file.")
parser.add_argument("-fo", "--fileout",type=str, default="Vd.nc",\
  help="Name of output NETCDF file. Default=Vd.nc")
parser.add_argument("-vn", "--vnames",type=str, nargs=3, default=['u','v','w'],\
  help="Names of the vector components in (x,y,z)-order. Default = ['u','v','w'].")
parser.add_argument('-rs',"--reynoldsStresses", action="store_true", default=False,\
  help="Compute and write magnitude of Reynolds stresses.")
parser.add_argument("-nt", "--ntimeskip", type=int, default=0,\
  help="Skip <nt> number of time steps. Default = 0.")
parser.add_argument('-m',"--mags", action="store_true", default=False,\
  help="Compute and write magnitudes of each decomposition component.")
parser.add_argument("-c", "--coarse", type=int, default=1,\
  help="Coarsening level. Int > 1. Default = 1.")
args = parser.parse_args()
writeLog( parser, args )

#==========================================================#
# Initial renaming operations and variable declarations. 

filename   = args.filename
fileout    = args.fileout
vnames     = args.vnames
nt         = args.ntimeskip
magsOn     = args.mags
rsOn       = args.reynoldsStresses
cl         = abs(int(args.coarse))

'''
Establish two boolean variables which indicate whether the created variable is an
independent or dependent variable in function createNetcdfVariable().
'''
parameter = True;  variable = False

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# Read in data.
dataDict = read3dDataFromNetCDF( filename , vnames , cl )
u = dataDict.pop(vnames[0])
v = dataDict.pop(vnames[1])
w = dataDict.pop(vnames[2])

# Coords and time:
x  = dataDict.pop('x'); y = dataDict.pop('y'); z = dataDict.pop('z')
time = dataDict.pop('time')
Nt = len(time)
dataDict = None

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# Create a NETCDF output dataset (dso) for writing out the data.
dso = netcdfOutputDataset( fileout )
ft = 'f4'
units = 'm s^(-1)'


# Create the output independent variables right away and empty memory.
tv = createNetcdfVariable( dso, time,'time', Nt,'s',ft,('time',), parameter )
time = None  

xv = createNetcdfVariable( dso, x , 'x' , len(x) , 'm', ft, ('x',) , parameter )
x = None

yv = createNetcdfVariable( dso, y , 'y' , len(y) , 'm', ft, ('y',) , parameter )
y = None

zv = createNetcdfVariable( dso, z , 'z' , len(z) , 'm', ft, ('z',) , parameter )
z = None

## u components  ##
uda, utilde, up = decomp3( u )
u = None

utilde = cleanValues(utilde, 'utilde')
udo = createNetcdfVariable(dso, uda   , 'uda'   , 1 , units, ft,('uda',) , parameter )
uto = createNetcdfVariable(dso, utilde, 'utilde', 1 , units, ft,('z','y','x',) , variable )


if( magsOn ):
  Udamag    = uda**2
  Upmag     = up**2     
  Utildemag = utilde**2

utilde = None

## v components  ##
vda, vtilde, vp = decomp3( v )
v = None

vtilde = cleanValues(vtilde, 'vtilde')
vdo = createNetcdfVariable(dso, vda   , 'vda'   , 1 , units, ft,('vda',) , parameter )
vto = createNetcdfVariable(dso, vtilde, 'vtilde', 1 , units, ft,('z','y','x',) , variable )


if( magsOn ):
  Udamag    += vda**2
  Upmag     += vp**2     
  Utildemag += vtilde**2

vtilde = None


## w components  ##
wda, wtilde, wp = decomp3( w )
w = None

wtilde = cleanValues(wtilde, 'wtilde')
wdo = createNetcdfVariable(dso, wda   , 'wda'   , 1 , units, ft,('wda',) , parameter )
wto = createNetcdfVariable(dso, wtilde, 'wtilde', 1 , units, ft,('z','y','x',) , variable )


if( magsOn ):
  Udamag    += wda**2
  Upmag     += wp**2
  Utildemag += wtilde**2 

  Udamag **=(0.5); Upmag **=(0.5); Utildemag **=(0.5)
  Upmag = cleanValues( Upmag, '|Up|' ); Utildemag = cleanValues( Utildemag, '|Utilde|')

  Udamag = np.concatenate( (Udamag, np.array([np.nanmean(Utildemag), np.nanmean(Upmag) ])) )

  Udo = createNetcdfVariable(dso, Udamag   , 'Uda'   , 3 , units, ft,('Uda',) , parameter )
  Uto = createNetcdfVariable(dso, Utildemag, 'Utilde', 1 , units, ft,('z','y','x',) , variable )
  Upo = createNetcdfVariable(dso, Upmag    , 'Up'    , Nt, units, ft,('time','z','y','x',) , variable )
  Udamag = Upmag = Utildemag = None

if( rsOn ):
  normRS, normDevRS, TKE = normReynodsStressTensor( up, vp, wp )
  
  normRS    = cleanValues( normRS   , 'normRS' )
  normDevRS = cleanValues( normDevRS, 'normDevRS' )
  TKE       = cleanValues( TKE      , 'TKE' )
  
  mNRS  = np.nanmean( normRS )
  mNdRS = np.nanmean( normDevRS )
  mTKE  = np.nanmean( TKE )
  
  mNRS = np.array([mNRS, mNdRS, mTKE])
  
  Nro = createNetcdfVariable(dso, mNRS     , 'mNRS'  , 3 ,'m^2 s^-2', ft,('mNRS',) , parameter )
  RSo = createNetcdfVariable(dso, normRS   , 'nRS'   , 1,'m^2 s^-2', ft,('z','y','x',) , variable, zlib=False, fill_value=0. )
  dRSo= createNetcdfVariable(dso, normDevRS, 'nDevRS', 1,'m^2 s^-2', ft,('z','y','x',) , variable )
  TKo = createNetcdfVariable(dso, TKE      , 'TKE'   , 1,'m^2 s^-2', ft,('z','y','x',) , variable )


up = cleanValues( up, 'up' ); vp = cleanValues( vp, 'vp' ); wp = cleanValues( wp, 'wp' )
upo = createNetcdfVariable(dso, up    , 'up'    , Nt, units, ft,('time','z','y','x',) , variable ); up = None
vpo = createNetcdfVariable(dso, vp    , 'vp'    , Nt, units, ft,('time','z','y','x',) , variable ); vp = None
wpo = createNetcdfVariable(dso, wp    , 'wp'    , Nt, units, ft,('time','z','y','x',) , variable ); wp = None


# - - - - Done , finalize the output - - - - - - - - - -
netcdfWriteAndClose( dso )
