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
  qtilde  = np.mean( q , axis=0 ) # mean at the moment
  qp      = q - qtilde
  qda     = np.mean( qtilde )     # double average
  qtilde -= qda
  
  return qda, qtilde, qp

#==========================================================#
parser = argparse.ArgumentParser(prog='vectorDecompNetCdf.py')
parser.add_argument('-f', '--filename', type=str, required=True,\
  help="Name of input NETCDF file.")
parser.add_argument("-fo", "--fileout",type=str, default="Vd.nc",\
  help="Name of output NETCDF file. Default=Vd.nc")
parser.add_argument("-vn", "--vnames",type=str, nargs=3, default=['u','v','w'],\
  help="Names of the vector components in (x,y,z)-order. Default = ['u','v','w'].")
parser.add_argument("-nt", "--ntimeskip", type=int, default=0,\
  help="Skip <nt> number of time steps. Default = 0.")
parser.add_argument("-c", "--coarse", type=int, default=1,\
  help="Coarsening level. Int > 1. Default = 1.")
args = parser.parse_args()
writeLog( parser, args )

#==========================================================#
# Initial renaming operations and variable declarations

filename   = args.filename
fileout    = args.fileout
vnames     = args.vnames
nt         = args.ntimeskip
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

udo = createNetcdfVariable(dso, uda   , 'uda'   , 1 , units, ft,('uda',) , parameter )
uto = createNetcdfVariable(dso, utilde, 'utilde', Nt, units, ft,('z','y','x',) , variable )
upo = createNetcdfVariable(dso, utilde, 'up'    , Nt, units, ft,('time','z','y','x',) , variable )

## v components  ##
vda, vtilde, vp = decomp3( v )
v = None

vdo = createNetcdfVariable(dso, vda   , 'vda'   , 1 , units, ft,('vda',) , parameter )
vto = createNetcdfVariable(dso, vtilde, 'vtilde', Nt, units, ft,('z','y','x',) , variable )
vpo = createNetcdfVariable(dso, vtilde, 'vp'    , Nt, units, ft,('time','z','y','x',) , variable )

## w components  ##
wda, wtilde, wp = decomp3( w )
w = None

wdo = createNetcdfVariable(dso, wda   , 'wda'   , 1 , units, ft,('wda',) , parameter )
wto = createNetcdfVariable(dso, wtilde, 'wtilde', Nt, units, ft,('z','y','x',) , variable )
wpo = createNetcdfVariable(dso, wtilde, 'wp'    , Nt, units, ft,('time','z','y','x',) , variable )

# - - - - Done , finalize the output - - - - - - - - - -
netcdfWriteAndClose( dso )
