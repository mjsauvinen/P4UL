#!/usr/bin/env python3
from netcdfTools import *
import sys
import argparse
import numpy as np
''' 
Description:

Author: Mikko Auvinen
        Finnish Meteorological Institute
'''
#==========================================================#

def idBounds( cmin, cmaxL, cminL, dc, nc ):
  ic1 = np.round((cminL-cmin)/dc , decimals=0).astype(int)
  ic2 = np.round((cmaxL-cmin)/dc , decimals=0).astype(int)+1
  ic2 = np.minimum( nc, ic2 )
  return ic1, ic2

#==========================================================#
#==========================================================#
parser = argparse.ArgumentParser(prog='concat3dNetCDF.py')
parser.add_argument("-f", "--filenames",type=str, nargs='+', default=None, \
  help="Name of the input NETCDF file.")
parser.add_argument("-fo", "--fileout",type=str, default="CONCAT.nc",\
  help="Name of the output NETCDF file.")
parser.add_argument("-sn", "--scalars",type=str, nargs='+', default=None,\
  help="(Optional) Scalars to be included.")
parser.add_argument("-so", "--scalarsOnly", action="store_true", default=False,\
  help="Only scalars, skip wind vector components.")
parser.add_argument("-nt", "--ntimeskip", type=int, default=0,\
  help="Skip <nt> number of time steps.")
args = parser.parse_args()
#==========================================================#
# Initial renaming operations and variable declarations

filenames  = args.filenames
fileout    = args.fileout
scalars    = args.scalars
ntskip     = args.ntimeskip
scalarsOnly= args.scalarsOnly

parameter = True;  variable  = False
dtol      = 1E-5  # distance tolerance

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

# Open output NetCDF file 
dso = netcdfOutputDataset( fileout )

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

if(scalarsOnly):
  vstr = scalars[0]
else:
  vstr = 'u'

if( filenames is not None ):
  Nf = len(filenames)
else:
  sys.exit(' No input NETCDF files given. Exiting ...')

dsL   = list()
xL    = list()
yL    = list() 
zL    = list()
timeL = list()

uL    = list()
vL    = list()
wL    = list()
sL    = list()

tv = None

xMaxL = np.zeros(Nf); yMaxL = np.zeros(Nf); zMaxL = np.zeros(Nf)
xMinL = np.zeros(Nf); yMinL = np.zeros(Nf); zMinL = np.zeros(Nf)

xMaxL[:] = yMaxL[:] = zMaxL[:] = -9999.
xMinL[:] = yMinL[:] = zMinL[:] =  9999.

dxo = dyo = dzo = None
nto  = None

for n,fn in enumerate(filenames):
  ds, varList, paramList = netcdfDataset(fn)
  dsL.append( ds )
  
  time, time_dims = read1DVariableFromDataset('time',vstr, ds, ntskip, 0, 1 ) # All values.
  print(' time dim = {} '.format(time_dims))
  if( nto is None ):
    nto = len(time)
  else:
    print(' Checking time array dimensions ...')
    if( nto != len(time) ): sys.exit(' Time dimensions do not match ')
  
  
  x, x_dims = read1DVariableFromDataset( 'x',vstr, ds, 0, 0, 1 )
  y, y_dims = read1DVariableFromDataset( 'y',vstr, ds, 0, 0, 1 )
  z, z_dims = read1DVariableFromDataset( 'z',vstr, ds, 0, 0, 1 )
  
  # A check should be entered 
  dx=np.mean(x[1:]-x[:-1]); dy=np.mean(y[1:]-y[:-1]); dz=np.mean(z[1:]-z[:-1])
  print('dx={:.3f}; dy={:.3f}; dz={:.3f}'.format(dx,dy,dz))  
  if( dxo is None ):
    dxo = dx; dyo = dy; dzo = dz
  else:
    print(' Checking that resolutions match ...')
    if( np.abs(dx-dxo)>dtol ): sys.exit(' dx do not match ')
    if( np.abs(dy-dyo)>dtol ): sys.exit(' dy do not match ')
    if( np.abs(dz-dzo)>dtol ): sys.exit(' dz do not match ')
  
  # Store coord arrays and record min and max values
  timeL.append(time)
  xL.append(x); xMaxL[n]=np.max(x); xMinL[n]=np.min(x)
  yL.append(y); yMaxL[n]=np.max(y); yMinL[n]=np.min(y)
  zL.append(z); zMaxL[n]=np.max(z); zMinL[n]=np.min(z)
  
  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # 
# Global coordinate bounds
xMax = np.max(xMaxL); xMin = np.min(xMinL)
yMax = np.max(yMaxL); yMin = np.min(yMinL)
zMax = np.max(zMaxL); zMin = np.min(zMinL)

# Create global coords
xc = np.arange( xMin, xMax+dx, dx )
yc = np.arange( yMin, yMax+dy, dy )
zc = np.arange( zMin, zMax+dz, dz )

nt = len(time)
nx = len(xc); ny = len(yc); nz = len(zc)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

tv = createNetcdfVariable( dso, time,'time', nt ,'s','f4',('time',), parameter )
xv = createNetcdfVariable( dso, xc   , 'x' , nx , 'm', 'f4', ('x',)   , parameter )
yv = createNetcdfVariable( dso, yc   , 'y' , ny , 'm', 'f4', ('y',)   , parameter )
zv = createNetcdfVariable( dso, zc   , 'z' , nz , 'm', 'f4', ('z',)   , parameter )

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
if( not scalarsOnly ):
  print('Concatenating u arrays ...')
  allocated = False 
  
  for n in range(Nf):
    u, u_dims = read3DVariableFromDataset( 'u', dsL[n], ntskip, 0, 0, 1 ) # All values.

    if( not allocated ):
      if( np.ma.is_masked(u) ): uc = np.ma.zeros( (nt, nz, ny, nx) )
      else:                     uc = np.zeros(    (nt, nz, ny, nx) )
      allocated = True
      print(' concat shape = {} '.format((nt, nz, ny, nx)))

    print(' u dims = {} '.format(u_dims))
    i1, i2 = idBounds(xMin, xMaxL[n], xMinL[n], dx, nx )
    print(' i1, i2 = {}, {} '.format(i1,i2))
  
    j1, j2 = idBounds(yMin, yMaxL[n], yMinL[n], dy, ny )
    print(' j1, j2 = {}, {} '.format(j1,j2))
  
    k1, k2 = idBounds(zMin, zMaxL[n], zMinL[n], dz, nz )
    print(' k1, k2 = {}, {} '.format(k1,k2))
  
    uc[:,k1:k2, j1:j2, i1:i2] = u[:,:,:,:]
    u = None
  
  uv = createNetcdfVariable(dso,uc,'u', nt, 'm/s', 'f4',('time','z','y','x',), variable)
  uc = None
  
  # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
  print('Concatenating v arrays ...')
  allocated = False
  
  for n in range(Nf):
    v, v_dims = read3DVariableFromDataset( 'v', dsL[n], ntskip, 0, 0, 1 ) # All values.
    
    if( not allocated ):
      if( np.ma.is_masked(v) ): vc = np.ma.zeros( (nt, nz, ny, nx) )
      else:                     vc = np.zeros(    (nt, nz, ny, nx) )
      allocated = True
    
    i1, i2 = idBounds( xMin, xMaxL[n], xMinL[n], dx, nx )
    j1, j2 = idBounds( yMin, yMaxL[n], yMinL[n], dy, ny )
    k1, k2 = idBounds( zMin, zMaxL[n], zMinL[n], dz, nz )
  
    vc[:,k1:k2, j1:j2, i1:i2] = v[:,:,:,:]
    v = None

  vv = createNetcdfVariable(dso,vc,'v', nt, 'm/s', 'f4',('time','z','y','x',), variable)
  vc = None
  # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
  print('Concatenating w arrays ...')
  allocated = False 
  
  for n in range(Nf):
    w, w_dims = read3DVariableFromDataset( 'w', dsL[n], ntskip, 0, 0, 1 ) # All values.
    
    if( not allocated ):
      if( np.ma.is_masked(w) ): wc = np.ma.zeros( (nt, nz, ny, nx) )
      else:                     wc = np.zeros(    (nt, nz, ny, nx) )
      allocated = True
    
    i1, i2 = idBounds( xMin, xMaxL[n], xMinL[n], dx, nx )
    j1, j2 = idBounds( yMin, yMaxL[n], yMinL[n], dy, ny )
    k1, k2 = idBounds( zMin, zMaxL[n], zMinL[n], dz, nz )
  
    wc[:,k1:k2, j1:j2, i1:i2] = w[:,:,:,:]
    w = None

  wv = createNetcdfVariable(dso,wc,'w', nt, 'm/s', 'f4',('time','z','y','x',), variable)
  wc = None
  # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

if( scalars ):
  sv = list()
  for sn in scalars:
    print('Concatenating {} arrays ...'.format(sn))
    allocated = False
    
    for n in range(Nf):
      s, s_dims = read3DVariableFromDataset( sn, dsL[n], ntskip, 0, 0, 1 ) # All values.
      
      if( not allocated ):
        if( np.ma.is_masked(s) ): sc = np.ma.zeros( (nt, nz, ny, nx) )
        else:                     sc = np.zeros(    (nt, nz, ny, nx) )
        allocated = True      
      
      i1, i2 = idBounds( xMin, xMaxL[n], xMinL[n], dx, nx )
      j1, j2 = idBounds( yMin, yMaxL[n], yMinL[n], dy, ny )
      k1, k2 = idBounds( zMin, zMaxL[n], zMinL[n], dz, nz )
  
      sc[:,k1:k2, j1:j2, i1:i2] = s[:,:,:,:]
      s = None

    sv.append(createNetcdfVariable(dso,sc,sn,nt,'-','f4',('time','z','y','x',),variable))
    sc = None

# - - - - Done , finalize the output - - - - - - - - - -
netcdfWriteAndClose( dso )
