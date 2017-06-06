#!/usr/bin/env python
from netcdfTools import *
from mapTools import *
from utilities import writeLog
import sys
import argparse
import numpy as np
''' 
Description:


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''

#============= functions ==================================#

def domainBoundsAndResolution( xf, yf ):
  xb = np.array([ np.min(xf), np.max(xf) ])
  yb = np.array([ np.min(yf), np.max(yf) ])
  dx = (xb[1]-xb[0])/float(len(xf)-1)
  dy = (yb[1]-yb[0])/float(len(yf)-1)
  
  return xb, yb, dx, dy

#============ Main     ====================================#
parser = argparse.ArgumentParser(prog='syncMaskWithNetCDF.py')
parser.add_argument("-fn", "--fileNetCDF",type=str, help="Name of input NETCDF file.")
parser.add_argument("-fm", "--fileMask",type=str, help="Name of input 2D Mask file.")
parser.add_argument("-d", "--decomp", action="store_true", default=False, \
    help="Decomposed into mean (V_m) and fluctuating (V^prime) components.") 
parser.add_argument("-dd", "--decompOnly", help="Output V_m and V^prime components only.",\
  action="store_true", default=False)
parser.add_argument("-c", "--coarse", help="Coarsening level for the NETCDF data. Int > 1.",\
  type=int, default=1) 
args = parser.parse_args()
writeLog( parser, args )
#==========================================================#
# Initial renaming operations and variable declarations

fnc   = args.fileNetCDF
fmsk  = args.fileMask
fout  = fnc.strip('.nc')+'-Msk.nc'
cl    = abs(int(args.coarse))


# Boolean switch for the decomposition option.
decompOn = args.decomp or args.decompOnly

'''
Establish two boolean variables which indicate whether the created variable is an
independent or dependent variable in function createNetcdfVariable().
'''
parameter = True;  variable  = False

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
''' 
Create a NETCDF input dataset (ds), and its associated lists of dependent (varList)
and independent (dimList) variables. 
'''
ds, varList, paramList = netcdfDataset(fnc)

# Create a NETCDF output dataset (dso) for writing out the data.
dso = netcdfOutputDataset( fout )

'''
Read cell center coordinates and time.
Create the output independent variables right away and empty memory.
'''

time, time_dims = read1DVariableFromDataset('time', ds, paramList, 0, 0, 1 ) # All values.
tv = createNetcdfVariable( dso, time,'time', len(time),'s','f4',('time',), parameter )
time = None  

x, x_dims = read1DVariableFromDataset( 'x',ds, paramList, 0, 0, cl ) # All values.
print(' x_dims = {} '.format(x_dims))
x[np.isnan(x)] = 0.  # Special treatment.
xv = createNetcdfVariable( dso, x   , 'x'   , len(x)   , 'm', 'f4', ('x',)   , parameter )

y, y_dims = read1DVariableFromDataset( 'y',ds, paramList, 0, 0, cl )
print(' y_dims = {} '.format(y_dims))
y[np.isnan(y)] = 0.  # Special treatment.
yv = createNetcdfVariable( dso, y   , 'y'   , len(y)   , 'm', 'f4', ('y',)   , parameter )

# Determine the NETCDF domain bounds and resolution.
xb, yb, dx, dy = domainBoundsAndResolution( x, y )
x = None; y = None # Clear memory ASAP.

z, z_dims = read1DVariableFromDataset( 'z',ds, paramList, 0, 0, cl )
print(' z_dims = {} '.format(z_dims))
zv = createNetcdfVariable( dso, z   , 'z'   , len(z)   , 'm', 'f4', ('z',)   , parameter )
z = None

# - - - - First, read u-component - - - - - - - - - -
u, u_dims = read3DVariableFromDataset( 'u', ds, varList, 0, 0, cl ) # All values.
print(' u_dims = {} '.format(u_dims))
yx_dims = np.array(u_dims[2:])
z_dim   = u_dims[1]; t_dim = u_dims[0]


'''
At this point the mask raster data can be treated because 
it needs one scalar NETCDF variable to determine the required
index bounds and coarsening level.
'''


# Read the mask raster info.
Rdict = readNumpyZTile(fmsk)
R = Rdict['R']
R_dims = np.array(np.shape(R))
ROrig = Rdict['LocalOrig']
dPx = Rdict['dPx']
Rdict = None
dr = entry2Int( dPx )  # Resolution as a single number
clr = int( dx/dr )     # Raster to NETCDF coarsening factor
print(' Orig mask dims = {} '.format(R_dims))

# We need y_max value for the Raster data to determine the reversed j-indecies. 
ybr_max = R_dims[0]*dr
print(' ybr_max = {}'.format(ybr_max))

# Determine the index range for the raster data to match the NETCDF (sub)domain.
# NOTE: dy is subtracted to make first index 0-based.
irx = np.array([ int(xb[0]-dy)        , int(xb[1])         ])/ dr  # xb[0]:=min, xb[1]:=max
jry = np.array([ int(ybr_max-yb[1]-dy), int(ybr_max-yb[0]) ])/ dr

print(' irx = {}, iry = {}'.format(irx, jry))

# Create sub-region of the raster domain. This should match the NETCDF yx-domain.
Rsub = R[jry[0]:jry[1]:clr, irx[0]:irx[1]:clr]
Rsub_dims = np.shape( Rsub )
if( not (yx_dims==(Rsub_dims)).all ):
  print(' xy-dimensions do not match: nc={} vs. r={}. Exiting ...'.format(yx_dims, Rsub_dims))
  sys.exit(1)


# Create mask array m(z,y,x)
m = np.zeros( u_dims[1:], 'uint8') # u_dims[1:] := (z_dim, y_dim, x_dim)

# Copy raster data onto each z-plane. NOTE: y-direction is reversed.
for i in xrange(z_dim):
  m[i,:,:] = Rsub[::-1,:]

# The mask data R, by default, may contain values 0 and >0. It has to be converted into
# a proper mask data [0,1]:
m[m>0] = 1 

mv = createNetcdfVariable( dso, m, 'mask', 1,      ' ',   'i4',('z','y','x',) , variable )
m = None

# To finalize, the NETCDF variables need to be copied to the new file.

uv = createNetcdfVariable( dso, u, 'u', t_dim, 'm/s', 'f4',('time','z','y','x',) , variable )
u  = None

v, v_dims = read3DVariableFromDataset( 'v', ds, varList, 0, 0, cl ) # All values.
vv = createNetcdfVariable( dso, v, 'v', t_dim, 'm/s', 'f4',('time','z','y','x',) , variable )
v  = None

w, w_dims = read3DVariableFromDataset( 'w', ds, varList, 0, 0, cl ) # All values.
wv = createNetcdfVariable( dso, w, 'w', t_dim, 'm/s', 'f4',('time','z','y','x',) , variable )
w  = None

# - - - - Done , finalize the output - - - - - - - - - -

netcdfWriteAndClose( dso )
