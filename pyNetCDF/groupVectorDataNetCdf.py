#!/usr/bin/env python
from netcdfTools import *
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

#==========================================================#
parser = argparse.ArgumentParser(prog='groupVectorDataNetCDF.py')
parser.add_argument("-f", "--filename",type=str, help="Name of the input NETCDF file.")
parser.add_argument("-fo", "--fileout",type=str, help="Name of the output NETCDF file.", default="U-VECTOR.nc")
parser.add_argument("-sn", "--scalar",type=str, help="(Optional) Scalar to be included.", default=None)
parser.add_argument("-d", "--decomp", help="Decomposed into mean (V_m) and fluctuating (V^prime) components.",\
  action="store_true", default=False) 
parser.add_argument("-dd", "--decompOnly", help="Output V_m and V^prime components only.",\
  action="store_true", default=False)
parser.add_argument("-nt", "--ntimeskip", type=int, help="Skip <nt> number of time steps.",\
  default=0)
parser.add_argument("-c", "--coarse", type=int, help="Coarsening level. Int > 1.",\
  default=1) 
args = parser.parse_args() 
#==========================================================#
# Initial renaming operations and variable declarations

filename   = args.filename
fileout    = args.fileout
scalarName = args.scalar 
nt         = args.ntimeskip
cl         = abs(int(args.coarse))

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
ds, varList, paramList = netcdfDataset(filename)

# Create a NETCDF output dataset (dso) for writing out the data.
dso = netcdfOutputDataset( fileout )


'''
Read cell center coordinates and time.
Create the output independent variables right away and empty memory.
'''
time, time_dims = read1DVariableFromDataset('time', ds, nt, 0, 1 ) # All values.
tv = createNetcdfVariable( dso, time,'time', len(time),'s','f4',('time',), parameter )
time = None  

x, x_dims = read1DVariableFromDataset( 'x',ds, 0, 1, cl ) # Exclude the last value.
xv = createNetcdfVariable( dso, x   , 'x'   , len(x)   , 'm', 'f4', ('x',)   , parameter )
x = None

y, y_dims = read1DVariableFromDataset( 'y',ds, 0, 1, cl ) # Exclude the last value.
print(' y_dims = {} '.format(y_dims))
y[np.isnan(y)] = 0.  # Special treatment.
yv = createNetcdfVariable( dso, y   , 'y'   , len(y)   , 'm', 'f4', ('y',)   , parameter )
y = None

z, z_dims = read1DVariableFromDataset( 'zu_3d',ds, 1, 0, cl ) # Exclude the first value.
zv = createNetcdfVariable( dso, z   , 'z'   , len(z)   , 'm', 'f4', ('z',)   , parameter )
z = None


'''
Read in the velocity components.
PALM netCDF4: 
  u(time, zu_3d, y, xu)
  v(time, zu_3d, yv, x)
  w(time, zw_3d, y, x)  
  Interpolate u0 -> uc and output the values right away. Empty memory asap.
'''

# - - - - First, u-component - - - - - - - - - -
u0, u0_dims = read3DVariableFromDataset( 'u', ds, nt, 0, 0, cl ) # All values.

''' 
New, cell-center dimension lengths: 
Number of times remains the same, but coord. lengths 
are reduced by one due to interpolation.
'''
cc_dims  = np.array( u0_dims )  # Change to numpy array for manipulation
cc_dims[1:] -= 1   # Reduce the coord. dimensions by one. Note: time = uc_dims[0].

uc = np.zeros( cc_dims )
uc, um = interpolatePalmVectors( u0, u0_dims, 'i' , decompOn ); u0 = None

if( not args.decompOnly ):
  uv = createNetcdfVariable( dso, uc, 'u', cc_dims[0], 'm/s', 'f4',('time','z','y','x',) , variable )
  if( not decompOn ): uc = None
if( decompOn ):
  up = vectorPrimeComponent( uc, um ); uc = None
  upv = createNetcdfVariable( dso, up, 'up', cc_dims[0], 'm/s', 'f4',('time','z','y','x',) , variable )
  up = None
  umv = createNetcdfVariable( dso, um, 'um', cc_dims[0], 'm/s', 'f4',('z','y','x',) , variable )
  um = None

# - - - - Second, v-component - - - - - - - - - -

v0, v0_dims = read3DVariableFromDataset( 'v', ds, nt, 0, 0, cl ) # All values.

vc = np.zeros( cc_dims )
vc, vm = interpolatePalmVectors( v0, v0_dims, 'j' , decompOn ); v0 = None

if( not args.decompOnly ):
  vv = createNetcdfVariable( dso, vc, 'v', cc_dims[0], 'm/s', 'f4',('time','z','y','x',) , variable )
  if( not decompOn ): vc = None
if( decompOn ):
  vp = vectorPrimeComponent( vc, vm ); vc = None
  vpv = createNetcdfVariable( dso, vp, 'vp', cc_dims[0], 'm/s', 'f4',('time','z','y','x',) , variable )
  vp = None
  vmv = createNetcdfVariable( dso, vm, 'vm', cc_dims[0], 'm/s', 'f4',('z','y','x',) , variable )
  vm = None

# - - - - Third, w-component - - - - - - - - - -

w0, w0_dims = read3DVariableFromDataset( 'w', ds, nt, 0, 0, cl ) # All values.

wc = np.zeros( cc_dims )
wc, wm = interpolatePalmVectors( w0, w0_dims, 'k' , decompOn ); w0 = None

if( not args.decompOnly ):
  wv = createNetcdfVariable( dso, wc, 'w', cc_dims[0], 'm/s', 'f4',('time','z','y','x',) , variable )
  if( not decompOn ): wc = None  # ASAP
if( decompOn ):
  wp = vectorPrimeComponent( wc, wm )
  wpv = createNetcdfVariable( dso, wp, 'wp', cc_dims[0], 'm/s', 'f4',('time','z','y','x',) , variable )
  wp = None 
  wmv = createNetcdfVariable( dso, wm, 'wm', cc_dims[0], 'm/s', 'f4',('z','y','x',) , variable )
  wm = None

if( scalarName ):
  s0, s0_dims = read3DVariableFromDataset( scalarName, ds, nt, 0, 0, cl ) # All values.
  sc_dims  = np.array( s0_dims )  # Change to numpy array for manipulation
  sc_dims[1:] -= 1   # Reduce the coord. dimensions by one. Note: time = sc_dims[0].
  sc = np.zeros( sc_dims )
  sc, sm = interpolatePalmVectors( s0, s0_dims, 'i' , decompOn ); s0 = None
  if( not args.decompOnly ):
    sv = createNetcdfVariable( dso, sc, scalarName, sc_dims[0], ' ', 'f4',('time','z','y','x',) , variable )
    if( not decompOn ): sc = None
  if( decompOn ):
    sp = vectorPrimeComponent( sc, sm ); sc = None
    spv = createNetcdfVariable( dso, sp, scalarName+'p', sc_dims[0], ' ', 'f4',('time','z','y','x',) , variable )
    sp = None
    smv = createNetcdfVariable( dso, sm, scalarName+'m', sc_dims[0], ' ', 'f4',('z','y','x',) , variable )
    sm = None

# - - - - Done , finalize the output - - - - - - - - - -
netcdfWriteAndClose( dso )

