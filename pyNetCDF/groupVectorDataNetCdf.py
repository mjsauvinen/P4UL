#!/usr/bin/env python3
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
def replaceValues( q, qa, qb ):
  if( qa[0] != 1e9 ):
    idr = (q > qa[0])
    q[idr] = qa[1]
  if( qb[0] != -1e9 ):
    idr = (q < qb[0])
    q[idr] = qb[1]
    
  return q

#==========================================================#
parser = argparse.ArgumentParser(prog='groupVectorDataNetCDF.py')
parser.add_argument("-f", "--filename",type=str,\
  help="Name of the input NETCDF file.")
parser.add_argument("-fo", "--fileout",type=str, default="U-VECTOR.nc",\
  help="Name of the output NETCDF file.")
parser.add_argument("-sn", "--scalars",type=str, nargs='+', default=None,\
  help="(Optional) Scalars to be included.")
parser.add_argument("-sx", "--suffix",type=str, default='',\
  help="Potential suffix to be appended to variable names. Example: '_xy'. ")
parser.add_argument("-va", "--replValuesAbove", nargs=2, type=float, default=[1.e9,0.0],\
  help="Entry <Max> <val> := replace values above given threshold <Max> by <val>. Default= 1.e9 0.0")
parser.add_argument("-vb", "--replValuesBelow", nargs=2, type=float, default=[-1.e9,0.0],\
  help="Entry <Min> <val> := replace values below given threshold <Min> by <val>. Default= -1.e9 0.0")
parser.add_argument("-d", "--decomp", action="store_true", default=False,\
  help="Decomposed into mean (V_m) and fluctuating (V^prime) components.")
parser.add_argument("-dd", "--decompOnly",action="store_true", default=False,\
  help="Output V_m and V^prime components only.")
parser.add_argument("-nt", "--ntimeskip", type=int, default=0,\
  help="Skip <nt> number of time steps.")
parser.add_argument("-c", "--coarse", type=int, default=1,\
  help="Coarsening level. Int > 1.")
parser.add_argument("-kc", "--kcopy",action="store_true", default=False,\
  help="Copy in z-direction, do not interpolate.")
args = parser.parse_args()
#==========================================================#
# Initial renaming operations and variable declarations

filename   = args.filename
fileout    = args.fileout
scalars    = args.scalars
suffix     = args.suffix
ntskip     = args.ntimeskip
cl         = abs(int(args.coarse))
kcopy      = args.kcopy
va         = args.replValuesAbove
vb         = args.replValuesBelow

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
time, time_dims = read1DVariableFromDataset('time','u'+suffix, ds, ntskip, 0, 1 ) # All values.
tv = createNetcdfVariable( dso, time,'time', len(time),'s','f4',('time',), parameter )
time = None

x, x_dims = read1DVariableFromDataset( 'x','w'+suffix, ds, 0, 1, cl ) # Exclude the last value.
xv = createNetcdfVariable( dso, x   , 'x'   , len(x)   , 'm', 'f4', ('x',)   , parameter )
x = None

y, y_dims = read1DVariableFromDataset( 'y','u'+suffix, ds, 0, 1, cl ) # Exclude the last value.
print(' y_dims = {} '.format(y_dims))
y[np.isnan(y)] = 0.  # Special treatment.
yv = createNetcdfVariable( dso, y   , 'y'   , len(y)   , 'm', 'f4', ('y',)   , parameter )
y = None

if( kcopy ): xk = 0
else:        xk = 1
z, z_dims = read1DVariableFromDataset('z','u'+suffix, ds, xk, 0, cl )
zv = createNetcdfVariable( dso, z   , 'z'   , len(z)   , 'm', 'f4', ('z',)   , parameter )
print(' z_dims = {} '.format(z_dims))
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
u0, u0_dims = read3DVariableFromDataset( 'u'+suffix, ds, ntskip, 0, 0, cl ) # All values.
u0 = replaceValues(u0, va, vb)

''' 
New, cell-center dimension lengths: 
Number of times remains the same, but coord. lengths 
are reduced by one due to interpolation.
'''
cc_dims  = np.array( u0_dims )  # Change to numpy array for manipulation
if( kcopy ): cc_dims[2:] -= 1   # Reduce the x, y coord. dimensions by one. Note: time = uc_dims[0].
else:        cc_dims[1:] -= 1   # Reduce all coord. dimensions by one.
print(' u0_dims = {}, cc_dims = {} '.format(u0_dims,cc_dims))


uc = np.zeros( cc_dims )
uc, um = interpolatePalmVectors( u0, cc_dims, 'i' , decompOn ); u0 = None

if( not args.decompOnly ):
  uv = createNetcdfVariable( dso, uc, 'u', cc_dims[0], 'm/s', 'f4',('time','z','y','x',) , variable )
  if( not decompOn ): uc = None
if( decompOn ):
  up = vectorPrimeComponent( uc, um ); uc = None
  upv = createNetcdfVariable( dso, up, 'up', cc_dims[0], 'm/s', 'f4',('time','z','y','x',) , variable )
  up = None
  umv = createNetcdfVariable( dso, um, 'um', cc_dims[0], 'm/s', 'f4',('z','y','x',) , variable )
  um = None


# - - - - Third, w-component - - - - - - - - - -

w0, w0_dims = read3DVariableFromDataset( 'w'+suffix, ds, ntskip, 0, 0, cl ) # All values.
w0 = replaceValues(w0, va, vb)

wc = np.zeros( cc_dims )
if( kcopy ):
  wc, wm = interpolatePalmVectors( w0, cc_dims, 'kc' , decompOn ); w0 = None
else:
  wc, wm = interpolatePalmVectors( w0, cc_dims, 'k' , decompOn ); w0 = None

if( not args.decompOnly ):
  wv = createNetcdfVariable( dso, wc, 'w', cc_dims[0], 'm/s', 'f4',('time','z','y','x',) , variable )
  if( not decompOn ): wc = None  # ASAP
if( decompOn ):
  wp = vectorPrimeComponent( wc, wm )
  wpv = createNetcdfVariable( dso, wp, 'wp', cc_dims[0], 'm/s', 'f4',('time','z','y','x',) , variable )
  wp = None
  wmv = createNetcdfVariable( dso, wm, 'wm', cc_dims[0], 'm/s', 'f4',('z','y','x',) , variable )
  wm = None


# - - - - Second, v-component - - - - - - - - - -
'''
v0, v0_dims = read3DVariableFromDataset( 'v'+suffix, ds, ntskip, 0, 0, cl ) # All values.
v0 = replaceValues(v0, va, vb)

vc = np.zeros( cc_dims )
vc, vm = interpolatePalmVectors( v0, cc_dims, 'j' , decompOn ); v0 = None

if( not args.decompOnly ):
  vv = createNetcdfVariable( dso, vc, 'v', cc_dims[0], 'm/s', 'f4',('time','z','y','x',) , variable )
  if( not decompOn ): vc = None
if( decompOn ):
  vp = vectorPrimeComponent( vc, vm ); vc = None
  vpv = createNetcdfVariable( dso, vp, 'vp', cc_dims[0], 'm/s', 'f4',('time','z','y','x',) , variable )
  vp = None
  vmv = createNetcdfVariable( dso, vm, 'vm', cc_dims[0], 'm/s', 'f4',('z','y','x',) , variable )
  vm = None
'''

# - - - - Fouth, possible scalars - - - - - - - - - -

if( scalars ):
  for sn in scalars:
    s0, s0_dims = read3DVariableFromDataset( sn+suffix, ds, ntskip, 0, 0, cl ) # All values.
    s0 = replaceValues(s0, va, vb)
    
    sc_dims  = np.array( s0_dims )  # Change to numpy array for manipulation
    if( kcopy ): sc_dims[2:] -= 1   # Reduce the x, y dimensions by one. Note: time = sc_dims[0].
    else:        sc_dims[1:] -= 1   # Reduce all coord. dimensions by one.
    sc = np.zeros( sc_dims )
    sc, sm = interpolatePalmVectors( s0, s0_dims, 'i' , decompOn ); s0 = None
    if( not args.decompOnly ):
      sv = createNetcdfVariable( dso, sc, sn, sc_dims[0], ' ', 'f4',('time','z','y','x',) , variable )
      if( not decompOn ): sc = None
    if( decompOn ):
      sp = vectorPrimeComponent( sc, sm ); sc = None
      spv = createNetcdfVariable( dso, sp, sn+'p', sc_dims[0], ' ', 'f4',('time','z','y','x',) , variable )
      sp = None
      smv = createNetcdfVariable( dso, sm, sn+'m', sc_dims[0], ' ', 'f4',('z','y','x',) , variable )
      sm = None

# - - - - Done , finalize the output - - - - - - - - - -
netcdfWriteAndClose( dso )
