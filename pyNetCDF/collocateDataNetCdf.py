#!/usr/bin/env python3
from netcdfTools import *
import sys
import argparse
import numpy as np
''' 
Description: 
  Script to interpolate PALM 3d netcdf data (vectors and scalars) onto scalar grid (i.e. cell centers).
  This script is destined to replace groupVectorDataNetCDF.py, whose operation is too limited.

Author: Mikko Auvinen
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
parser = argparse.ArgumentParser(prog='collocateDataNetCdf.py')
parser.add_argument("-f", "--filename",type=str,\
  help="Name of the input NETCDF file.")
parser.add_argument("-fo", "--fileout",type=str, default="V_3D.nc",\
  help="Name of the output NETCDF file. Default = V_3D.nc ")
parser.add_argument("-vn", "--vectorNames", type=str, nargs=3, default=['u','v','w'],\
  help="Names of the vector components. Default = u v w")
parser.add_argument("-sn", "--scalarNames",type=str, nargs='+', default=None,\
  help="(Optional) Scalars to be included.")
parser.add_argument("-dn", "--derivNames",type=str, nargs='+', default=None,\
  help="(Optional) Names of derived coordinates to output to same file.")
parser.add_argument("-va", "--replValuesAbove", nargs=2, type=float, default=[1.e9,0.0],\
  help="Entry <Max> <val> := replace values above given threshold <Max> by <val>. Default= 1.e9 0.0")
parser.add_argument("-vb", "--replValuesBelow", nargs=2, type=float, default=[-1.e9,0.0],\
  help="Entry <Min> <val> := replace values below given threshold <Min> by <val>. Default= -1.e9 0.0")
parser.add_argument("-d", "--decomp", action="store_true", default=False,\
  help="Decomposed into mean (V_m) and fluctuating (V^prime) components.")
parser.add_argument("-sx", "--suffix",type=str, default='',\
  help="Potential suffix to be appended to variable names. Example: '_xy'. ")
helpk2z='''Conversion factor for converting vertical coord from index (like ku_above_surf) 
to elevation in meters. Example usage: --k2z 2 means vertical resolution is 2 m.'''
parser.add_argument("-k2z", "--k2z",type=float, default=None,\
  help=helpk2z)
parser.add_argument("-dd", "--decompOnly",action="store_true", default=False,\
  help="Output V_m and V^prime components only.")
parser.add_argument("-nt", "--ntimeskip", type=int, default=0,\
  help="Skip <nt> number of time steps.")
parser.add_argument("-c", "--coarse", type=int, default=1,\
  help="Coarsening level (int). Ex: c=2 means every other data point.")
parser.add_argument("-kc", "--kcopy",action="store_true", default=False,\
  help="Copy in z-direction, do not interpolate.")
args = parser.parse_args()
#==========================================================#
# Initial renaming operations and variable declarations

filename    = args.filename
fileout     = args.fileout
vn          = args.vectorNames
sn         = args.scalarNames
dn         = args.derivNames
ntskip     = args.ntimeskip
cl         = abs(int(args.coarse))
kcopy      = args.kcopy
va         = args.replValuesAbove
vb         = args.replValuesBelow
k2z        = args.k2z
suffix     = args.suffix

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
ds, vD, uD = netcdfDataset2(filename) # vD: variableDict, uD: unitDict

# Create a NETCDF output dataset (dso) for writing out the data.
dso = netcdfOutputDataset( fileout )

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# Check that vector components are found
for vi in vn:
  if( vi not in vD.keys() ):
    sys.exit(' Vector component {} not found from variable list: {}'.format(vi, vD.keys()))

# Same for scalars and derived coords
if( sn ):
  for si in sn:
    if( si not in vD.keys() ):
      sys.exit(' Scalar {} not found from variable list: {}'.format(si, vD.keys()))

if( dn ):
  for di in dn:
    if( di not in vD.keys() ):
      sys.exit(' Derived coord {} not found from variable list: {}'.format(di, vD.keys()))

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# Let's construct the desired parameter list, which may not be (time, z, y, x) always
vx = vn[0]; vy = vn[1]
tn = vD[vx][0]; zn = vD[vx][1]; yn = vD[vx][2]; xn = vD[vy][3]
zunit = uD[zn] # Store unit for z
if('z' == zn[0]): zn = 'z'

'''
Read cell center coordinates and time. 
In PALM the output cell center grid is constructed using (x, y, zu_3d) from the following:
  u('time', 'zu_3d', 'y' , 'xu')  <- take y and zu_3d from here
  v('time', 'zu_3d', 'yv', 'x')   <- take x from here
  w('time', 'zw_3d', 'y' , 'x')
Create the output independent variables right away and empty memory.
'''
time, time_dims = read1DVariableFromDataset(tn, vn[0] , ds, ntskip, 0, 1 ) # All values.
tv = createNetcdfVariable( dso, time, tn, len(time), uD[tn],'f4', (tn,), parameter )

x, x_dims = read1DVariableFromDataset( xn, vn[1], ds, 0, 1, cl ) # Exclude the last value.
xv = createNetcdfVariable( dso, x  , xn , len(x)   , uD[xn],'f4', (xn,), parameter )

y, y_dims = read1DVariableFromDataset(yn, vn[0], ds, 0, 1, cl ) # Exclude the last value.
yv = createNetcdfVariable( dso, y  , yn , len(y)   , uD[yn],'f4', (yn,), parameter )

if( kcopy ): xk = 0
else:        xk = 1
z, z_dims = read1DVariableFromDataset(zn, vn[0], ds, xk, 0, cl ) # Exclude determined by xk
if(('k' in zn)  and (k2z is not None) ):
  print(' NOTE: Converting vertical coord from index {} to z [m] using dz={} m.'.format(zn,k2z))
  z *= k2z;  zn = 'z'; zunit = 'meters'
zv = createNetcdfVariable( dso, z  , zn , len(z)   , zunit, 'f4', (zn,), parameter )

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# Include additional (derived) coordinates into the output file.
if( dn ):
  for di in dn:
    dc = ds.variables[di][:]
    dc_dims = np.shape( dc )
    if(   len( dc_dims ) == 1 ): 
      dc = dc[:-1]
    elif( len( dc_dims ) == 2 ): 
      dc = dc[:-1,:-1]
    else: 
      print(' Only 2d coords allowed. Skipping {} with dims={}.'.format(di,dc_dims))
      continue
    dv = createNetcdfVariable( dso, dc, di, None, uD[di], 'f4', vD[di], variable )
    dc = None

time = None; x = None; y = None; z = None

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
'''
Read in the velocity components.
Example PALM netCDF4 format: 
  u(time, zu_3d, y, xu)
  v(time, zu_3d, yv, x)
  w(time, zw_3d, y, x)  
  Interpolate u0 -> uc and output the values right away. Empty memory asap.
'''

# - - - - First, u-component - - - - - - - - - -
u0, u0_dims = read3DVariableFromDataset( vn[0], ds, ntskip, 0, 0, cl ) # All values.
u0 = replaceValues(u0, va, vb)

''' 
New, cell-center dimension lengths: 
Number of times remains the same, but coord. lengths 
are reduced by one due to interpolation.
'''
cc_dims  = np.array( u0_dims )  # Change to numpy array for manipulation
if( kcopy ): cc_dims[2:] -= 1   # Reduce the x, y coord. dimensions by one. Note: time=cc_dims[0].
else:        cc_dims[1:] -= 1   # Reduce all coord. dimensions by one.
print(' u0_dims = {}, cc_dims = {} '.format(u0_dims,cc_dims))


uc = np.zeros( cc_dims )
uc, um = interpolatePalmVectors( u0, cc_dims, 'i' , decompOn ); u0 = None

if( not args.decompOnly ):
  uv = createNetcdfVariable( dso, \
    uc, vn[0].replace(suffix,''), None, uD[vn[0]], 'f4',(tn,zn,yn,xn,) , variable )
  if( not decompOn ): uc = None

if( decompOn ):
  up = vectorPrimeComponent( uc, um ); uc = None
  
  upv = createNetcdfVariable( dso, \
    up, vn[0].replace(suffix,'')+'p', None, uD[vn[0]], 'f4',(tn,zn,yn,xn,), variable )
  up = None
  
  umv = createNetcdfVariable( dso, \
    um, vn[0].replace(suffix,'')+'m', None, uD[vn[0]], 'f4',(zn,yn,xn,), variable )
  um = None


# - - - - Third, w-component - - - - - - - - - -

w0, w0_dims = read3DVariableFromDataset( vn[2], ds, ntskip, 0, 0, cl ) # All values.
w0 = replaceValues(w0, va, vb)

wc = np.zeros( cc_dims )
if( kcopy ):
  wc, wm = interpolatePalmVectors( w0, cc_dims, 'kc' , decompOn ); w0 = None
else:
  wc, wm = interpolatePalmVectors( w0, cc_dims, 'k' , decompOn ); w0 = None

if( not args.decompOnly ):
  wv = createNetcdfVariable( dso, \
    wc, vn[2].replace(suffix,''), None, uD[vn[2]], 'f4',(tn,zn,yn,xn,) , variable )
  if( not decompOn ): wc = None  # ASAP
if( decompOn ):
  wp = vectorPrimeComponent( wc, wm )
  wpv = createNetcdfVariable( dso, \
    wp, vn[2].replace(suffix,'')+'p', None, uD[vn[2]], 'f4',(tn,zn,yn,xn,) , variable )
  wp = None
  
  wmv = createNetcdfVariable( dso, \
    wm, vn[2].replace(suffix,'')+'m', None, uD[vn[2]], 'f4',(zn,yn,xn,) , variable )
  wm = None


# - - - - Second, v-component - - - - - - - - - -
v0, v0_dims = read3DVariableFromDataset( vn[1], ds, ntskip, 0, 0, cl ) # All values.
v0 = replaceValues(v0, va, vb)

vc = np.zeros( cc_dims )
vc, vm = interpolatePalmVectors( v0, cc_dims, 'j' , decompOn ); v0 = None

if( not args.decompOnly ):
  vv = createNetcdfVariable( dso, \
    vc, vn[1].replace(suffix,''), None, uD[vn[1]], 'f4',(tn,zn,yn,xn,) , variable )
  if( not decompOn ): vc = None

if( decompOn ):
  vp = vectorPrimeComponent( vc, vm ); vc = None
  vpv = createNetcdfVariable( dso, \
    vp, vn[1].replace(suffix,'')+'p', None, uD[vn[1]],'f4',(tn,zn,yn,xn,), variable )
  vp = None
  
  vmv = createNetcdfVariable( dso, \
    vm, vn[1].replace(suffix,'')+'m', None, uD[vn[1]],'f4',(zn,yn,xn,) , variable )
  vm = None

# - - - - Fouth, possible scalars - - - - - - - - - -

if( sn ):
  for si in sn:
    # Let's construct the desired parameter list, which may not be (time, z, y, x) always
    sunit = uD[si]

    s0, s0_dims = read3DVariableFromDataset( si, ds, ntskip, 0, 0, cl ) # All values.
    s0 = replaceValues(s0, va, vb)
    
    sc_dims  = np.array( s0_dims )  # Change to numpy array for manipulation
    if( kcopy ): sc_dims[2:] -= 1   # Reduce the x, y dimensions by one. Note: time = sc_dims[0].
    else:        sc_dims[1:] -= 1   # Reduce all coord. dimensions by one.
    sc = np.zeros( sc_dims )
    #sc, sm = interpolatePalmVectors( s0, s0_dims, 'i' , decompOn ); s0 = None
    sc = s0[:,1:,:-1,:-1].copy(); s0 = None
    
    if( not args.decompOnly ):
      sv = createNetcdfVariable(dso, \
        sc, si.replace('suffix',''), None, uD[si],'f4',(tn,zn,yn,xn,), variable )
      if( not decompOn ): sc = None
    
    if( decompOn ):
      sm = np.mean( sc, axis=0 )
      sp = vectorPrimeComponent( sc, sm ); sc = None
      spv = createNetcdfVariable( dso, \
        sp, si.replace('suffix','')+'p', None, uD[si], 'f4',(tn,zn,yn,xn,) , variable )
      sp = None
      
      smv = createNetcdfVariable( dso, \
        sm, si.replace('suffix','')+'m', None, uD[si], 'f4',(zn,yn,xn,) , variable )
      sm = None



# - - - - Done , finalize the output - - - - - - - - - -
netcdfWriteAndClose( dso )
