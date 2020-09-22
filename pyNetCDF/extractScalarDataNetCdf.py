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
parser = argparse.ArgumentParser(prog='extractScalarDataNetCDF.py')
parser.add_argument("-f", "--filename",type=str,\
  help="Name of the input NETCDF file.")
parser.add_argument("-fo", "--fileout",type=str, default="Scalar.nc", \
  help="Name of the output NETCDF file.")
parser.add_argument("-s", "--scalars",type=str, nargs='+', required=True,\
  help="Name of the NETCDF scalar in the file. e.g. e, p, pt")
parser.add_argument("-d", "--decomp", action="store_true", default=False,\
  help="Decomposed into mean (V_m) and fluctuating (V^prime) components.")
parser.add_argument("-cp", "--copyOnly", action="store_true", default=False,\
  help="Copy the scalar variable to a new file without interpolation.")
parser.add_argument("-nt", "--ntimeskip", type=int, default=0,\
  help="Skip <nt> number of time steps.")
parser.add_argument("-c", "--coarse", type=int, help="Coarsening level. Int > 1.", default=1) 
args = parser.parse_args() 
#==========================================================#
# Initial renaming operations and variable declarations

filename = args.filename
fileout  = args.fileout
cl       = abs(int(args.coarse))
nt       = args.ntimeskip
scalarNames = args.scalars
decompOn = args.decomp
copyOnly = args.copyOnly

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
sN = scalarNames[0] 

time, time_dims = read1DVariableFromDataset('time', sN, ds, nt, 0, 1 ) # All values.
print(' time_dims = {} '.format(time_dims))
tv = createNetcdfVariable( dso, time,'time', len(time),'s','f4',('time',), parameter )
time = None  

if( copyOnly ): ex = 0
else:           ex = 1

x, x_dims = read1DVariableFromDataset('x', sN, ds, 0, ex, cl ) # Scalar grid's x, exclude the last value.
xv = createNetcdfVariable( dso, x, 'x', len(x), 'm', 'f4', ('x',), parameter )
x = None

y, y_dims = read1DVariableFromDataset('y', sN, ds, 0, ex, cl ) # Scalar grid's y, exclude the last value.
yv = createNetcdfVariable( dso, y, 'y', len(y), 'm', 'f4', ('y',), parameter )
y = None

z, z_dims = read1DVariableFromDataset('z', sN, ds, ex, 0, cl ) # Scalar grid's z, exclude the first value.
zv = createNetcdfVariable( dso, z, 'z', len(z), 'm', 'f4', ('z',), parameter )
z = None

# - - - - Scalar components - - - - - - - - - -
sv = []
for sname in scalarNames:
  
  s0, s0_dims = read3DVariableFromDataset( sname, ds,  nt, 0, 0, cl ) # All values.
  print(' Ref: z.shape = {}, y.shape = {}, x.shape = {} '.format(z_dims,y_dims,x_dims) )
  print(' Orig: s0.shape = {} '.format(s0.shape) )
  
  idx = np.isnan( s0 ); s0[idx] = 0.
  idx = ( np.abs(s0) > 10.**9 ); s0[idx] = 0.
  

  if( not copyOnly ):
    st_dims  = np.array( s0_dims )  # Change to numpy array for manipulation
    st_dims[1:] -= 1   # Reduce the coord. dimensions by one. Note: time = uc_dims[0].

  if( decompOn and (not copyOnly) ):
    #s = np.zeros( st_dims )   # This is not really needed.
    s, sm = interpolatePalmVectors( s0, s0_dims, 'k' , decompOn ); #s0 = None
    sp = vectorPrimeComponent( s, sm )
    sv.append(createNetcdfVariable( dso, sp, sname+'p', st_dims[0],'','f4',('time','z','y','x',), variable ))
  else:
    # Take the portion that matches the coords.
    if( copyOnly ): s = s0
    else: s = s0[:,1:,:-1,:-1].copy(); s0 = None
    
  s_dims = np.shape(s)
  print(' Final shape = {} '.format(s_dims) )
  sv.append( createNetcdfVariable(dso, s, sname, s_dims[0],'[]','f4',('time','z','y','x',), variable) )

netcdfWriteAndClose( dso )
