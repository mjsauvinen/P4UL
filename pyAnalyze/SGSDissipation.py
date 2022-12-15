#!/usr/bin/env python3
import argparse
import sys
from netcdfTools import *

'''Calculate SGS dissipation
  
Calculate the SGS dissipation in a PALM run according to the Deardorf SGS
model.

Author:
Jukka-Pekka Keskinen
Finnish Meteorological Insitute
11/2022

'''

#=argument parser=============================================================#

parser = argparse.ArgumentParser(
    prog='SGSDissipation.py',
    description="Calculate the SGS dissipation (ϵ) using Deardorff's "
    "approach.")
parser.add_argument('-f', '--filename',type=str, 
                    help='Name of the input data file. It has to contain '
                    'SGS TKE (e) and local SGS eddy diffusivity (K_m). ')
parser.add_argument('-fo', '--fileout',type=str, help='Name of output file.',
                    default = 'epsilon.nc')
parser.add_argument('-n', '--missval',type=float, help='Value for missing '
                    'values in output file. Default = NaN.', default = None)

args = parser.parse_args()

#=inputs######================================================================#

ds, vD, uD = netcdfDataset2(args.filename)

for i in ['e', 'km']:
    if (i not in vD.keys() ):
        sys.exit(
            '{} not found from variable list: {}'.format(
                i, vD.keys()))

for i in ['x', 'y', 'z', 'time']:
    if (i not in ds.dimensions.keys() ):
        sys.exit(
            '{} not found from dimension list: {}'.format(
                i, ds.dimensions.keys()))
        
Cm = 0.1   # SGS model constant
        

#=calculations================================================================#

print('   Calculating SGS dissipation.')

x = ds['x'][:].data
y = ds['y'][:].data
z = ds['z'][:].data
t = ds['time'][:].data

dx = x[1] - x[0]
dy = y[1] - y[0]
dz = np.append(z[0],z[1:] - z[:-1])
e = ds['e'][:,:,:,:].data
e[np.isclose(e,-9999.0)] = np.nan
km = ds['km'][:,:,:,:]
km[np.isclose(km,-9999.0)] = np.nan

muoto = e.shape
delta = np.broadcast_to(
    np.reshape(
        np.broadcast_to(
            np.reshape(
                np.minimum((dx*dy*dz)**(1/3),1.8*z),
                (muoto[1],1)),
            (muoto[1],muoto[2])),
        (muoto[1],muoto[2],1)), muoto)

epsi = (0.19*Cm*(e**2)/km + 0.74*(e**(3/2))/delta)

if args.missval != None:
    epsi[np.isnan(epsi)] = args.missval

#=output======================================================================#

   
dso = netcdfOutputDataset( args.fileout )

tv = createNetcdfVariable( 
    dso, t, 'time', len(t), uD['time'], 'f4', ('time',), True )

xv = createNetcdfVariable( 
    dso, x, 'x' , len(x), uD['x'], 'f4', ('x',), True )

yv = createNetcdfVariable( 
    dso, y, 'y' , len(y), uD['y'], 'f4', ('y',), True )

zv = createNetcdfVariable( 
    dso, z, 'z' , len(z), uD['z'], 'f4', ('z',), True )

Ev = createNetcdfVariable(dso, epsi , 'ee' , None , 'm2/s3', 'f4',
                          ('time', 'z', 'y', 'x', ), False,
                          fill_value =-9999.0 )

netcdfWriteAndClose( dso )

#                                     |_|/                                    #
#===================================|_|_|\====================================#
#                                     |                                       #
