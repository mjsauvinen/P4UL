#!/usr/bin/env python3
import argparse
import sys
from netcdfTools import *

'''Calculate Reynolds stresses
  
Calculate the Reynolds stresses from time averaged velocity fields.

Author:
Jukka-Pekka Keskinen
Finnish Meteorological Insitute
11/2022

'''

#=argument parser=============================================================#

parser = argparse.ArgumentParser(
    prog='reynoldsStressNetCdf2.py',
    description="Calculate Reynolds stresses (properly) using time averaged "
    "velocity fields (u, v, w, uu, vv, ww uv, uw, vw). Collocated fields are "
    "expected. ")
parser.add_argument('-f', '--files',type=str, nargs='+',
                    help='Name of the input data files. These should contain '
                    'some or all of the time-averaged quantities u, v, w, uu, '
                    'vv, ww, uv, uw, and vw. These are used based on '
                    'availability. These fields need to be collocated.')
parser.add_argument('-fo', '--fileout',type=str, help='Name of output file.',
                    default = 'R.nc')

args = parser.parse_args()

#=inputs######================================================================#

vels=dict()
first = True

for i in args.files:
    ds, vD, uD = netcdfDataset2(i)
    for j in ['u', 'v', 'w', 'uu', 'vv', 'ww', 'uv', 'uw', 'vw']:
        if ( j not in vels and j in vD.keys() ):
           vels[j] = ds[j][:,:,:,:].data
           if first:
               x = ds['x'][:].data
               udx = uD['x']
               y = ds['y'][:].data
               udy = uD['y']
               z = ds['z'][:].data
               udz = uD['z']
               t = ds['time'][:].data
               udt = uD['time']
               first = False
               
#=calculations and output======================================================#

for i in ['u', 'v', 'w']:
    for j in ['u', 'v', 'w']:
        if i in vels and j in vels and i+j in vels:
            vels['R'+i+j] = vels[i+j]-vels[i]*vels[j]


dso = netcdfOutputDataset( args.fileout )

tv = createNetcdfVariable( 
    dso, t, 'time', len(t), udt, 'f4', ('time',), True )

xv = createNetcdfVariable( 
    dso, x, 'x' , len(x), udx, 'f4', ('x',), True )

yv = createNetcdfVariable( 
    dso, y, 'y' , len(y), udy, 'f4', ('y',), True )

zv = createNetcdfVariable( 
    dso, z, 'z' , len(z), udz, 'f4', ('z',), True )

for i in ['u', 'v', 'w']:
    for j in ['u', 'v', 'w']:
        if i in vels and j in vels and i+j in vels:
            vels['R'+i+j] = vels[i+j]-vels[i]*vels[j]
            createNetcdfVariable(
                dso, vels['R'+i+j] , 'R'+i+j, None , 'm2/s2', 'f4',
                ('time', 'z', 'y', 'x', ), False )

netcdfWriteAndClose( dso )


#                                     |_|/                                    #
#===================================|_|_|\====================================#
#                                     |                                       #
