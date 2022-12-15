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
parser.add_argument('-i', '--invert',action="store_true", default=False,
                    help='Output also the inverse of the Reynolds stress '
                    'tensor. If partial input data is given, symmetry of the '
                    'Reynolds stress tensor is assumed.')
parser.add_argument('-n', '--missval',type=float, help='Value for missing '
                    'values in output file. Default = NaN.', default = None)

args = parser.parse_args()

#=inputs######================================================================#

vels=dict()
first = True

for i in args.files:
    ds, vD, uD = netcdfDataset2(i)
    for j in ['u','v','w','uu','vv','ww','uv','uw','vw','vu','wu','wv']:
        if ( j not in vels and j in vD.keys() ):
           vels[j] = ds[j][:,:,:,:].data
           vels[j][np.isclose(vels[j],-9999.0)] = np.nan
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

dso = netcdfOutputDataset( args.fileout )

tv = createNetcdfVariable( 
    dso, t, 'time', len(t), udt, 'f4', ('time',), True )

xv = createNetcdfVariable( 
    dso, x, 'x' , len(x), udx, 'f4', ('x',), True )

yv = createNetcdfVariable( 
    dso, y, 'y' , len(y), udy, 'f4', ('y',), True )

zv = createNetcdfVariable( 
    dso, z, 'z' , len(z), udz, 'f4', ('z',), True )

print(' Calculating the Reynolds stress tensor.')

for i in ['u', 'v', 'w']:
    for j in ['u', 'v', 'w']:
        if i in vels and j in vels and i+j in vels and 'R'+i+j not in vels:
            vels['R'+i+j] = vels[i+j]-vels[i]*vels[j]

            if args.missval != None:
                vels['R'+i+j][np.isnan(vels['R'+i+j])] = args.missval

            createNetcdfVariable(
                dso, vels['R'+i+j] , 'R'+i+j, None , 'm2/s2', 'f4',
                ('time', 'z', 'y', 'x', ), False, fill_value =-9999.0  )

if args.invert:
    print(' Inverting the Reynolds stress tensor.')
    
    # Check if Reynolds stress tensor input is symmetrical
    if 'Ruv' in vels and 'Rvu' in vels:        
        symmetric = False
    elif 'Ruw' in vels and 'Rwu' in vels:
        symmetric = False
    elif 'Rvw' in vels and 'Rwv' in vels:
        symmetric = False
    else:
        symmetric = True

    # Fill missing elements
    for i in ['u', 'v', 'w']:
        for j in ['u', 'v', 'w']:
            if 'R'+i+j not in vels:
                if 'R'+j+i in vels:
                    vels['R'+i+j] = vels['R'+j+i]
                else:
                    vels['R'+i+j] = 0.0
                    print(' *** Warning: Using R'+i+j+'=0 in the calculation '
                          'of the inverse Reynolds stress tensor.')

                
    # Analytical inverse from Wikipedia:
    # https://en.wikipedia.org/wiki/Invertible_matrix#Inversion_of_3_%C3%97_3_matrices
    vels['Luu'] = vels['Rvv']*vels['Rww']-vels['Rwv']*vels['Rvw']
    vels['Luv'] = vels['Rwv']*vels['Ruw']-vels['Ruv']*vels['Rww']
    vels['Luw'] = vels['Ruv']*vels['Rvw']-vels['Rvv']*vels['Ruw']
    det = ( vels['Ruu']*vels['Luu']
            + vels['Rvu']*vels['Luv']
            + vels['Rwu']*vels['Luw'] )

    # Set non-invertible tensors to nan.
    det[np.isclose(det,0.0)] = np.nan       
    vels['Luu'] /= det
    vels['Luv'] /= det
    vels['Luw'] /= det
    if not symmetric:
        vels['Lvu'] = ( vels['Rwu']*vels['Rvw'] - vels['Rvu']*vels['Rww'] )/det
        vels['Lwu'] = ( vels['Rvu']*vels['Rwv'] - vels['Rwu']*vels['Rvv'] )/det
        vels['Lwv'] = ( vels['Rwu']*vels['Ruv'] - vels['Ruu']*vels['Rwv'] )/det

    vels['Lvv'] = ( vels['Ruu']*vels['Rww'] - vels['Rwu']*vels['Ruw'] )/det
    vels['Lvw'] = ( vels['Rvu']*vels['Ruw'] - vels['Ruu']*vels['Rvw'] )/det        
    vels['Lww'] = ( vels['Ruu']*vels['Rvv'] - vels['Rvu']*vels['Ruv'] )/det

    for i in ['u', 'v', 'w']:
        for j in ['u', 'v', 'w']:
            if 'L'+i+j in vels:
                if args.missval != None:
                    vels['L'+i+j][np.isnan(vels['L'+i+j])] = args.missval

                
                createNetcdfVariable(
                    dso, vels['L'+i+j] , 'L'+i+j, None , 's2/m2', 'f4',
                    ('time', 'z', 'y', 'x', ), False, fill_value =-9999.0  )

netcdfWriteAndClose( dso )


#                                     |_|/                                    #
#===================================|_|_|\====================================#
#                                     |                                       #
