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
parser.add_argument('-t', '--tolerance',type=float, help='If inverting the '
                    'Reynolds stress tensor, the tolerance will be used on '
                    'the determinant to determine if it can be inverted.',
                    default = 1e-8)
parser.add_argument('-r', '--regularisation',choices=['nan','mean'],type=str,
                    help='If inverting the Reynolds stress tensor, specify '
                    'non-invertible tensors are regularised. nan: Set to NaN.'
                    ' mean: Discard non-invertible tensors and create a new '
                    'one as a mean of self and up to 26 nearest neighbours.')
args = parser.parse_args()

#=inputs######================================================================#

vels=dict()
first = True

for fi in args.files:
    ds, vD, uD = netcdfDataset2(fi)
    for vstr in ['u', 'v', 'w', 'uu', 'vv', 'ww', 'uv', 'uw', 'vw', 'vu', 'wu', 'wv']:
        if ( vstr not in vels and vstr in vD.keys() ):            
           vels[vstr] = ds[vstr][:,:,:,:]
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

for vi in ['u', 'v', 'w']:
    for vj in ['u', 'v', 'w']:
        if vi in vels and vj in vels and vi+vj in vels and 'R'+vi+vj not in vels:
            vels['R'+vi+vj] = vels[vi+vj]-vels[vi]*vels[vj]

if args.invert:
    print(' Inverting the Reynolds stress tensor.')
    print('  Tolerance for non-invertible tensors: '+str(args.tolerance))
    
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
    for vi in ['u', 'v', 'w']:
        for vj in ['u', 'v', 'w']:
            if 'R'+vi+vj not in vels:
                if 'R'+vj+vi in vels:
                    vels['R'+vi+vj] = vels['R'+vj+vi]
                else:
                    vels['R'+vi+vj] = 0.0
                    print(' *** Warning: Using R'+vi+vj+'=0 in the calculation '
                          'of the inverse Reynolds stress tensor.')

                
    # Analytical inverse from Wikipedia:
    # https://en.wikipedia.org/wiki/Invertible_matrix#Inversion_of_3_%C3%97_3_matrices
    vels['Luu'] = vels['Rvv']*vels['Rww']-vels['Rwv']*vels['Rvw']
    vels['Luv'] = vels['Rwv']*vels['Ruw']-vels['Ruv']*vels['Rww']
    vels['Luw'] = vels['Ruv']*vels['Rvw']-vels['Rvv']*vels['Ruw']
    det = ( vels['Ruu']*vels['Luu']
            + vels['Rvu']*vels['Luv']
            + vels['Rwu']*vels['Luw'] )

    # Regularisation
    if args.regularisation == 'mean':
        print('  Regularising Reynolds stresses using mean filter.')
        # Find elements that need regularisation
        Ri = np.isclose(det,0.0,atol=args.tolerance)
        for vi in ['u', 'v', 'w']:
            for vj in ['u', 'v', 'w']:
                Ra = vels['R'+vi+vj].data
                Ra[vels['R'+vi+vj].mask] = np.nan 
                Ra = np.pad(Ra , ((0, 0), (1,1), (1,1), (1,1)) , 'constant', constant_values=np.nan)
                for dk in np.argwhere(Ri):
                    vels['R'+vi+vj][dk[0],dk[1],dk[2],dk[3]] = np.nanmean(
                        Ra[dk[0],(dk[1]-1+1):(dk[1]+2+1),(dk[2]-1+1):(dk[2]+2+1),(dk[3]-1+1):(dk[3]+2+1)])
                del Ra

        # Recalculate earlier inverses.
        vels['Luu'] = vels['Rvv']*vels['Rww']-vels['Rwv']*vels['Rvw']
        vels['Luv'] = vels['Rwv']*vels['Ruw']-vels['Ruv']*vels['Rww']
        vels['Luw'] = vels['Ruv']*vels['Rvw']-vels['Rvv']*vels['Ruw']
        det = ( vels['Ruu']*vels['Luu']
                + vels['Rvu']*vels['Luv']
                + vels['Rwu']*vels['Luw'] )                                
    else:
        # Set non-invertible tensors to nan.
        print('  Non-invertible Reynolds stress tensors are set to nan.')
        Ri = np.isclose(det,0.0,atol=args.tolerance)
        det[Ri] = np.nan       
        for vi in ['u', 'v', 'w']:
            for vj in ['u', 'v', 'w']:
                vels['R'+vi+vj][Ri] = np.nan

    vels['Luu'] /= det
    vels['Luv'] /= det
    vels['Luw'] /= det

    vels['Lvv'] = ( vels['Ruu']*vels['Rww'] - vels['Rwu']*vels['Ruw'] )/det
    vels['Lvw'] = ( vels['Rvu']*vels['Ruw'] - vels['Ruu']*vels['Rvw'] )/det        
    vels['Lww'] = ( vels['Ruu']*vels['Rvv'] - vels['Rvu']*vels['Ruv'] )/det
    
    if not symmetric:
        vels['Lvu'] = ( vels['Rwu']*vels['Rvw'] - vels['Rvu']*vels['Rww'] )/det
        vels['Lwu'] = ( vels['Rvu']*vels['Rwv'] - vels['Rwu']*vels['Rvv'] )/det
        vels['Lwv'] = ( vels['Rwu']*vels['Ruv'] - vels['Ruu']*vels['Rwv'] )/det
    elif symmetric:
        for i in ['vu', 'wu', 'wv']:
            if 'R'+i in vels:
                del vels['R'+i]
    else:
        sys.exit(' Something went wrong with the symmetries. Exiting.')


    for vi in ['u', 'v', 'w']:
        for vj in ['u', 'v', 'w']:
            if 'L'+vi+vj in vels:
                createNetcdfVariable(
                    dso, vels['L'+vi+vj] , 'L'+vi+vj, None , 's2/m2', 'f4',
                    ('time', 'z', 'y', 'x', ), False )
    

for vi in ['u', 'v', 'w']:
    for vj in ['u', 'v', 'w']:
        if 'R'+vi+vj in vels:
            createNetcdfVariable(
                dso, vels['R'+vi+vj] , 'R'+vi+vj, None , 'm2/s2', 'f4',
                ('time', 'z', 'y', 'x', ), False )


                
netcdfWriteAndClose( dso )


#                                     |_|/                                    #
#===================================|_|_|\====================================#
#                                     |                                       #
