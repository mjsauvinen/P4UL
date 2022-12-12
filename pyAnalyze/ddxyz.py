#!/usr/bin/env python3
import argparse
import sys
from netcdfTools import *

'''Calculate derivatives

Calculate derivatives for a given field. 

NB! Boundaries are not treated properly.

Author:
Jukka-Pekka Keskinen
Finnish Meteorological Insitute
12/2022

'''

#=argument parser=============================================================#

parser = argparse.ArgumentParser(
    prog='ddxyz.py',
    description="Calculate derivatives for a given field.")
parser.add_argument('-f', '--filename',type=str, 
                    help='Name of the input data file. ', required=True)
parser.add_argument('-fo', '--fileout',type=str, help='Name of output file.',
                    default = 'ddx.nc')
parser.add_argument('-v', '--variable', type=str, nargs='+', help = 'Names of'
                    ' variables to be derivated.', required=True)
parser.add_argument('-x', '--ddx', action="store_true", default=False,
                    help = 'Calculate x derivative.')
parser.add_argument('-y', '--ddy', action="store_true", default=False,
                    help = 'Calculate y derivative.')
parser.add_argument('-z', '--ddz', action="store_true", default=False,
                    help = 'Calculate z derivative.')
parser.add_argument('-n', '--missToNan',action="store_true", default=False,
                    help='Set PALM missing values (-9999.0) to numpy in the'
                    'calculation of the derivatives. Default: set to 0.0.')
parser.add_argument('-cy', '--cyclicy', action="store_true", default=False,
                    help = 'Cyclic boundaries in y direction.')
parser.add_argument('-cx', '--cyclicx', action="store_true", default=False,
                    help = 'Cyclic boundaries in x direction.')
args = parser.parse_args()

#=inputs######================================================================#

ds, vD, uD = netcdfDataset2(args.filename)

for i in args.variable:
    if (i not in vD.keys() ):
        sys.exit(
            '{} not found from variable list: {}'.format(
                i, vD.keys()))

for i in ['x', 'y', 'z']:
    if (i not in ds.dimensions.keys() ):
        if ('dd'+i) in args:
            sys.exit(
                '{} not found from dimension list: {}'.format(
                    i, ds.dimensions.keys()))
        
#=calculations================================================================#

print(' Calculating derivatives.')

outddx = {}

if args.cyclicx:
    x = ds['x'][:].data
else:
    x = ds['x'][1:-1].data

if args.cyclicy:
    y = ds['y'][:].data
else:
    y = ds['y'][1:-1].data
    
# Assume uniform mesh in x and y directions
dx = x[1]-x[0]
dy = y[1]-y[0]


for i in args.variable:
    if len(vD[i]) == 4:
        A = ds[i][:,:,:,:].data
        if args.missToNan:
            A[np.isclose(A,-9999.0)]=np.nan
        else:
            A[np.isclose(A,-9999.0)]=0.0

        # Set lower boundary as zero. 
        A = np.concatenate((np.zeros((A.shape[0],1,A.shape[2],A.shape[3])),A),axis=1)
        # Apply cyclic boundaries if applicable
        if args.cyclicx:
            A = np.insert(A,0,A[:,:,:,-1],axis=3)
            A = np.insert(A,A.shape[3],A[:,:,:,1],axis=3)
        if args.cyclicy:
            A = np.insert(A,0,A[:,:,-1,:],axis=2)
            A = np.insert(A,A.shape[2],A[:,:,1,:],axis=2)
        z = np.insert(ds['z'][:].data,0,0.0)
        muoto = A[:,1:-1,1:-1,1:-1].shape
        if args.ddx:
            print('  d'+i+'dx')
            outddx['d'+i+'dx'] = ((A[:,1:-1,1:-1,2:] - A[:,1:-1,1:-1,:-2])
                                  / dx)
        if args.ddy:
            print('  d'+i+'dy')
            outddx['d'+i+'dy'] = ((A[:,1:-1,2:,1:-1] - A[:,1:-1,:-2,1:-1]) 
                                  / dy )
        if args.ddz:
            print('  d'+i+'dz')
            # In z direction, the mesh spacing can be nonuniform.
            # Calculate ddz with (4.3.7) from Hirsch.
            outddx['d'+i+'dz'] = ((( A[:,2:,1:-1,1:-1]
                                     - A[:,1:-1,1:-1,1:-1])
                                   * np.broadcast_to(
                                       np.reshape(
                                           np.broadcast_to(
                                               np.reshape(
                                                   (z[1:-1] - z[0:-2])
                                                   / (z[2:] - z[1:-1]),
                                                   (muoto[1],1)),
                                               (muoto[1],muoto[2])),
                                           (muoto[1],muoto[2],1)),
                                       muoto)
                                   + ( A[:,1:-1,1:-1,1:-1]
                                       - A[:,0:-2,1:-1,1:-1])
                                   * np.broadcast_to(
                                       np.reshape(
                                           np.broadcast_to(
                                               np.reshape(                                       
                                                   (z[2:] - z[1:-1])
                                                   / (z[1:-1] - z[0:-2]),
                                                (muoto[1],1)),
                                               (muoto[1],muoto[2])),
                                           (muoto[1],muoto[2],1)),
                                       muoto))
                                  / np.broadcast_to(
                                      np.reshape(
                                          np.broadcast_to(
                                              np.reshape(
                                                  z[2:] - z[:-2],
                                                  (muoto[1],1)),
                                              (muoto[1],muoto[2])),
                                          (muoto[1],muoto[2],1)),
                                      muoto))
             
        A = None
    else:
        print('* Warning: Variable '+i+' does not have four coordinates.')
        print('           Derivatives will not be calculated.')

    
#=output======================================================================#

dso = netcdfOutputDataset( args.fileout )

tv = createNetcdfVariable( 
    dso, ds['time'][:].data, 'time', len(ds['time'][:].data), uD['time'],
    'f4', ('time',), True )

xv = createNetcdfVariable( 
    dso, x, 'x' , len(x), uD['x'], 'f4', ('x',), True )

yv = createNetcdfVariable( 
    dso, y, 'y' , len(y), uD['y'], 'f4', ('y',), True )

zv = createNetcdfVariable( 
    dso, z[1:-1], 'z' , len(z[1:-1]), uD['z'],
    'f4', ('z',), True )

for i in outddx:
    iorig = i[1:-2]
    if '/' in uD[iorig]:
        nu = uD[iorig]+'m'
    else:
        nu = uD[iorig]+'/m'
        
    createNetcdfVariable(
        dso, outddx[i] , i , None , nu,
        'f4', ('time', 'z', 'y', 'x', ), False )

netcdfWriteAndClose( dso )
        

#                                     |_|/                                    #
#===================================|_|_|\====================================#
#                                     |                                       #
