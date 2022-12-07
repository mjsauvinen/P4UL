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
                    'calculation of Q. Default: set to 0.0.')
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

for i in args.variable:
    if len(vD[i]) == 4:
        A = ds[i][:,:,:,:].data
        if args.missToNan:
            A[np.isclose(A,-9999.0)]=np.nan
        else:
            A[np.isclose(A,-9999.0)]=0.0

        muoto = A[:,1:-1,1:-1,1:-1].shape
        if args.ddx:
            print('  d'+i+'dx')
            outddx['d'+i+'dx'] = ((A[:,1:-1,1:-1,2:] - A[:,1:-1,1:-1,:-2])
                                  / np.broadcast_to(
                                      ds['x'][2:] - ds['x'][:-2],
                                      muoto))
        if args.ddy:
            print('  d'+i+'dy')
            outddx['d'+i+'dy'] = ((A[:,1:-1,2:,1:-1] - A[:,1:-1,:-2,1:-1]) 
                                  / np.broadcast_to(np.reshape(
                                      ds['y'][2:] - ds['y'][:-2],
                                      (muoto[2],1)),
                                                    muoto))
        if args.ddz:
            print('  d'+i+'dz')
            outddx['d'+i+'dz'] = ((A[:,2:,1:-1,1:-1] - A[:,:-2,1:-1,1:-1]) 
                                  / np.broadcast_to(
                                      np.reshape(
                                          np.broadcast_to(
                                              np.reshape(
                                                  ds['z'][2:] - ds['z'][:-2],
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
    dso, ds['x'][1:-1].data, 'x' , len(ds['x'][1:-1].data), uD['x'],
    'f4', ('x',), True )

yv = createNetcdfVariable( 
    dso, ds['y'][1:-1].data, 'y' , len(ds['y'][1:-1].data), uD['y'],
    'f4', ('y',), True )

zv = createNetcdfVariable( 
    dso, ds['z'][1:-1].data, 'z' , len(ds['z'][1:-1].data), uD['z'],
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
