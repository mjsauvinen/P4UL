#!/usr/bin/env python3
import argparse
import sys
from netcdfTools import *

'''Calculate derivatives

Calculates derivatives for a given field. 

Author:
Jukka-Pekka Keskinen
Finnish Meteorological Insitute
12/2022

'''

#=functions===================================================================#

def one_sided_ddxy(ddxyarray, Mmm, Mm, M, Mp, Mpp, Amm, Am, A, Ap, App, d):
    # Calculate incorrectly calculated near-mask points using one-sided 
    # differenced ormulas.

    # Remove values calculated using masked points.
    ddxyarray[M] = -9999.0
    ddxyarray[Mm] = -9999.0
    ddxyarray[Mp] = -9999.0

    # Calculate new values with one-sided derivatives.                
    Mkn2 = 0
    Mkn1 = 0
    Mkn0 = 0

    # Point n+1 in mask, 2nd order backward difference
    Mk = ( Mp & ~M & ~Mm & ~Mmm )
    Mkn2 += np.count_nonzero(Mk[0,:,:,:])
    ddxyarray[Mk] = (3*A[Mk] - 4*Am[Mk] + Amm[Mk])/(2*d)
    Mk = None
                
    # Point n+1 in mask, 1st order backward difference
    Mk = ( Mp & ~M & ~Mm & Mmm )
    Mkn1 += np.count_nonzero(Mk[0,:,:,:])
    ddxyarray[Mk] = (A[Mk] - Am[Mk])/d
    Mk = None

    # Point n+1 in mask, not enough points for 1st order difference
    Mk = ( Mp & ~M & Mm)
    Mkn0 += np.count_nonzero(Mk[0,:,:,:])
    ddxyarray[Mk] = 0.0
    Mk = None
                
    # Point n-1 in mask, 2nd order forward difference
    Mk = ( Mm & ~M & ~Mp & ~Mpp )
    Mkn2 += np.count_nonzero(Mk[0,:,:,:])
    ddxyarray[Mk] = (-3*A[Mk] + 4*Ap[Mk] - App[Mk])/(2*d)
    Mk = None
                
    # Point n-1 in mask, 1st order forward difference
    Mk = ( Mm & ~M & ~Mp & Mpp )
    Mkn1 += np.count_nonzero(Mk[0,:,:,:])
    ddxyarray[Mk] = (Ap[Mk] - A[Mk])/d
    Mk = None

    # Point n-1 in mask, not enough points for 1st order difference
    Mk = ( Mm & ~M & Mp)
    Mkn0 += np.count_nonzero(Mk[0,:,:,:])
    ddxyarray[Mk] = 0.0
    Mk = None

    print('   Close to boundaries and topography, one-sided '
          'differences were used.')
    print('   2nd order: '+str(Mkn2)+' points ('
          + str(np.round(100*Mkn2/np.size(M),1))+' %)')
    print('   1st order: '+str(Mkn1)+' points ('
          + str(np.round(100*Mkn1/np.size(M),1))+' %)')
    print('   0th order: '+str(Mkn0)+' points ('
          + str(np.round(100*Mkn0/np.size(M),1))+' %)')

    return ddxyarray

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
                    help='Set PALM missing values (-9999.0) to nan in the'
                    'calculation of the derivatives. Default: set to 0.0.')
parser.add_argument('-cy', '--cyclicy', action="store_true", default=False,
                    help = 'Cyclic boundaries in y direction.')
parser.add_argument('-cx', '--cyclicx', action="store_true", default=False,
                    help = 'Cyclic boundaries in x direction.')
parser.add_argument('-a', '--all', action="store_true", default=False,
                    help = 'Calculate derivatives on all grid points ie. do '
                    'not skip boundaries. One-sided differences will be used '
                    'close to boundaries. Missing values will be set to ???.')
parser.add_argument('-zb', '--zerobottom', action="store_true", default=False,
                    help = 'Use zero for values below the grid.')
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

if args.cyclicx or args.all:
    x = ds['x'][:].data
else:
    x = ds['x'][1:-1].data

if args.cyclicy or args.all:
    y = ds['y'][:].data
else:
    y = ds['y'][1:-1].data

z = ds['z'][:].data

# Assume uniform mesh in x and y directions
dx = x[1]-x[0]
dy = y[1]-y[0]


for i in args.variable:
    if len(vD[i]) == 4:
        A = ds[i][:,:,:,:].data
        #        if args.missToNan or args.all:
        #            A[np.isclose(A,-9999.0)]=np.nan
        #        else:
        #            A[np.isclose(A,-9999.0)]=0.0

            
        if args.all:
            # True when missing value
            M = ds[i][:,:,:,:].mask
            if M.shape==():
                sys.exit('** Error: Variable '+i+' does not have a mask.')   
            #        elif args.zerobottom:
            # Set lower boundary as zero.
            #A = np.concatenate(
            #                (np.zeros((A.shape[0],1,A.shape[2],A.shape[3])),A),axis=1)                
            #z = np.insert(ds['z'][:].data,0,0.0)
        # Apply cyclic boundaries if applicable
            #        if args.cyclicx:
            #A = np.insert(A,0,A[:,:,:,-1],axis=3)
            #A = np.insert(A,A.shape[3],A[:,:,:,1],axis=3)
            #       if args.cyclicy:
            #A = np.insert(A,0,A[:,:,-1,:],axis=2)
            #A = np.insert(A,A.shape[2],A[:,:,1,:],axis=2)
       
        if args.ddx:
            print('  d'+i+'dx')

            # Initialise output with missing values.
            outddx['d'+i+'dx'] = -9999.0*np.ones(A.shape)
            
            # Calculates derivatives for inner points. This is incorrect close
            # to walls due to contamination by masked points.
            outddx['d'+i+'dx'][:,:,:,1:-1] = ((A[:,:,:,2:] - A[:,:,:,:-2])
                                              / 2*dx)
            
            if args.all:
                # Prepare difference masks
                temp = np.pad(M, ((0, 0), (0,0), (0,0), (2,2)),
                             'constant', constant_values=True)
                Mmm = temp[:,:,:,:-4] # n-2
                Mm = temp[:,:,:,1:-3] # n-1
                Mp = temp[:,:,:,3:-1] # n+1
                Mpp = temp[:,:,:,4:]  # n+2
                temp = None

                # Prepare Values
                temp = np.pad(A, ((0, 0), (0,0), (0,0), (2,2)),
                              'constant', constant_values=-9999.0)
                Amm = temp[:,:,:,:-4] # n-2
                Am = temp[:,:,:,1:-3] # n-1
                Ap = temp[:,:,:,3:-1] # n+1
                App = temp[:,:,:,4:]  # n+2
                temp = None

                
                outddx['d'+i+'dx'] = one_sided_ddxy(outddx['d'+i+'dx'],
                                                    Mmm, Mm, M, Mp, Mpp,
                                                    Amm, Am, A, Ap, App,
                                                    dx)

                Mmm = None
                Mm = None
                Mp = None
                Mpp = None
                Amm = None
                Am = None
                Ap = None
                App = None
                
        if args.ddy:
            print('  d'+i+'dy')
            # Initialise output with missing values.
            outddx['d'+i+'dy'] = -9999.0*np.ones(A.shape)
            
            # Calculates derivatives for inner points. This is incorrect close
            # to walls due to contamination by masked points.
            outddx['d'+i+'dy'][:,:,1:-1,:] = ((A[:,:,2:,:] - A[:,:,:-2,:]) 
                                  / 2*dy )
            
            if args.all:
                # Prepare difference masks
                temp = np.pad(M, ((0, 0), (0,0), (2,2), (0,0)),
                             'constant', constant_values=True)
                Mmm = temp[:,:,:-4,:] # n-2
                Mm = temp[:,:,1:-3,:] # n-1
                Mp = temp[:,:,3:-1,:] # n+1
                Mpp = temp[:,:,4:,:]  # n+2
                temp = None

                # Prepare Values
                temp = np.pad(A, ((0, 0), (0,0), (2,2), (0,0)),
                              'constant', constant_values=-9999.0)
                Amm = temp[:,:,:-4,:] # n-2
                Am = temp[:,:,1:-3,:] # n-1
                Ap = temp[:,:,3:-1,:] # n+1
                App = temp[:,:,4:,:]  # n+2
                temp = None

                outddx['d'+i+'dy'] = one_sided_ddxy(outddx['d'+i+'dy'],
                                                    Mmm, Mm, M, Mp, Mpp,
                                                    Amm, Am, A, Ap, App,
                                                    dy)

                Mmm = None
                Mm = None
                Mp = None
                Mpp = None
                Amm = None
                Am = None
                Ap = None
                App = None


        if args.ddz:
            print('  d'+i+'dz')
            # Initialise output with missing values.
            outddx['d'+i+'dz'] = -9999.0*np.ones(A.shape)
            
            muoto = A[:,1:-1,:,:].shape
            # In z direction, the mesh spacing can be nonuniform.
            # Calculate ddz with (4.3.7) from Hirsch.
            outddx['d'+i+'dz'][:,1:-1,:,:] = ((( A[:,2:,:,:]
                                     - A[:,1:-1,:,:])
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
                                   + ( A[:,1:-1,:,:]
                                       - A[:,0:-2,:,:])
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

            if args.all:
                
                # Prepare difference masks
                temp = np.pad(M, ((0, 0), (2,2), (0,0), (0,0)),
                             'constant', constant_values=True)
                Mmm = temp[:,:-4,:,:] # n-2
                Mm = temp[:,1:-3,:,:] # n-1
                Mp = temp[:,3:-1,:,:] # n+1
                Mpp = temp[:,4:,:,:]  # n+2
                temp = None

                # Prepare Values
                temp = np.pad(A, ((0, 0), (2,2), (0,0), (0,0)),
                              'constant', constant_values=-9999.0)
                Amm = temp[:,:-4,:,:] # n-2
                Am = temp[:,1:-3,:,:] # n-1
                Ap = temp[:,3:-1,:,:] # n+1
                App = temp[:,4:,:,:]  # n+2
                temp = None

                dz = z[1]-z[0]
                if np.any(~np.isclose(dz,z[1:]-z[0:-1])):
                    print('** Warning: Non-uniform mesh in z direction. '
                          'Accuracy will be degraded next ')
                    print('   to boundaries at following levels:')
                    for j in range(z.size):
                        if ~np.isclose(dz,z[j]):
                            print('   '+str(j))
                
                outddx['d'+i+'dz'] = one_sided_ddxy(outddx['d'+i+'dz'],
                                                    Mmm, Mm, M, Mp, Mpp,
                                                    Amm, Am, A, Ap, App,
                                                    dz)

                Mmm = None
                Mm = None
                Mp = None
                Mpp = None
                Amm = None
                Am = None
                Ap = None
                App = None

            
             
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
    dso, z, 'z' , len(z), uD['z'],
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
