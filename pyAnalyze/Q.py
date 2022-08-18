#!/usr/bin/env python3
import argparse
import sys
from netcdfTools import *

'''Calculate Q values
  
Calculate the Q values for the identification of vortex cores. The script
expects collocated velocity data. This can be achieved using
collocateDataNetCdf.py

Author:
Jukka-Pekka Keskinen
Finnish Meteorological Insitute
8/2022
'''

#=argument parser=============================================================#

parser = argparse.ArgumentParser(
    prog='Q.py',
    description='Calculate the Q values for the identification of vortices.')
parser.add_argument("-f", "--filename",type=str, 
                    help="Name of the input data file. It has to contain "
                    "all velocity components.")
parser.add_argument("-fo", "--fileout",type=str, 
                    help="Name of output file.")

args = parser.parse_args()

#=calculations================================================================#

ds, vD, uD = netcdfDataset2(args.filename)

for i in ['u', 'v', 'w']:
    if (i not in vD.keys() ):
        sys.exit(
            ' Vector component {} not found from variable list: {}'.format(
                i, vD.keys()))

muoto = ds['u'][:,1:-1,1:-1,1:-1].shape


dudx = ((ds['u'][:,1:-1,1:-1,2:].data - ds['u'][:,1:-1,1:-1,:-2].data) 
        / np.broadcast_to(ds['x'][2:].data - ds['x'][:-2].data, muoto))
dvdx = ((ds['v'][:,1:-1,1:-1,2:].data - ds['v'][:,1:-1,1:-1,:-2].data) 
        / np.broadcast_to(ds['x'][2:].data - ds['x'][:-2].data, muoto))
dwdx = ((ds['w'][:,1:-1,1:-1,2:].data - ds['w'][:,1:-1,1:-1,:-2].data) 
        / np.broadcast_to(ds['x'][2:].data - ds['x'][:-2].data, muoto))
dudy = ((ds['u'][:,1:-1,2:,1:-1].data - ds['u'][:,1:-1,:-2,1:-1].data) 
        / np.broadcast_to(np.reshape(
            ds['y'][2:].data - ds['y'][:-2].data, (muoto[2],1)), muoto))
dvdy = ((ds['v'][:,1:-1,2:,1:-1].data - ds['v'][:,1:-1,:-2,1:-1].data) 
        / np.broadcast_to(np.reshape(
            ds['y'][2:].data - ds['y'][:-2].data, (muoto[2],1)), muoto))
dwdy = ((ds['w'][:,1:-1,2:,1:-1].data - ds['w'][:,1:-1,:-2,1:-1].data) 
        / np.broadcast_to(np.reshape(
            ds['y'][2:].data - ds['y'][:-2].data, (muoto[2],1)), muoto))
dudz = ((ds['u'][:,2:,1:-1,1:-1].data - ds['u'][:,:-2,1:-1,1:-1].data) 
        / np.broadcast_to(np.reshape(np.broadcast_to(np.reshape(
            ds['z'][2:].data - ds['z'][:-2].data, (muoto[1],1)), 
          (muoto[1],muoto[2])), (muoto[1],muoto[2],1)), muoto))
dvdz = ((ds['v'][:,2:,1:-1,1:-1].data - ds['v'][:,:-2,1:-1,1:-1].data) 
        / np.broadcast_to(np.reshape(np.broadcast_to(np.reshape(
            ds['z'][2:].data - ds['z'][:-2].data, (muoto[1],1)), 
        (muoto[1],muoto[2])), (muoto[1],muoto[2],1)), muoto))
dwdz = ((ds['w'][:,2:,1:-1,1:-1].data - ds['w'][:,:-2,1:-1,1:-1].data) 
        / np.broadcast_to(np.reshape(np.broadcast_to(np.reshape(
            ds['z'][2:].data - ds['z'][:-2].data, (muoto[1],1)), 
        (muoto[1],muoto[2])), (muoto[1],muoto[2],1)), muoto))

Q = dudx**2 + dvdy**2 + dwdz**2 + 2*dudy*dvdx + 2*dudz*dwdx + 2*dvdz*dwdy

#=output======================================================================#

dso = netcdfOutputDataset( args.fileout )

tv = createNetcdfVariable( 
    dso, ds['time'][:], 'time', len(ds['time'][:]), uD['time'], 'f4', 
    ('time',), True )

xv = createNetcdfVariable( 
    dso, ds['x'][1:-1] , 'x' , len(ds['x'][1:-1]), uD['x'], 'f4', 
    ('x',), True )

yv = createNetcdfVariable( 
    dso, ds['x'][1:-1] , 'y' , len(ds['x'][1:-1]), uD['y'], 'f4', 
    ('y',), True )

zv = createNetcdfVariable( 
    dso, ds['z'][1:-1] , 'z' , len(ds['z'][1:-1]), uD['z'], 'f4', 
    ('z',), True )

Qv = createNetcdfVariable( 
    dso, Q , 'Q' , None , 'm2/s2', 'f4', ('time', 'z', 'y', 'x', ), False )

netcdfWriteAndClose( dso )
netcdfWriteAndClose( ds )
