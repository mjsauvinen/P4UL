#!/usr/bin/env python3
import argparse
import sys
from netcdfTools import *

'''Calculate Q values
  
Calculate the Q values (Hunt, Wray & Moin, 1988) for the identification of
vortex cores. The script expects collocated velocity data. This can be achieved
using collocateDataNetCdf.py

NB: Output dataset has different lenght in z, y, and x directions. Also, there
might be some space for improvement at the boundaries.

Author:
Jukka-Pekka Keskinen
Finnish Meteorological Insitute
8/2022

'''

#=argument parser=============================================================#

parser = argparse.ArgumentParser(
    prog='Q.py',
    description='Calculate the Q values for the identification of vortices.')
parser.add_argument('-f', '--filename',type=str, 
                    help='Name of the input data file. It has to contain '
                    'all velocity components.')
parser.add_argument('-fo', '--fileout',type=str, help='Name of output file.',
                    default = 'Q.nc')
parser.add_argument('-n', '--missToNan',action="store_true", default=False,
                    help='Set PALM missing values (-9999.0) to numpy in the'
                    'calculation of Q. Default: set to 0.0.')

args = parser.parse_args()

#=inputs######================================================================#

ds, vD, uD = netcdfDataset2(args.filename)

for i in ['u', 'v', 'w']:
    if (i not in vD.keys() ):
        sys.exit(
            ' Vector component {} not found from variable list: {}'.format(
                i, vD.keys()))

if (vD['u']!=vD['v'] or vD['u']!=vD['w']):
    sys.exit(' Velocity components are not in the same locations. Exiting.') 

muoto = ds['u'][:,1:-1,1:-1,1:-1].shape

x = ds['x'][:].data
y = ds['y'][:].data
z = ds['z'][:].data
t = ds['time'][:].data


#=calculations================================================================#

# In order to consume less memory, variables are emptied when they are not
# needed anymore.

print('  Calculating Q.')

u = ds['u'][:,:,:,:].data
if args.missToNan:
    u[np.isclose(u,-9999.0)]=np.nan
else:
    u[np.isclose(u,-9999.0)]=0.0

dudx = ((u[:,1:-1,1:-1,2:] - u[:,1:-1,1:-1,:-2])
        / np.broadcast_to(x[2:] - x[:-2], muoto))

Q = dudx**2
dudx = None

dudy = ((u[:,1:-1,2:,1:-1] - u[:,1:-1,:-2,1:-1]) 
        / np.broadcast_to(np.reshape(
            y[2:] - y[:-2], (muoto[2],1)), muoto))
dudz = ((u[:,2:,1:-1,1:-1] - u[:,:-2,1:-1,1:-1]) 
        / np.broadcast_to(np.reshape(np.broadcast_to(np.reshape(
            z[2:] - z[:-2], (muoto[1],1)), 
          (muoto[1],muoto[2])), (muoto[1],muoto[2],1)), muoto))

u = None

v = ds['v'][:,:,:,:].data
if args.missToNan:
    v[np.isclose(v,-9999.0)]=np.nan
else:
    v[np.isclose(v,-9999.0)]=0.0

dvdx = ((v[:,1:-1,1:-1,2:] - v[:,1:-1,1:-1,:-2]) 
        / np.broadcast_to(x[2:] - x[:-2], muoto))

Q = Q + 2*dudy*dvdx
dudy = None
dvdx = None 

dvdy = ((v[:,1:-1,2:,1:-1] - v[:,1:-1,:-2,1:-1]) 
        / np.broadcast_to(np.reshape(
            y[2:] - y[:-2], (muoto[2],1)), muoto))

Q = Q + dvdy**2
dvdy = None

dvdz = ((v[:,2:,1:-1,1:-1] - v[:,:-2,1:-1,1:-1]) 
        / np.broadcast_to(np.reshape(np.broadcast_to(np.reshape(
            z[2:] - z[:-2], (muoto[1],1)), 
        (muoto[1],muoto[2])), (muoto[1],muoto[2],1)), muoto))

v = None

w = ds['w'][:,:,:,:].data
netcdfWriteAndClose( ds )

if args.missToNan:
    w[np.isclose(w,-9999.0)]=np.nan
else:
    w[np.isclose(w,-9999.0)]=0.0

dwdx = ((w[:,1:-1,1:-1,2:] - w[:,1:-1,1:-1,:-2]) 
        / np.broadcast_to(x[2:] - x[:-2], muoto))

Q = Q + 2*dudz*dwdx
dudz = None
dwdx = None

dwdy = ((w[:,1:-1,2:,1:-1] - w[:,1:-1,:-2,1:-1]) 
        / np.broadcast_to(np.reshape(
            y[2:] - y[:-2], (muoto[2],1)), muoto))

Q = Q + 2*dvdz*dwdy
dvdz = None
dwdy = None

dwdz = ((w[:,2:,1:-1,1:-1] - w[:,:-2,1:-1,1:-1]) 
        / np.broadcast_to(np.reshape(np.broadcast_to(np.reshape(
            z[2:] - z[:-2], (muoto[1],1)), 
        (muoto[1],muoto[2])), (muoto[1],muoto[2],1)), muoto))

w = None
Q = Q + dwdz**2 
dwdz = None




#=output======================================================================#

dso = netcdfOutputDataset( args.fileout )

tv = createNetcdfVariable( 
    dso, t, 'time', len(t), uD['time'], 'f4', ('time',), True )

xv = createNetcdfVariable( 
    dso, x[1:-1], 'x' , len(x[1:-1]), uD['x'], 'f4', ('x',), True )

yv = createNetcdfVariable( 
    dso, y[1:-1], 'y' , len(y[1:-1]), uD['y'], 'f4', ('y',), True )

zv = createNetcdfVariable( 
    dso, z[1:-1], 'z' , len(z[1:-1]), uD['z'], 'f4', ('z',), True )

Qv = createNetcdfVariable( 
    dso, Q , 'Q' , None , 'm2/s2', 'f4', ('time', 'z', 'y', 'x', ), False )

netcdfWriteAndClose( dso )

#                                     |_|/                                    #
#===================================|_|_|\====================================#
#                                     |                                       #
   
