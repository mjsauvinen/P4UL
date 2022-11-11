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
    description='Calculate SGS dissipation (ϵ).')
parser.add_argument('-f', '--filename',type=str, 
                    help='Name of the input data file. It has to contain '
                    'SGS TKE (e) and local SGS eddy diffusivity (K_m) in '
                    'their native grids .')
parser.add_argument('-fo', '--fileout',type=str, help='Name of output file.',
                    default = 'epsilon.nc')
#parser.add_argument('-n', '--missToNan',action="store_true", default=False,
#                    help='Set PALM missing values (-9999.0) to numpy in the'
#                    'calculation of Q. Default: set to 0.0.')

args = parser.parse_args()

#=inputs######================================================================#

ds, vD, uD = netcdfDataset2(args.filename)

for i in ['e', 'km']:
    if (i not in vD.keys() ):
        sys.exit(
            '{} not found from variable list: {}'.format(
                i, vD.keys()))

Cm = 0.1   # SGS model constant
        

x = 
y = ds['y'][:].data
z = ds['zu_3d'][:].data
t = ds['time'][:].data


#=calculations================================================================#

# Miten paikallinen hilaväli määritellään? Ei ole ihan selvää. Ehkä symmetrisesti?
dx = ds['x'][1:].data - ds['x'][:-1].data
dx = np.hstack((dx[0],dx))
dy = ds['y'][1:].data - ds['y'][:-1].data
dy = np.hstack((dy[0],dy))
dz = ds['zu_3d'][1:].data - ds['zu_3d'][:-1].data
dz = np.hstack((dz[0],dz))

epsi = 0.19*Cm*(ds['e'][:,:,:,:].data)**2+0.74*(ds['e'][:,:,:,:].data)**(3/2)


#=output======================================================================#

dso = netcdfOutputDataset( args.fileout )

tv = createNetcdfVariable( 
    dso, t, 'time', len(t), uD['time'], 'f4', ('time',), True )

xv = createNetcdfVariable( 
    dso, x, 'x' , len(x[1:-1]), uD['x'], 'f4', ('x',), True )

yv = createNetcdfVariable( 
    dso, y, 'y' , len(y[1:-1]), uD['y'], 'f4', ('y',), True )

zv = createNetcdfVariable( 
    dso, z, 'zu_3d' , len(z[1:-1]), uD['zu_3d'], 'f4', ('z',), True )

Ev = createNetcdfVariable( 
    dso, epsi , 'SGS dissipation' , None , 'm2/s2', 'f4', ('time', 'zu_3d', 'y', 'x', ), False )

netcdfWriteAndClose( dso )

#                                     |_|/                                    #
#===================================|_|_|\====================================#
#                                     |                                       #
