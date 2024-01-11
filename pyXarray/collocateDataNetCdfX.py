#!/usr/bin/env python3
import sys
import argparse
import numpy as np
import xarray as xr
import time
from xarrayTools import *

# A script to collocate PALM output data to the scalar grid. Utilises xarray.

# Jukka-Pekka Keskinen
# Finnish Meteorological Insitute
# 2023

#==========================================================================================#
start = time.time()

parser = argparse.ArgumentParser(prog='collocateDataNetCdfX.py', description='Collocates '
                                 'staggered variables to collocated coordinates.'
                                 'Coordinate locations are preserved.')
parser.add_argument('-f', '--filename',type=str,
                    help='Name of the input netCDF file.')
parser.add_argument('-o', '--fileout',type=str,
                    help='Name of the output netCDF file.')
parser.add_argument('-v', '--variables', type=str, nargs='*',                    
                    help='Names variables to collocate. Default: all.')
parser.add_argument('-i', '--interpolation',type=str, help='Interpolation method. Supports'
                    ' those supported by scipy, e.g. linear, nearest, zero, slinear, '
                    'quadratic, cubic, etc. Default: linear.', default='linear')
parser.add_argument('-e', '--noExtrapolation',help='Do not extrapolate. Grid points that'
                    ' can not be interpolated will be set to Nan', default=False,
                    action='store_true')
parser.add_argument('-z', '--zeroBoundaries',help='Interpolate using 0.0 at obstacles.',
                    default=False, action='store_true')
parser.add_argument('-m', '--maskVariable',help='Use topography mask from a specific '
                    'variable. It makes sense to use a scalar variable here. Very useful '
                    'with --zeroBoundaries.', type=str)
parser.add_argument('-r', '--rename',help='Rename zu_3d axis to z.',
                    default=False, action='store_true')
parser.add_argument('-s', '--removeSides',help='Remove the furthest layers from x and y '
                    'directions. These will be subject to extrapolation.', default=False,
                    action='store_true')


#==========================================================================================#

args = parser.parse_args()

if args.noExtrapolation:
    extrapolation = np.nan
else:
    extrapolation = 'extrapolate'


    
with xr.open_dataset(args.filename) as F:
    for i in ['x', 'y', 'zu_3d']:
        if i not in F.coords:
            sys.exit('*** Collocated coordinates not available in input file. '
                     'Terminating.')
    if args.variables == None:
        varis = F.data_vars
    else:
        varis = args.variables
    for i in varis:
        if ((i in F.data_vars) and not isCollocated(F,i)):
            if args.zeroBoundaries:
                F[i] = F[i].fillna(0.0)            
            print(' Collocating '+i)
            # Assuming here that only one of the cooordinates needs to be
            # changed.
            if 'zw_3d' in F[i].coords:
                F[i] = F[i].interp(zw_3d=F['zu_3d'], method=args.interpolation,
                                   kwargs={'fill_value': extrapolation})
            elif 'xu' in F[i].coords:
                F[i] = F[i].interp(xu=F['x'], method=args.interpolation,
                                   kwargs={'fill_value': extrapolation})
            elif 'yv' in F[i].coords:
                F[i] = F[i].interp(yv=F['y'], method=args.interpolation,
                                   kwargs={'fill_value': extrapolation})            
            if ((args.maskVariable != None) and isCollocated(F,i) ):
                # Get mask from a given variable (makes sense if it's a scalar)
                # and apply it to the newly collocated variable. This should
                # make sure that all output variables have the same topography
                # mask.
                # print(' Enforcing original mask from '+args.maskVariable+' to '+i+'.')
                F[i] = F[i].where(~np.isnan(F[args.maskVariable].data))
        else:
            print(' Skipping '+i+'.')

    if args.rename:
        # Conform with P4UL convention.
        F= F.rename({'zu_3d':'z'})

    if args.removeSides:
        # The last layers in both x and y directions will be extrapolated. This
        # causes problems, such as negative variances, some times. To avoid
        # problematic values, the extrapolated regions can be removed.
        F = F.isel(x=slice(0,-1),y=slice(0,-1))

        
    F.to_netcdf(args.fileout)
            
print(' Script finished after '+str(int((time.time()-start)/60))+' min.')

#                                         |_|/                                             #
#=======================================|_|_|\=============================================#
#                                         |                                                #
