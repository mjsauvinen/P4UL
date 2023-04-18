#!/usr/bin/env python3
import sys
import argparse
import numpy as np
import xarray as xr
import time

# A script to collocate PALM output data to the scalar grid. Utilises xarray.

# Jukka-Pekka Keskinen
# Finnish Meteorological Insitute
# 2023

#=================================================================================#
start = time.time()

parser = argparse.ArgumentParser(prog='collocateDataNetCdfX.py')
parser.add_argument("-f", "--filename",type=str,
                    help="Name of the input netCDF file.")
parser.add_argument("-o", "--fileout",type=str,
                    help="Name of the output netCDF file.")
parser.add_argument("-v", "--variables", type=str, nargs='*',                    
                    help="Names variables to collocate. Default: all.")

#=================================================================================#

args = parser.parse_args()

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
        if i in F.data_vars:
            print(' Collocating '+i)
            # Assuming here that only one of the cooordinates needs to be
            # changed.
            if 'zw_3d' in F[i].coords:
                F[i] = F[i].interp(zw_3d=F['zu_3d'], kwargs={"fill_value": "extrapolate"})
            elif 'xu' in F[i].coords:
                F[i] = F[i].interp(xu=F['x'], kwargs={"fill_value": "extrapolate"})
            elif 'yv' in F[i].coords:
                F[i] = F[i].interp(yv=F['y'], kwargs={"fill_value": "extrapolate"})
        else:
            print(' Variable '+i+' not in input file. Skipping.')

    F.to_netcdf(args.fileout)
            
print(' Script finished after '+str(int((time.time()-start)/60))+' min.')

#                                     |_|/                                        #
#===================================|_|_|\========================================#
#                                     |                                           #  
