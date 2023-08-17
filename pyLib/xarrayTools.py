#!/usr/bin/env python3
import xarray as xr

# A collection on helper functions that assume a dependency to xarray or are
# expected to process xarray objects.

# Jukka-Pekka Keskinen
# Finnish Meteorological Insitute
# 18.4.2023

#==========================================================================================#

def isCollocated(XRD,varstr):
    iscol = ( ('zu_3d' in XRD[varstr].coords)
              and ('x' in XRD[varstr].coords)
              and ('y' in XRD[varstr].coords) )
    return iscol

#                                         |_|/                                             #
#=======================================|_|_|\=============================================#
#                                         |                                                #
