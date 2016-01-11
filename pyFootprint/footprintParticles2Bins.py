#!/usr/bin/env python
from footprintTools import *
from mapTools import readNumpyZTile, entry2Int
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
''' 
Description:


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''
# = # = # = # Function definitions # = # = # = # = # = # = # 

def binBounds( s, Ns ):
  s_min = np.min(s)
  s_max = np.max(s)
  ds = (s_max-s_min)/float(Ns)
  s1 = np.zeros( Ns )
  s2 = np.zeros( Ns )
  for i in xrange(Ns):
    s1[i] = s_min + i*ds
    s2[i] = s1[i] + ds
  # Fix the last bound to encompass the whole range:
  s2[-1] = s_max
  
  return s1, s2, ds

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 

def sourceAreaMask( fmsk , xO , yO ):
  # Read the mask raster info.
  R, R_dims, ROrig, dPx = readNumpyZTile(fmsk)
  dPx = entry2Int( dPx )  # Resolution as a single number
  
  idM = (R[::-1,:]>0)  # True where the mask is nonzero.
  R = None     # Clear memory
  
  ixO = (xO/dPx).astype(int) # indecies of particle origins.
  iyO = (yO/dPx).astype(int) 
  
  iPrt = np.zeros( np.shape(ixO), bool) # Tmp array of particles, all False.
  iPrt[:] = idM[iyO[:],ixO[:]]  # Set True only the ones where particles exist.
  
  print(' From {} particles {} are contained within the mask.'.format(len(xO),np.sum(iPrt)))
  
  idM = None; ixO = None; iyO = None # Clear.
  
  return iPrt

# = # = # = # End Function definitions  # = # = # = # = # = #

#========================================================== #
parser = argparse.ArgumentParser(prog='footprintParticles2Bins.py')
parser.add_argument("-f", "--filename",type=str, help="Name of the raw footprint input (.npz) file.")
parser.add_argument("-fo", "--fileout",type=str, help="Name (Stem) for the output footprint files.",\
  default='FP')
parser.add_argument("-fm", "--filemask",type=str, help="Name of the mask (.npz) file.", default=None)
parser.add_argument("-N","--Nxyz", help="Number of dividing partitions [Nx, Ny, Nz] of the target volume.",\
  type=int, nargs=3, default=[2,2,2] )
args = parser.parse_args() 
#writeLog( parser, args )
#========================================================== #

# Rename ... that's all.
filename = args.filename
fileout  = args.fileout
filemask = args.filemask 
Nxyz    = args.Nxyz

# = = = = = = = = = = = = = = = = = = = = = = = = = = = =   #

xO, yO, zO,\
  xt, yt, zt, \
    ut, vt, wt = readNumpyZFootprintRaw( filename )


Nx = Nxyz[0]; Ny = Nxyz[1]; Nz = Nxyz[2]
x1, x2, dx = binBounds( xt, Nx ); print(' x1={}, x2={}], dx={}'.format(x1,x2,dx))
y1, y2, dy = binBounds( yt, Ny ); print(' y1={}, y2={}], dy={}'.format(y1,y2,dy))
z1, z2, dz = binBounds( zt, Nz ); print(' z1={}, z2={}], dz={}'.format(z1,z2,dz))

if( filemask ):
  imsk = sourceAreaMask( filemask , xO , yO )
else:
  imsk = True


for k in xrange( Nz ):
  kx = (zt>=z1[k]) * (zt<z2[k])
  for j in xrange( Ny ):
    jx = (yt>=y1[j]) * (yt<y2[j])
    for i in xrange( Nx ):
      ix =  (xt>=x1[i]) * (xt<x2[i])

      idx = ix*jx*kx*imsk  # Boolean matrix for chosen indecies.
      fstr = fileout+'_ijk_{0:02d}.{1:02d}.{2:02d}'.format(i,j,k)
      writeNumpyZFootprintIJK( fstr, xO[idx], yO[idx], zO[idx], \
        xt[idx], yt[idx], zt[idx], \
          ut[idx], vt[idx], wt[idx], np.array([dx,dy,dz]) )
