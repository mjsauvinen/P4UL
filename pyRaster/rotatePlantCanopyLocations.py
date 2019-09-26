#!/usr/bin/env python3
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

from mapTools import saveTileAsNumpyZ

'''
Description:
Performs a non-interpolated rotation to individual plant canopy points in a given plant canopy raster map.
After the rotation performed by this script the extractDomainFromTile.py may be used with --noRotation argument.
Using rotation from extractDomainFromTile.py results in loss of individual trees due to grid discretization.

Author:
Sasu Karttunen <sasu.karttunen@helsinki.fi>
Institute for Atmospheric and Earth System Research (INAR) / Physics
Faculty of Science, University of Helsinki
'''

parser = argparse.ArgumentParser(prog='rotatePlantCanopyMap.py', \
                                 description='''Performs a non-interpolated rotation to individual plant canopy points in a given plant canopy raster map.''')
parser.add_argument("-f","--file", type=str, help="Name of the input plant raster map file.")
parser.add_argument("-fo","--fileout", type=str, help="Name of the output raster map file.")
parser.add_argument("-iP","--iPivot", help="Local pixel ids [N,E] for the pivot in the raster file.",\
  type=int,nargs=2,required=True)
parser.add_argument("-wd", "--windDir", type=float,default=270.0,\
  help="Wind direction (deg) --> Rotation angle around the pivot point. North wind = 0deg")
parser.add_argument("-p", "--printOn", help="Print the resulting raster data.",\
  action="store_true", default=False)
parser.add_argument("-pp", "--printOnly", help="Only print the resulting data. Don't save.",\
  action="store_true", default=False)
args = parser.parse_args()
# ------------------- #

ds=np.load(args.file)
R=ds["R"]
GlobOrig=ds["GlobOrig"]
dPx=[1.,1.]
pivot = np.array([[args.iPivot[0]],[args.iPivot[1]]])

# Calculate rotation in radians
theta = np.radians(270. - args.windDir)

# Construct a rotation matrix
c, s = np.cos(theta), np.sin(theta)
rotmat=np.array(((c,-s), (s, c))).T

# Create a new array for rotated points
rotatedMap = np.zeros(np.shape(R))

# Get list of unique ids, first one is 0 so let's skip that
ids = np.unique(R)[1:]

for id in ids:
  rotR = np.round(np.dot(rotmat, np.where(R==id)-pivot)).astype(int)
  rotR= rotR+pivot
  rotatedMap[rotR[0],rotR[1]] = id

# TODO: calculate new global origo
Rdict = {"R" : rotatedMap, "dPx" : dPx, "GlobOrig" : [0.,0.]}

if( not args.printOnly ):
  saveTileAsNumpyZ( args.fileout, Rdict)

if( args.printOn or args.printOnly ):
  plt.imshow(rotatedMap)
  plt.show()
# Calculate new global origo
# GlobOrig =
# np.savez("puutyyppi2c_rot.npz",R=arr,GlobOrig=[0.,0.],dPx=[1.,1.])
