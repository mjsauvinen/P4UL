#!/usr/bin/env python
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import csv

import plantCanopyTools as pct
import netcdfTools as nct


'''
Description:
Generates a 3D array of LAD values based on information about individual plants.

Author:
Sasu Karttunen <sasu.karttunen@helsinki.fi>
Institute for Atmospheric and Earth System Research (INAR) / Physics
Faculty of Science, University of Helsinki
'''

parser = argparse.ArgumentParser(prog='generatePlantCanopy.py', description='''Generates a 3D array of LAD values
                                 based on information about individual plants.''')
parser.add_argument("-pr","--plantmap", type=str, help="Name of the input plant raster map file.")
parser.add_argument("-pi","--plantinfo", type=str, help="Name of the input plant information file.")
parser.add_argument("-r", "--resolution", type=float, nargs=3, help="Resolution (x,y,z) of the output 3D array.")
parser.add_argument("-fo","--fileout", type=str, help="Name of the 3D output file.")
args = parser.parse_args()
# ------------------- #

dpx = args.resolution

# Read in the plant information
ds=np.load(args.plantmap)
R=ds["R"]
rdims=np.shape(R)
plant_info = {}
with open(args.plantinfo, 'r') as csvfile:
  reader = csv.reader(csvfile, delimiter=',')
  # Generate a dictionary with plant ids as keys and following list as a value:
  # 0: tree radius [m]
  # 1: lower bound of the crown [m]
  # 2: upper bound of the crown [m]
  # 3: leaf area density [m^2/m^3]
  # 4: alpha
  # 5: beta
  for row in reader:
    try:
      plant_id = int(row[0])
    except:
      continue
    if(plant_id in plant_info):
      print("Warning: plant id {} multiply defined in the plant information file.".format(plant_id))
    plant_data = np.array(map(lambda item: float(item), row[1:]))
    plant_info[plant_id] = plant_data

# Calculate grid point index for canopy upper boundary and initialize 3D LAD array
pch_index = int(np.ceil(np.amax(np.array(plant_info.values())[:,2])/dpx[2]))+1
print("Plant canopy height index (pch_index): {}".format(pch_index))
lad_3d = np.zeros((rdims[0],rdims[1],pch_index+1))

for i in xrange(rdims[0]):
  for j in xrange(rdims[1]):
    if (R[i,j]==0):
      continue
    plant = R[i,j]
    if plant not in plant_info:
      warnings.warn("Warning: plant id {} found from input raster but not defined in the plant info file. Skipping...".format(plant))
      continue

    datalist = plant_info[plant]
    radius = datalist[0]
    lb = datalist[1]
    ub = datalist[2]
    lad = datalist[3]
    alpha = datalist[4]
    beta = datalist[5]

    if(alpha==0 or beta==0):
      i_top = int(np.round(ub/dpx[2]))+3
      lad_3d[i,j,0:i_top] = lad
      continue

    plant_dist, zdims = pct.betaDistributionProfile(alpha,beta,lb,ub,dpx[2])
    plant_model = pct.constructTreeFromProfile(plant_dist, zdims, lad, radius, dpx[0])

    lad_3d[i-int(np.floor(radius)):i+int(np.floor(radius))+1,j-int(np.floor(radius)):j+int(np.floor(radius))+1,0:np.shape(plant_model)[2]]=plant_model


# Save into netCDF4 dataset
dso = nct.netcdfOutputDataset(args.fileout)
nPc = np.shape(lad_3d)
dPc = np.array([1.,1.,0.75])
xv = nct.createCoordinateAxis(dso, nPc, dPc, 0, 'x', 'f4', 'm', parameter=True)
yv = nct.createCoordinateAxis(dso, nPc, dPc, 1, 'y', 'f4', 'm', parameter=True)
zv = nct.createCoordinateAxis(dso, nPc, dPc, 2, 'z', 'f4', 'm', parameter=True)
lad_3d=np.rollaxis(lad_3d,2)
lad_3d=np.fliplr(lad_3d) # move the origo to left-bottom corner
nc_ds = nct.createNetcdfVariable(dso, lad_3d, "leaf_area_density", 0, 'm', 'f4', ('z', 'y', 'x'), parameter=False)
nct.netcdfWriteAndClose(dso)
