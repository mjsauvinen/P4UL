#!/usr/bin/env python3
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
parser.add_argument("-o", "--output", type=str, default="ascii", help="Output format of the LAD array: npz array [npz], ASCII array \
                    [ascii] or a NetCDF4 file [nc].")
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
print("Plant canopy height index (pch_index): {}".format(pch_index+1))
lad_3d = np.zeros((int(rdims[1]/dpx[1]),int(rdims[0]/dpx[0]),pch_index+1))

# Reverse the y-axis because of the top-left origo in raster
R=R[::-1,:]

for i in range(rdims[1]):
  for j in range(rdims[0]):
    if (R[j,i]==0):
      continue
    plant = R[j,i]
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
    i=int(i/dpx[0])
    j=int(j/dpx[0])

    if(alpha==0 or beta==0):
      i_top = int(np.round(ub/dpx[2]))+1
      lad_3d[i,j,0:i_top] = lad
      continue
    plant_dist, zdims = pct.betaDistributionProfile(alpha,beta,lb,ub,dpx[2])
    plant_model = pct.constructTreeFromProfile(plant_dist, zdims, lad, radius, dpx[0])

    lad_3d[i-int(np.floor(radius/dpx[1])):i+int(np.floor(radius/dpx[1]))+1,j-int(np.floor(radius/dpx[0])):j+int(np.floor(radius/dpx[0]))+1,0:np.shape(plant_model)[2]]+=plant_model

if(args.output=="ascii"):
  fx = open( args.fileout, 'w')
  # The first row of the file is number of vertical canopy layers
  nPc=np.shape(lad_3d)
  print(nPc)
  fx.write(str(nPc[2])+"\n")
  # Loop through the verical canopy columns
  for x in range(nPc[0]):
    for y in range(nPc[1]):
      if (np.all(lad_3d[x,y,:]==0)):
        # There is no need to write empty columns
        continue
      # Convert everything into strings and write
      # datatype x y col(:)
      lineStr = str(1)+","+str(x)+","+str(y)+","+",".join(map("{:.3g}".format,lad_3d[x,y,:]))+"\n"
      fx.write(lineStr)
  fx.close()
elif(args.output=="nc"):
  #Save into netCDF4 dataset
  dso = nct.netcdfOutputDataset(args.fileout)
  nPc = np.shape(lad_3d)
  xv = nct.createCoordinateAxis(dso, nPc, dpx, 0, 'x', 'f4', 'm', parameter=True)
  yv = nct.createCoordinateAxis(dso, nPc, dpx, 1, 'y', 'f4', 'm', parameter=True)
  zv = nct.createCoordinateAxis(dso, nPc, dpx, 2, 'z', 'f4', 'm', parameter=True)
  lad_3d=np.rollaxis(lad_3d,2)
  lad_3d=np.swapaxes(lad_3d,1,2)
  nc_ds = nct.createNetcdfVariable(dso, lad_3d, "leaf_area_density", 0, 'm', 'f4', ('z', 'y', 'x'), parameter=False)
  nct.netcdfWriteAndClose(dso)
elif(args.output=="npz"):
  np.savez(args.fileout,R=lad_3d,dPx=dpx,GlobOrig=ds["GlobOrig"])
else:
  raise ValueError("Unknown output format: {}".format(args.output))
