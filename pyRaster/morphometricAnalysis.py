#!/usr/bin/env python
import sys
import argparse
import numpy as np
from utilities import writeLog
from mapTools import *
'''
Description:
Estimate values for z_0 and z_d using morphometric methods described by Macdonald et al. (1998), Kanda et al. (2013) and Kent et al. (2017).


Author: Sasu Karttunen
        sasu.karttunen@helsinki.fi
        University of Helsinki
'''
#==========================================================#
parser = argparse.ArgumentParser(prog='morphometricAnalysis.py', description='''Estimate values for z_0 and z_d using morphometric methods described by Macdonald et al. (1998), Kanda et al. (2013) and Kent et al. (2017).''')
parser.add_argument("-bf","--buildings", type=str, help="Name of a raster map file of the buildings [npz].")
parser.add_argument("-wd","--windDir", type=float, default=0.0, help="Offset of wind direction from x-axis in degrees. Wind along x-axis is 0, crosswind is 90.")
parser.add_argument("-vf","--vegetation", type=str, help="Name of a 3-dimensional mask file of the vegetation [npz]. Used for Kent et al. (2017) method.")
parser.add_argument("--cdb", type=float, default=1.2, help="Drag coefficient for buildings.")
parser.add_argument("--p3d", type=float, default=0.2, help="Porosity of the vegetation.")
parser.add_argument("--alpha", type=float, default=4.43, help="Empirical constant for Macdonald et al (1998) method (used in all methods). The default is 4.43 (Hall et al. 1996).")
parser.add_argument("--beta", type=float, default=1.0, help="Empirical constant for Macdonald et al (1998) method (used in all methods). The default is 1.0 (Hall et al. 1996).")
args = parser.parse_args()
writeLog(parser, args)
#==========================================================#
# Define methods
def morphMacdonald(alpha, beta, planFrac, frontalFrac, h_av, cdb):
  # Morphometric method as described by Macdonald et al. (1998)
  mac_zd = (1.0+alpha**(-planFrac)*(planFrac-1.0))*h_av
  mac_z0 = ((1.0-(mac_zd/h_av))*np.exp(-(0.5*beta*(cdb/0.4**2)*(1.0-(mac_zd/h_av))*frontalFrac)**-0.5))*h_av
  return mac_z0, mac_zd

def morphKanda(h_av, h_max, std_h, planFrac, mac_z0):
  # Morphometric method as described by Kanda et al (2013)
  a0 = 1.29
  b0 = 0.36
  c0 = -0.17
  a1 = 0.71
  b1 = 20.21
  c1 = -0.77

  kanda_x = (std_h+h_av)/h_max
  kanda_zd = (c0*kanda_x**2.0+(a0*planFrac**b0-c0)*kanda_x)*h_av
  kanda_y = (planFrac*std_h)/h_av
  kanda_z0 = (b1*kanda_y**2.0+c1*kanda_y+a1)*mac_z0

  return  kanda_z0, kanda_zd

# Rename variables
alpha=args.alpha
beta=args.beta
cdb=args.cdb
p3d=args.p3d
windDir=args.windDir

# Read the building data.
bDict = readNumpyZTile(args.buildings)
bR = bDict['R']
bnPx = np.shape(bR)
bdPx = bDict['dPx']

# Calculate h_av, average roughness element height (2 m threshold to ignore park benches, cars and such)
h_av=np.mean(bR[np.where(bR>2.0)])
print("\nAverage height of the buildings: {:.2f} m".format(h_av))

# Maximum building height
h_max = np.amax(bR)
print("Maximum building height: {:.2f} m".format(h_max))

# Standard deviation of the building height
std_h = np.std(bR[np.where(bR>2.0)])
print("Standard deviation of the height of the buildings: {:.2f} m".format(std_h))

# Calculate plan area fraction planFrac=A_p/A_T
# Using unit areas and disregard resolution
planFracb=float(len(bR[np.where(bR>2.0)]))/np.prod(bnPx)
print("Plan area fraction of the builings: {:.4f}".format(planFracb))

# Frontal area fraction of buildings
frontalAreaXb, frontalAreaYb = frontalAreas(bR)
frontalFracXb = frontalAreaXb*np.cos(np.deg2rad(windDir))/np.prod(bnPx)
frontalFracYb = frontalAreaYb*np.sin(np.deg2rad(windDir))/np.prod(bnPx)
frontalFracb = frontalFracXb+frontalFracYb
print("Frontal area fraction of the buildlings: {:.4f}".format(frontalFracb))

mac_z0, mac_zd = morphMacdonald(alpha, beta, planFracb, frontalFracb, h_av, cdb)
print("\nMac z_0: {:.2f} m".format(mac_z0))
print("Mac z_d: {:.2f} m".format(mac_zd))

kanda_z0, kanda_zd = morphKanda(h_av,h_max,std_h,planFracb,mac_z0)
print("\nKanda z_0: {:.2f} m".format(kanda_z0))
print("Kanda z_d: {:.2f} m".format(kanda_zd))

# Inclusion of vegetation
# Plan area and frontal area fraction of vegetation
if (args.vegetation):
  # Inclusion of vegetation following Kent et al. (2017)
  print("\nIncluding vegetation:")
  vDict = readNumpyZTile(args.vegetation)
  vR = vDict['R']
  vnPx = np.shape(vR)
  vdPx = vDict['dPx']

  planFracvb = (float(len(bR[np.where(bR>2.0)]))+float(len(vR[np.where(vR>2.0)]))*(1.-p3d))/np.prod(bnPx)
  pv = (-1.251*p3d**2+0.489*p3d+0.803)/cdb
  frontalAreaXv, frontalAreaYv = frontalAreas(vR)
  frontalFracXv = frontalAreaXv*np.cos(np.deg2rad(windDir))/np.prod(vnPx)
  frontalFracYv = frontalAreaYv*np.sin(np.deg2rad(windDir))/np.prod(vnPx)
  frontalFracvb = frontalFracb+(frontalFracXv+frontalFracYv)*pv

  mac_z0, mac_zd = morphMacdonald(alpha, beta, planFracvb, frontalFracvb, h_av, cdb)
  kanda_z0, kanda_zd = morphKanda(h_av,h_max,std_h,planFracvb,mac_z0)

  print("\nMac z_0: {:.2f} m".format(mac_z0))
  print("Mac z_d: {:.2f} m".format(mac_zd))

  print("\nKanda z_0: {:.2f} m".format(kanda_z0))
  print("Kanda z_d: {:.2f} m".format(kanda_zd))
