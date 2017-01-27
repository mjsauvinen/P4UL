#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse
from utilities import filesFromList
from plotTools import userLabels, plotBar
''' 
Description:


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''

#==========================================================#
parser = argparse.ArgumentParser()
parser.add_argument("strKey", help="Search string for collecting files.",\
    nargs='?', default=None)
parser.add_argument("-r", "--rho", help="Density of fluid (kg m^-3).", type=float,\
    default=0.653)
parser.add_argument("--debug", help="Activate debug switch.", action="store_true",\
    default=False)
parser.add_argument("--finn", help="Tekstit suomeksi.", action="store_true",\
    default=False)
parser.add_argument("--legend", help="User Legend.", type=str,\
    default="Percentage Deviation from Mean Flow Rate") 
args = parser.parse_args()
#==========================================================#

if( not args.strKey ): 
  args.strKey = raw_input(" Enter search string: ")
  if( not args.strKey ): sys.exit(1)

rho = args.rho  # Rename
  
while 1:

  fileNos, fileList = filesFromList( "*"+args.strKey+"*" )
  
  fn = fileNos[0]
  dat = np.loadtxt(fileList[fn])
  xc    = dat[:,0]
  Umean = dat[:,1]
  width = dat[:,4]
  A = np.pi * width**2/4.
  mdot = rho * Umean * A
  mtot = np.sum(mdot)
  mdotAvg = np.mean(mdot)
  m_ratios = (mdot-mdotAvg)/mtot * 100.
  error = 0.0 * m_ratios  # This is an estimate, based on paraview data.
  
  
  if( args.finn ):
    titleStr = "Massavirran jakautuminen. \nPoikkeamat keskiarvosta (%) vs. X-koord."
    xlabelStr = "X-koord. [m]"; ylabelStr = " Poikkeamat [%] " 
  else:
    titleStr = "Mass Flow Distribution. \nDeviations (%) from Mean Value vs. X-coord."  
    xlabelStr = "X-coord. [m]"; ylabelStr = " Deviation [%] " 

  
  if( args.debug) : titleStr += " M_tot = {:g} kg/s".format(mtot)
  plotStr = [ titleStr, xlabelStr, ylabelStr]

  fig = plt.figure(num=1, figsize=(9.5,8.))
  legendStr = args.legend
  plotBar(fig, xc, m_ratios, legendStr, plotStr, width, error)
  
  #pfig = plt.figure(num=1, figsize=(7.,8.));
  #for fn in fileNos:
  #  pfig = plotXX( pfig, fileList[fn], args.log )

  #if(args.labels):
  #  pfig = userLabels( pfig )
  #plt.grid(True)
  plt.legend(loc=0)

  plt.show()