#!/usr/bin/python3
import sys
import argparse
from utilities import filesFromList
from plotTools import userLabels, extractFromCSV, addToPlot
import matplotlib.pyplot as plt
''' 
A script for processing vertical mean profiles extracted from Palm LES analysis with Paraview.

Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='profilesFromCSV.py')
parser.add_argument("strKey", help="Search string for collecting files.",nargs='?',\
  default=".csv")
parser.add_argument("--yx", help="Reverse axes: plot(x,y) --> plot(y,x)", action="store_true",\
  default=False)
parser.add_argument("-y", "--yvars", help="Y-variables in CSV-file", nargs='+',\
  default=["um","vm","wm"])
parser.add_argument("-x", "--xvar", help="X-variable in CSV-file", type=str,\
  default="arc_length")
parser.add_argument("--labels", help="User specified labels.", action="store_true",\
  default=False)
parser.add_argument("-a", "--appendv", help="Variable in CSV-file to append data.", type=str,\
  default="Points:0")
parser.add_argument("-af", "--appendf", help="Multiplication factor for appending variable.", type=float,\
  default=0.0025)
args = parser.parse_args()
#==========================================================#

# Renaming ...
strKey  = args.strKey
varList = args.yvars
xvar    = args.xvar
avar    = args.appendv
af      = args.appendf

nY = len(args.yvars)

# Append the x- and a-variables into the variable list. 
varList.append(xvar)
varList.append(avar)

fileNos, fileList = filesFromList( "*"+strKey+"*" )

pfig = plt.figure(num=1, figsize=(12,10))

for fn in fileNos:
  dat = extractFromCSV( fileList[fn] , varList )
  a = dat[-1]; x = dat[-2]; y = dat[:-2]
  dat = None
  
  for i in range(nY):
    astr ='_'+str(int(a[0]/100.))
    y[i]+=(a*af)
    axLabels = ["{}(z)".format(varList[i]),"{} +x/100 [m/s]".format(varList[i]),"z [m]"]
    varLabel = varList[i]+astr
    pfig = addToPlot(pfig,y[i], x, varLabel, axLabels, False )
    
if(args.labels): pfig = userLabels( pfig )
  
plt.legend(loc=0)
plt.grid(True)
plt.show()
