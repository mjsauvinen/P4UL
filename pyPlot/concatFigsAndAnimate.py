#!/usr/bin/env python3
import sys
import argparse
import numpy as np
import subprocess as sb
from utilities import filesFromList
''' 
Description: 
        


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='reformatNumbersInFilenames.py')
parser.add_argument("-f1", "--fig1", type=str, help="Search string for figures (1).")
parser.add_argument("-f2", "--fig2", type=str, help="Search string for figures (2).")
parser.add_argument("-p", "--printOn", help="Print the renaming commands on screen.",\
  action="store_true", default=False) 
parser.add_argument("-pp", "--printOnly", help="Print only the commands on screen.",\
  action="store_true", default=False) 
parser.add_argument("-na", "--noAnim", help="Do not make an animation.",\
  action="store_true", default=False) 
args = parser.parse_args()
#==========================================================#



# Renaming for convenience ....
printOn   = args.printOn
printOnly = args.printOnly
noAnimation = args.noAnim

fileNos1, fileList1 = filesFromList( "*"+args.fig1+"*" )
fileNos2, fileList2 = filesFromList( "*"+args.fig2+"*" )

nfigs = min( len(fileList1), len(fileList2) )

for i in range(nfigs):
  cmd = 'montage -mode concatenate -tile 1x2 {0} {1} output_{2:04d}.png'\
    .format(fileList1[i], fileList2[i], i )
  if( printOn or printOnly ):
    print(cmd)
  if( not printOnly ):
    sb.call(cmd, shell=True)


if( not noAnimation ):
  cmd = 'convert output_* anim.gif '
  if( printOn or printOnly ):
    print(cmd)
  if( not printOnly ):
    sb.call(cmd, shell=True)

