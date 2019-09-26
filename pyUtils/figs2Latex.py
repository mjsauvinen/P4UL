#!/usr/bin/env python3
import sys
import glob
import argparse
''' 
Description: A script to quickly generate LaTex commands for adding figures to a document.

Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''

#==========================================================# 
parser = argparse.ArgumentParser()
parser.add_argument("-e","--fext", type=str, default="jpg",\
  help="Figure file extension (used for searching figures). Default=jpg")
args = parser.parse_args()    
#==========================================================# 
searchStr = "*"+args.fext

cmd = "\\includegraphics[width=0.8\\textwidth]{{./Pics/{}}}"

fileList = []
files = glob.glob(searchStr)        # obtain the list of files 
fileList.extend(files)              # Lists are iterable
fileList.sort()                     # Sort the list alphabetically

for f in fileList:
  print( cmd.format(f) )