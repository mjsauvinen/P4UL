#!/usr/bin/env python
import sys
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-e","--fext", help="No. of spatial dimensions.", default="jpg")
args = parser.parse_args()    

searchStr = "*."+args.fext

cmd = "\\includegraphics[width=0.8\\textwidth]{{./Pics/{}}}"



fileList = []
files = glob.glob(searchStr)        # obtain the list of files 
fileList.extend(files)              # Lists are iterable
fileList.sort()                     # Sort the list alphabetically

for f in fileList:
  print cmd.format(f)