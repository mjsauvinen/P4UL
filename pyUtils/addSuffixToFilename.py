#!/usr/bin/env python
import subprocess as sb
import argparse
import sys
from utilities import filesFromList, inputIfNone
'''
Description: Script for appending suffixes at the end of files. This is rarely necessary, but sometimes programs do not recognize the file types if they lack the correct identifier.


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi
        University of Helsinki &
        Finnish Meteorological Institute
'''


#==========================================================#
parser = argparse.ArgumentParser()
parser.add_argument("strKey", help="Search string for collecting files.",\
  nargs='?', default=None)
parser.add_argument("-s", "--suffix",type=str, help="File type (suffix) identifier.")
parser.add_argument("-pp", "--printOnly", help="Only print the commands. Don't execute.",\
  action="store_true", default=False)
args = parser.parse_args()
#==========================================================#
# Rename ...
strKey = args.strKey
suffix = args.suffix.strip('.')
printOnly = args.printOnly

strKey = inputIfNone( strKey , " Enter search string: " )

fileNos, fileList = filesFromList( "*"+strKey+"*" )

for fn in fileNos:
  s0 = fileList[fn]
  s1 = s0.split('.')[0]
  s2 = s1+'.'+suffix
  cmd = 'mv {}  {} '.format(s0,s2)
  print(cmd)
  if( not printOnly ): 
    try: sb.call(cmd, shell=True)
    except: pass


print('Done!')