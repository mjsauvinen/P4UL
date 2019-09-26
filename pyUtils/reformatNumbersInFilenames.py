#!/usr/bin/env python3
import sys
import argparse
import numpy as np
import subprocess as sb
from utilities import filesFromList, reformatNumberInWord
''' 
Description: 
        


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='reformatNumbersInFilenames.py')
parser.add_argument("strKey", help="Search string for collecting files.",default=None)
parser.add_argument("-sp", "--sepChar", type=str, default=".", \
  help="Separator character in file names. (default=\".\") ")
parser.add_argument("-p", "--printOn", help="Print the renaming commands on screen.",\
  action="store_true", default=False) 
parser.add_argument("-pp", "--printOnly", help="Print only the commands on screen.",\
  action="store_true", default=False) 
#==========================================================#

try: 
  args = parser.parse_args()
except:
  args = parser.parse_args([0])   # Send strKey=0
  args.strKey = raw_input(" Enter search string: ")
  if( not args.strKey ): sys.exit(1)
  
# Renaming for convenience ....
sepChar   = args.sepChar
printOn   = args.printOn
printOnly = args.printOnly

fileNos, fileList = filesFromList( "*"+args.strKey+"*" )

for filename in fileList:
  filename_new = reformatNumberInWord(filename, sepChar)
  cmd = 'mv {} {}'.format(filename, filename_new)
  if( printOn or printOnly ):
    print(cmd)
  if( not printOnly ):
    sb.call(cmd, shell=True)
  
  

    
      
