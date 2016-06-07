#!/usr/bin/env python

import os, sys
import argparse
import subprocess
from txtTools import openIOFile

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
parser = argparse.ArgumentParser(prog='executeCommandsInFile.py')
parser.add_argument("-f", "--filename", help="Name of file containing commands to be run.", type=str)
args = parser.parse_args() 
# = = = = = = = = = Main Program = = = = = = = = = = = = = #
filename = args.filename


try:
  inFile = openIOFile( filename, 'r' )
except:
  sys.exit("Cannot open file {}. Exiting ...".format(filename))


for line in inFile:
  print " Executing: "+line
  subprocess.call( line , shell=True )
