#!/usr/bin/env python3

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

inFile = openIOFile( filename, 'r' )

for line in inFile:
  print(" Executing: {}".format(line))
  subprocess.run( line , shell=True )
