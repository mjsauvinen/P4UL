#!/usr/bin/env python

import os, sys
import subprocess
from txtTools import openIOFile

# = = = = = = = = = Function definitions = = = = = = = = = #

# = = = = = = = = = Main Program = = = = = = = = = = = = = #

try:
	fileName =  sys.argv[1]
except:
	print "Usage:", sys.argv[0], " [ file which contains executable commands ] "
	sys.exit(1)

inFile = openIOFile( fileName, 'r' )

for line in inFile:
	print " Executing: "+line
	subprocess.call( line , shell=True )
