#!/usr/bin/env python3
import numpy as np
import sys
import argparse
from utilities import filesFromList, writeLog

#==========================================================#
parser = argparse.ArgumentParser(prog='cumulateColumnsFromFiles.py')
parser.add_argument("strKey", help="Search string for collecting files.",\
  nargs='?', default=None)
parser.add_argument('-x', '--prefix', help='Prefix for file names. Default=CML_ ',\
  type=str, default='CX_')
parser.add_argument('-j', '--jcols', type=int, nargs='+',\
  help='Columns to cumulate with data format [x y1 y2 ... yn]. Ex: 1 2 4')
parser.add_argument('-jx', '--xcols', type=int, nargs='+', default=[0],\
  help='Columns to include without cumulations. Default=0 ')
parser.add_argument('-nr', '--nrows', type=int, default=None,\
  help='Number of rows in the cumulative plot array.')
args = parser.parse_args()
writeLog( parser, args )
#==========================================================# 
# Rename and cast for convenience ...
strKey = args.strKey
prefix = args.prefix
jcols  = args.jcols
xcols  = args.xcols
nr     = args.nrows

#print(' args = {} '.format(args))

allcols= list(); allcols.extend(xcols); allcols.extend(jcols)
print('All columns to extract: {}'.format(allcols))

if( not  strKey ): 
  strKey = input(" Enter search string: ")
  if( not strKey ): sys.exit(1)

# Gather the desired files:
fileNos, fileList = filesFromList( strKey+"*" )

# Process the files:
for fn in fileNos:
  
  d = np.loadtxt( fileList[fn], usecols=allcols )
  
  ddims = np.shape( d )
  
  nr  = min( nr , ddims[0] )
  ncols = len( allcols )
  
  cdims = (nr, ncols )
  cdat = np.zeros( cdims )
  
  idelta = np.round(ddims[0]/nr, decimals=0).astype(int)
  
  njx = 0
  for n,j in enumerate(xcols):
    njx += 1    
    for i in range(nr+1):
      irow = i*idelta; irow = min( irow, ddims[0]-1 )
      i = min(i, nr-1)
      cdat[i,n] = d[irow,j]
  
  for n,j in enumerate(jcols):
    for i in range(nr+1):
      irow = i*idelta; irow = min( irow, ddims[0] )
      i = min(i, nr-1)
      cdat[i,n+njx] = np.sum( d[:irow,j] )
  
  d = None
  fileout = prefix+fileList[fn]
  print(' Writing out file: {} '.format( fileout ) )
  np.savetxt(fileout, cdat[:,allcols], fmt='%3.6e')
  cdat = None

print(' All done! ')
