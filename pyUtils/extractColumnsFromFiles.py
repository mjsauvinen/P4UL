#!/usr/bin/env python3
import numpy as np
import sys
import argparse
from utilities import filesFromList, writeLog
from numTools import rotation_by_euler_angles
from mapTools import applyFilter

#==========================================================#
parser = argparse.ArgumentParser(prog='extractColumnsFromFiles.py')
parser.add_argument("strKey", help="Search string for collecting files.",\
  nargs='?', default=None)
parser.add_argument('-x', '--prefix', help='Prefix for file names. Default=CX_ ',\
  type=str, default='CX_')
parser.add_argument('-c', '--cols', type=int, help='Columns to extract. Ex: 0,1,2,4',\
  nargs='+')
parser.add_argument('-s', '--scale', help='Scaling factors for columns. Ex: 1,10,10,10',\
  nargs='+', type=float, default=[1.])
parser.add_argument('-vc', '--vcols', type=int, help='Vector columns for transformation (x,y,z).',\
  nargs=3, default=[None,None,None])
EulerHelp=''''Euler angles (around z-, y- and x-axes) in [deg] for applying rotation of 3D vector values, 
which are specified via -vc or --vcols arguments. 
Note! The rotation is applied in the following order: z, y, x. '''
parser.add_argument('-ea', '--euler', type=float, help=EulerHelp, nargs=3, default=[None,None,None])
parser.add_argument("-lf", "--localFilter", type=int, default=None,\
  help="Apply local filter with a given integer width to all columns.")

args = parser.parse_args()
writeLog( parser, args )
#==========================================================# 
# Rename and cast for convenience ...
strKey = args.strKey
prefix = args.prefix
Icols  = args.cols
vcols  = args.vcols
scale  = np.array(args.scale)  # Should be numpy array.
euler  = args.euler
fltN   = args.localFilter      # Local filter width 

if( fltN is None ):
  filterOn = False
else:
  filterOn = True

rotationOn = False
# Check whether vector rotation will be performed.
if( (vcols.count(None) == 0) and (euler.count(None) == 0) ):
  rotationOn = True
  eulerAngles = np.array(euler)*(np.pi/180.)
else:
  rotationOn = False

# Assemble Tcols -- target columns -- for vector rotation.
if( rotationOn ):
  Tcols = np.zeros( 3 , int )
  for i in range(3):
    vi = vcols[i]
    try: Tcols[i] = Icols.index(vi) # where the vector indecies are in Icols.
    except: sys.exit(' Error: Vector indices not in --cols. Exiting ...')


if( not  strKey ): 
  strKey = input(" Enter search string: ")
  if( not strKey ): sys.exit(1)

sc = np.ones( np.shape(Icols) , float )
if( len(scale) == 1 ):
  sc[:] = scale[0]
elif( len(scale) == len(Icols) ):
  sc[:] = scale[:]
else:
  sys.exit(' Error: incompatible array lengths for columns and scaling factors. Exiting.')
  

# Gather the desired files:
fileNos, fileList = filesFromList( strKey+"*" )

# Process the files:
for fn in fileNos:
  
  dat = np.loadtxt( fileList[fn], usecols=Icols )
  
  for j in Icols:
    if j not in range(dat.shape[1]):
      sys.exit(' Error: {} not in range({}). Exiting ...'.format(j,dat.shape[1]) )
  
  # Perform scaling
  dat[:,Icols] *= sc  
  
  # Carry out rotation
  if( rotationOn ):
    v0 = np.transpose( dat[:,Tcols] ) # The array must be tranposed to make its shape (3,n)
    v1 = rotation_by_euler_angles( v0 , eulerAngles )
    dat[:,Tcols] = np.transpose(v1); v0 = None # Copy the values back and clear memory.
    
  if( filterOn ):
    for j in Icols:
      dat[:,j] = applyFilter( dat[:,j], ['local', fltN ] )
  
  
  fileout = prefix+fileList[fn]
  print(' Writing out file: {} '.format( fileout ) )
  np.savetxt(fileout, dat[:,Icols], fmt='%3.6e')
  dat = None


print(' All done! ')

