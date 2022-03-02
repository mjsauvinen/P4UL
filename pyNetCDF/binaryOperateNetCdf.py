#!/usr/bin/env python3
from netcdfTools import *
import sys
import argparse
import numpy as np

'''
Description:

'''
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
def checkVariables(vnames, vDict):
  # Check that variables are found
  for vi in vnames:
    if( vi not in vDict.keys() ):
      sys.exit(' Variable {} not found from variable list: {}'.format(vi, vDict.keys()))

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
def unaryOpr( v, opStr ):

  if( opStr is None ):
    return v

  if( opStr == '^2'):
    np.power(v, 2, out=v)
  elif( opStr == 'abs' ):
    np.abs(v, out=v)
  elif( opStr == 'sqrt'):
    np.sqrt(v, out=v)

  return v

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

def binaryOpr( v1, v2, opStr ):
  if(   opStr == '+' ):
    return v1+v2
  elif( opStr == '-' ):
    return v1-v2
  elif( opStr == '/' ):
    return v1/(v2+1.e-6)
  elif( opStr == '*' ):
    return v1*v2
  else:
    sys.exit('Unrecognized binary operator {}: Exiting ...'.format(opStr))
  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

def unaryOprUnit( unit1 , opStr ):

  if( opStr is None ):
    return unit1

  if( opStr == '^2'):
    unitout =  '({})^2'.format(unit1)
  elif( opStr == 'abs' ):
    unitout = unit1
  elif( opStr == 'sqrt'):
    unitout = '({})^(1/2)'.format(unit1)

  return unitout

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

def binaryOprUnit( unit1, unit2, opStr ):
  if(   opStr == '+' ):
    return unit1
  elif( opStr == '-' ):
    return unit1
  elif( opStr == '/' ):
    return '{}({})^-1'.format(unit1,unit2)
  elif( opStr == '*' ):
    if( unit1 == unit2 ): unitout = '({})^2'.format(unit1)
    else:                 unitout = '{} {}'.format(unit1,unit2)
    return unitout  


#==========================================================#
parser = argparse.ArgumentParser(prog='binaryOperateNetCdf.py')
parser.add_argument("-f1", "--filename1", metavar='FILE1', type=str,
  help="First NetCDF file.")
parser.add_argument("-f2", "--filename2", metavar='FILE2', type=str,
  help="Second NetCDF file.")
parser.add_argument("-fo", "--fileout", type=str, required=True, 
  help="Name of the output netCDF file.")
parser.add_argument("-vn1", "--varNames1", metavar='VN1', type=str, nargs='+',\
  help="Names of variables from f1 dataset.")
parser.add_argument("-vn2", "--varNames2", metavar='VN2', type=str, nargs='+',\
  help="Names of variables from f2 dataset.")
parser.add_argument("-vno", "--varOutNames", metavar='VNO', type=str, nargs='+',\
  help="Names of output variables. Their number must match with VN1.") 
parser.add_argument("-dn", "--derivNames",type=str, nargs='+', metavar='DN', default=None,\
  help="(Optional) Names of derived coordinates to output to same file.")
parser.add_argument("-s1", "--scale1", type=float, default=1.0, 
  help="Scale factor for file 1 dataset.") 
parser.add_argument("-s2", "--scale2", type=float, default=1.0, 
  help="Scale factor for file 2 dataset.")
parser.add_argument("-op", "--binaryOperator", type=str, choices=['+','-','/','*'],
  help="Binary operator: v1 <op> v2.") 
parser.add_argument("-uop1", "--unaryOperator1", type=str, choices=['^2','abs','sqrt'], default=None, 
  help="Unary operator for file 1 dataset.")
parser.add_argument("-uop2", "--unaryOperator2", type=str, choices=['^2','abs','sqrt'], default=None,
  help="Unary operator for file 2 dataset.")
args = parser.parse_args()
#==========================================================#
fn1 = args.filename1
fn2 = args.filename2
fileout = args.fileout
vn1     = args.varNames1
vn2     = args.varNames2
vno     = args.varOutNames
dn      = args.derivNames 
s1      = args.scale1
s2      = args.scale2
uop1    = args.unaryOperator1
uop2    = args.unaryOperator2
biop    = args.binaryOperator

'''
Establish two boolean variables which indicate whether the created variable is an
independent or dependent variable in function createNetcdfVariable().
'''
parameter = True; variable = False

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

N1 = len(vn1); N2 = len(vn2); NO = len(vno)
if( N1 == N2 ):
  pass
elif( (N1 > N2) and (N2 == 1) ):
  pass
else:
  sys.exit(' Incompatible number of variables: N1={} & N2={}. If N1 > N2, then N2 == 1 is required.'.format(N1,N2))

if( N1 != NO ):
  sys.exit(' The number of output variable names NO={} must match N1={}. Exiting ...'.format(NO,N1)) 

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

ds1, v1D, u1D = netcdfDataset2(fn1) # vD: variableDict, uD: unitDict
ds2, v2D, u2D = netcdfDataset2(fn2) 

checkVariables( vn1, v1D )
checkVariables( vn2, v2D )

vstr = vn1[0]  # We can use the first variable as the coords should match.
tn = v1D[vstr][0]; zn = v1D[vstr][1]; yn = v1D[vstr][2]; xn = v1D[vstr][3]


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# Create a NETCDF output dataset (dso) for writing out the data.
dso = netcdfOutputDataset( fileout )

# Read 
dD1 = dict()
dD2 = dict()

time, time_dims = read1DVariableFromDataset(tn, vstr , ds1, 0, 0, 1 ) # All values.
tv = createNetcdfVariable( dso, time, tn, len(time), u1D[tn],'f4', (tn,), parameter )

x, x_dims = read1DVariableFromDataset( xn, vstr, ds1, 0, 0, 1 )
xv = createNetcdfVariable( dso,  x  , xn , len(x)   , u1D[xn],'f4', (xn,), parameter )

y, y_dims = read1DVariableFromDataset( yn, vstr, ds1, 0, 0, 1 ) 
yv = createNetcdfVariable( dso,  y  , yn , len(y)   , u1D[yn],'f4', (yn,), parameter )

z, z_dims = read1DVariableFromDataset( zn, vstr, ds1, 0, 0, 1 )
zv = createNetcdfVariable( dso,  z  , zn , len(z)   , u1D[zn],'f4', (zn,), parameter )

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# Include additional (derived) coordinates into the output file.
if( dn ):
  for di in dn:
    if( di in v1D.keys() ):    
      dc = ds1.variables[di][:]
      uD = u1D[di]
      vD = v1D[di]
    elif( di in v2D.keys() ):
      dc = ds2.variables[di][:]
      uD = u2D[di]
      vD = v2D[di]
    else:
      sys.exit('Error: {} not found variable lists. Exiting ...'.format(di))
    
    dc_dims = np.shape( dc )
    dv = createNetcdfVariable( dso, dc, di, None, uD, 'f4', vD, variable )
    dc = None

time = None; x = None; y = None; z = None

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

for vi in vn1:
  vt , _ = read3DVariableFromDataset( vi, ds1, 0, 0, 0, 1 ) # All values.
  if( s1 != 1.0 ): vt *= s1
  vt = unaryOpr( vt, uop1 )
  dD1[vi] = vt
  
for vi in vn2:
  vt , _ = read3DVariableFromDataset( vi, ds2, 0, 0, 0, 1 ) # All values.
  if( s2 != 1.0 ): vt *= s2
  vt = unaryOpr( vt, uop2 )
  dD2[vi] = vt

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

for i in range( N1 ):
  if( N2 != 1 ):
    j = i
    vo = binaryOpr( dD1.pop(vn1[i]), dD2.pop(vn2[j]) , biop )
  else:          
    j = 0
    vo = binaryOpr( dD1.pop(vn1[i]), dD2[vn2[0]]     , biop )
  
  unit12 = binaryOprUnit( u1D[vn1[i]] , u2D[vn2[0]] , biop )
  
  vv = createNetcdfVariable( dso, vo, vno[i], None, unit12, 'f4',(tn,zn,yn,xn,) , variable )
  vo = None


netcdfWriteAndClose(dso)
