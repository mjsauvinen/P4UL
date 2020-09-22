#!/usr/bin/env python3
import sys
import math
import numpy as np
from functools import reduce

# =*=*=*=* FUNCTION DEFINITIONS *=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def euler_rotation_matrix(theta_z=0, theta_y=0, theta_x=0):
  ''' 
  Generates a matrix for euler angle rotation:
  * First around z, then around y and finally around x axes.
  * Inputs angles are in radians.
  * Returns a (3,3) matrix.
  * Multiplication with vector shape (3,1) or (3,n)
  * The direction of rotation is given by the right-hand rule
  >> M  = euler_rotation_matrix(alpha, beta, gamma)
  >> v  = np.random.random( (3,1) )
  >> vt = np.dot(M, v)
  '''
  T = []
  if( theta_z ):
      cz = math.cos(theta_z)
      sz = math.sin(theta_z)
      T.append(np.array([[cz, -sz, 0.], [sz, cz, 0.], [0., 0., 1.]]))
  if( theta_y ):
      cy = math.cos(theta_y)
      sy = math.sin(theta_y)
      T.append(np.array([[cy, 0., sy], [0., 1., 0.], [-sy, 0., cy]]))
  if( theta_x ):
      cx = math.cos(theta_x)
      sx = math.sin(theta_x)
      T.append(np.array([[1., 0., 0.], [0., cx, -sx], [0., sx, cx]]))
  if T:
      return reduce(np.dot, T[::-1])
  return np.eye(3)

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def rotation_by_euler_angles( v , ea ):
  ''' 
  v: input array whose shape should be (3,n)
  ea: euler angles, whose shape should be (3,), that is, len(ea) == 3.    
  '''
  debug = True
  
  if( len(ea) != 3 ): 
    sys.exit(' Error! Wrong length of euler angles data. Exiting ...')
  if( np.shape(v)[0] != 3 ):
    sys.exit(' Error! Vector array is not of shape (3,n), but {}. Exiting ...'.format(np.shape(v)))
  
  # Obtain the transformation matrix T and perform matrix vector product.
  Tx = euler_rotation_matrix( ea[0], ea[1], ea[2] )
  if(debug): print(' T = {} '.format(Tx))
  
  return np.dot(Tx,v)

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def rotateRasterCoordNE( Rne , theta ):
  T = np.array([ [np.cos(theta), np.sin(theta)] , [-np.sin(theta), np.cos(theta)] ] )
  return np.matmul( T, Rne )
  
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


