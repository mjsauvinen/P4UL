import numpy as np
import sys
'''
Description:
Functions for generating and manipulating plant canopy (trees etc.) data.

Author:
Sasu Karttunen <sasu.karttunen@helsinki.fi>
Institute for Atmospheric and Earth System Research (INAR) / Physics
Faculty of Science, University of Helsinki
'''

def betaDistributionProfile(alpha,beta,lb,ub,dz):
  from scipy.stats import beta as betadist
  """
  Return a beta probability density function with dimensions.

  Parameters
  ----------
    alpha (float) : alpha parameter
    beta (float) : beta parameter
    lb (float) : lowest nonzero level
    ub (float) : highest nonzero level
    dz (float) : resolution

  Returns
  -------
    dist (ndarray) : 1D array containing density function values
    z (ndarray) : 1D array containing dimensions
  """

  # Set the scale and loc so the x(lb) is the first
  # and x(ub) is the last nonzero level
  scale = ub-lb+2.*dz
  loc=lb-dz
  z = np.arange(0.0,ub+2.*dz,dz)
  dist = betadist.pdf(z,alpha,beta,loc=loc,scale=scale)

  return dist, z

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def constructTreeFromProfile(dist,z,lad,r,dx):
  '''
  Return a 3D array of a plant that corresponds to a given distribution.

  Parameters
  ----------
    dist (ndarray) : 1D array containing the vertical density profile
    z (ndarray) : 1D array containing the vertical profile dimensions
    lad (float) : constant leaf area density value to be used
    r (float) : maximum radius of the plant crown
    dx (float) : horizontal resultion of the grid

  Returns
  -------
    plant_3d (ndarray) : 3D array (2*r,2*r,len(z)-1) containing a 3D model of the plant
  '''
  n=int(2.*(r/dx))
  plant_3d = np.empty((n,n,len(z)))
  # Construct a profile for plant crown radius
  dist_sqrt = np.sqrt(dist)
  r_pr = dist_sqrt*(r/np.amax(dist_sqrt))

  # Create 2D slices and construct a 3D model of the plant
  hcenter = (n-dx)/2.
  y,x = np.ogrid[-hcenter:2.*(r/dx)-hcenter, -hcenter:2.*(r/dx)-hcenter]
  for k in np.arange(0,len(z)):
    # Calculate a circle mask with a radius of r_pr[z]
    r_k = r_pr[k]/dx
    cmask = x*x + y*y < r_k*r_k

    # Create a slice and use mask to fill it with LAD
    h_slice = np.zeros((n,n))
    h_slice[cmask] = lad
    plant_3d[:,:,k] = h_slice

  return plant_3d

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
