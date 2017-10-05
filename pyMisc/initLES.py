#!/usr/bin/env python
import sys
import argparse
import numpy as np
from plotTools import addToPlot
import matplotlib.pyplot as plt

''' 
Description: 
        Rough order-of-magnitude analyzer for ABL flows.


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='initLES.py')
parser.add_argument("-u", "--umag",type=float, help="Velocity Scale [m/s].")
parser.add_argument("-l", "--lengthscale",type=float, help="Length Scale: Atmospheric Boundary Layer Height [m].")
parser.add_argument("-nu", "--viscosity",type=float, help="Kinematic Viscosity (mu/rho) [m^2/s^2].")
parser.add_argument("-z0", "--z0",type=float, help="(Optinal) Roughness length [m].", default=0.)
parser.add_argument("-d", "--dheight",type=float, help="(Optinal) Displacement height [m].", default=0.)
parser.add_argument("-us", "--ustar",type=float, help="(Optinal) Friction velocity.", default=1.)
parser.add_argument("-p", "--printOn", help="Print the velocity profile.",\
  action="store_true", default=False)
args = parser.parse_args() 
#writeLog( parser, args )
#==========================================================#

# Rename ... that's all
z0 = args.z0
d  = args.dheight
us = args.ustar
printOn = args.printOn


mnx = np.array([ 0.4 , 1.1])
mxn = np.array([ 1. , 1.])

u  = args.umag        # * mnx    
l  = args.lengthscale # * mxn
nu = args.viscosity   # * mnx

# Eddy turnover time: 
t = l/u

# Invicid energy containing eddies: Rate of viscous dissipation of kinetic energy. 
epsilon = u**3/l

# Kolmogorov velocity scale:
v = ( nu*epsilon )**(0.25)

# Kolmogorov length scale:
eta = (nu**3/epsilon)**(0.25)

# Turbulence Reynolds number.
Re_t = u*l/nu        

# Reynolds number of dissipative eddies --> should be unity.
Re_d = v*eta/nu

dStr = '''
# ================================================= #
# ===== LES order-of-magnitude data =============== # 
# ================================================= #

Velocity scale: 
  u  = {:.1e} m/s 

Length scale: ~Atmospheric boundary layer height.
  l  = {:.1e} m 

Kinematic Viscosity (mu/rho):
  nu = {:.1e} m^2/s^2

# ================================================= #

Eddy turnover time: 
  T = {:.1e} s 

Rate of viscous dissipation of kinetic energy: 
  epsilon = {:.1e} m^2/s^3

Kolmogorov velocity scale:  
  v = {:.1e} m/s

Kolmogorov length scale:    
  eta = {:.1e} m

Turbulence Reynolds number: 
  Re_t = {:.1e}

Re of dissipative eddies (should be unity):   
  Re_d = {:.1e}
'''

if( isinstance(u,np.ndarray) ):
  for i in xrange(len(u)):
    print(dStr.format(u[i],l[i],nu[i],t[i],epsilon[i],v[i],eta[i],Re_t[i],Re_d[i]))
elif( isinstance(u, np.float) ):
  print(dStr.format(u,l,nu,t,epsilon,v,eta,Re_t,Re_d))


# The log-law profile for the velocity
z = np.linspace(0.,l,int(l/8))
kappa = 0.41

# Um := mean(U)/u*
Um = us/kappa * np.log( np.maximum( (z-d), 1e-5)/z0 )

if( printOn ):
  ufig = plt.figure(1)
  ufig = addToPlot(ufig, Um,z,'<U>', ["Log-law Velocity profile","<U>","z"], logOn=False)
  plt.show()

dUm = Um[1:]-Um[:-1]
dz = z[1:]-z[:-1]
zm = (z[1:]+z[:-1])/2.
dUmdz = dUm/dz

print('u(z):\n'+' '.join('{:.2f},'.format(uv) for uv in Um[::2]))
print('z :\n'+' '.join('{},'.format(zi) for zi in z[::2].astype(int)))
print('v(z):\n'+' '.join('{},'.format(v) for v in np.zeros( len(z[::2]))))

print('dudz(z):\n'+' '.join('{:.3f},'.format(duz) for duz in dUmdz[::4]))
print('dudz(z):\n'+' '.join('{},'.format(zim) for zim in zm[::4].astype(int)))

