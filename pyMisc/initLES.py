#!/usr/bin/env python
import sys
import argparse
import numpy as np
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
args = parser.parse_args() 
#writeLog( parser, args )
#==========================================================#

mnx = np.array([ 0.4 , 1.1])
mxn = np.array([ 1. , 1.])

u  = mnx*args.umag 
l  = mxn*args.lengthscale
nu = mnx*args.viscosity 

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

for i in xrange(len(u)):
  print dStr.format(u[i],l[i],nu[i],t[i],epsilon[i],v[i],eta[i],Re_t[i],Re_d[i])



