#!/usr/bin/env python3

from math import *
from numpy import *

Omega =  input(" Rotational Velocity (s^-1): Omega = ")
radius   = input(" Characteristic Radius (m): radius = ")
Ux       = input(" Axial Velocity Scale (m/s): Ux = ")
D         = input(" Characteristic Duct Diameter (m): D = ")
rho      = input(" Fluid density (kg/m^3): rho = ")
mu       = input(" Viscosity (kg m^-1 s^-1): mu = ")
y_plus  = input(" y+: y_plus = ")

Ut  = Omega * radius

Re_r  = rho * Ut * radius / mu        # Rotational Reynolds number
Re_x  = rho * Ux * D / mu             # Axial Duct Reynolds number
Cf_r   = 0.455 / log( 0.06 * Re_r )**2            # Flat plate correlation
Cf_x   = 0.455 / log( 0.06 * Re_x )**2

ke_r   = 0.5 * rho * Ut**2       # Kinetic energy
ke_x   = 0.5 * rho * Ux**2

tau_wr = ke_r * Cf_r
tau_wx = ke_x * Cf_x

u_tau_r = sqrt( tau_wr / rho )
u_tau_x = sqrt( tau_wx / rho )

y1_r     =  y_plus * mu / ( rho * u_tau_r ) 
y1_x     =  y_plus * mu / ( rho * u_tau_x ) 


f = open('yDist.dat' ,'w')             #'w' = for writing
f.write(' y+ =  %9.3e \n'  %  y_plus )
f.write(' Fluid density (kg/m^3): rho =  %9.3e \n'  %  rho )
f.write(' Viscosity (kg m^-1 s^-1): mu =  %9.3e \n'  %  mu  )
f.write(' Rotational Info: \n ')
f.write(' Rotational Velocity (s^-1): Omega = %9.3e \n'  %  Omega  )
f.write(' Tangential Velocity Scale (m/s): Ut = %9.3e \n'  %  Ut  )
f.write(' Radius (m): radius =  %9.3f \n'  %  radius )
f.write(' Rotational Reynolds no.: Re_r =  %9.3f  \n'  %  Re_r )
f.write(' First cell height (m): y1_r =  %9.3e \n'  %  y1_r )
f.write(' Axial Info: \n ')
f.write(' Axial Velocity Scale (m/s): Ux = %9.3e \n'  %  Ux  )
f.write(' Characteristic Duct Diameter (m) : D =  %9.3f  \n'  %  D )
f.write(' Axial Reynolds no.: Re_x =  %9.3f  \n'  %  Re_x )
f.write(' First cell height (m): y1_x =  %9.3e \n'  %  y1_x )
f.close()


