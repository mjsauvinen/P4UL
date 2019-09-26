#!/usr/bin/env python3

from math import *
from numpy import *

U    = input(" Velocity Scale (m/s): U = ")
rho = input(" Fluid density (kg/m^3): rho = ")
mu  = input(" Viscosity: mu = ")
L     = input(" Length Scale (m): L = ")
y_plus = input(" y+: y_plus = ")

Re  = rho * U * L / mu        # Reynolds number
Cf   = 0.455 / log( 0.06 * Re )**2 

ke   = 0.5 * rho * U**2       # Kinetic energy
tau_w = ke * Cf

u_tau = sqrt( tau_w / rho )
y1     =  y_plus * mu / ( rho * u_tau ) 


f = open('yDist.dat' ,'w')             #'w' = for writing
f.write(' Velocity Scale (m/s): U = %9.3e \n'  %  U  )
f.write(' Fluid density (kg/m^3): rho =  %9.3e \n'  %  rho )
f.write(' Viscosity: mu =  %9.3e \n'  %  mu  )
f.write(' Length Scale (m): L =  %9.3f \n'  %  L )
f.write(' Reynolds no.: Re =  %9.3f  \n'  %  Re )
f.write(' y+ =  %9.3e \n'  %  y_plus )
f.write(' First cell height (m): y1 =  %9.3e \n'  %  y1 )
f.close()


