#!/usr/bin/env python

from math import *


Uin = input(" Inlet Velocity (m/s): Uin = ")
iTu = input(" Turbulent Intensity (%): iTu = ")
lengthScale = input(" Length Scale (m): lenghtScale = ")

iTu = iTu / 100.0
Cmu = 0.09

k = 1.5*(iTu*Uin)**2
epsilon = (Cmu**0.75) * (k**1.5) / lengthScale
omega  = epsilon / ( Cmu * k + 1.0E-9 )

f = open('turbInit.dat' ,'w')             #'w' = for writing
f.write(' k_in = %9.3e \n'  %  k  )
f.write(' epsilon_in =  %9.3e \n'  %  epsilon )
f.write(' omega_in =  %9.3e \n'  %  omega )
f.write(' Inlet Velocity (m/s) =  %7.3f \n'  %  Uin )
f.write(' Turbulent Intensity =  %7.3f  \n'  %  iTu )
f.write(' Length Scale (m) =  %9.3e \n'  %  lengthScale )
f.close()


