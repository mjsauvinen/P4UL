#!/usr/bin/env python3 
import numpy as np
import sys


try:
  Umag = float(sys.argv[1])
except:
  Umag = input(' Velocity magnitude (m/s): ') 
  
n = 32
deg = np.arange(0,n,2)
alpha = deg*(np.pi/180.)

Ux = Umag*np.cos(alpha)
Uy = Umag*np.sin(alpha)

for i in range(len(deg)):
  print '%d deg :  %f , %f '%(deg[i], Ux[i], Uy[i])
  
print ' ==================== '

for i in range(len(deg)):
  print '%d deg : liftDir (%f %f 0); dragDir (%f %f 0); '%(deg[i],-np.sin(alpha[i]),np.cos(alpha[i]) ,\
  np.cos(alpha[i]),np.sin(alpha[i]))
