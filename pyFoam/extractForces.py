#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import numpy as np
import pylab as pl
from txtTools import openIOFile

# =*=*=*=* FUNCTION DEFINITIONS *=*=*=*=*=*=*=*=*=*=*=*

def isolateValues( line , stripChars ):
  v = []
  sl = line.split()
  for i in xrange(len(sl)):
    for sc in stripChars:
      sl[i] = sl[i].strip(sc)

  for s in sl:
    v.append(float(s))
    
  return v

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

try:
  factor = sys.argv[1]
except:
  factor = 1.

factor = float(factor)

f = openIOFile('forces.dat', 'r')
oc = openIOFile('forces.cmp', 'w')
ot = openIOFile('forces.tot', 'w')

lines = f.readlines()
spr = ['(',')']
Fx  = np.zeros(4,float)

for l in lines[1:]:
  x = np.array(isolateValues(l,spr)) 
  if( len(x) == 13 ):
    x.tofile(oc,sep=" \t"); oc.write("\n")

    Fx[0] = x[0]
    for i in xrange(1,len(Fx)):
      Fx[i]=factor*(x[i]+x[i+3]) # Pressure + Viscous
    Fx.tofile(ot, sep=" \t"); ot.write("\n")
    
f.close(); oc.close(); ot.close()