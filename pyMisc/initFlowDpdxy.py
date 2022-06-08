#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import argparse

#==========================================================#
parser = argparse.ArgumentParser(prog='flow_init_dpdxy.py')
parser.add_argument("-wd", "--winddir", type=float, default=270.,
  help="Wind direction in degrees relative to grid: 0 equals north (-y),"
                    " 180 south (+y), 270 west (+x). Default = 270.")
parser.add_argument("-s", "--scale", type=float, default=1.,\
  help="Scaling factor [s] for the wind values: Uout = s*Uin.")
parser.add_argument("-f", "--filename", type=str, default="umag_profil.dat",
                    help="File with reference wind profile.")
args = parser.parse_args()
#==========================================================#

scale  = args.scale 
alpha  = (270. - args.winddir ) * np.pi/180.

# = = = = = = = = = = = = = = = 

z, Umag = np.loadtxt(args.filename, usecols=(0,1), unpack=True)
dpmag   = 0.0018  # Turku setup 

u = scale * Umag * np.cos( alpha )
v = scale * Umag * np.sin( alpha )

pdx = -dpmag * np.cos( alpha )
pdy = -dpmag * np.sin( alpha )


ustr = ' '.join('{:.1f},'.format(ui) for ui in u[::6])
vstr = ' '.join('{:.1f},'.format(vi) for vi in v[::6])
zstr = ' '.join('{:.1f},'.format(zi) for zi in z[::6])

u_output = 'u_profile  = {}'.format(ustr)
v_output = 'v_profile  = {}'.format(vstr)
z_output = 'uv_heights = {}'.format(zstr)

dp_output = 'dpdxy =  {:.8f}, {:.8f},'.format( pdx, pdy )

print('        !- wd = {} deg'.format(args.winddir))
print('        {}\n'.format(dp_output))
print('        {}'.format(u_output))
print('        {}'.format(v_output))
print('        {}'.format(z_output))
print('        ! - - - - - - - !\n')
