#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import argparse

#==========================================================#
parser = argparse.ArgumentParser(prog='flow_init_dpdxy.py')
parser.add_argument("-wd", "--winddir", type=float, default=270.,
  help="Wind direction in degrees relative to grid: 0 equals north (-y),"
                    " 180 south (+y), 270 west (+x). Default = 270.")
parser.add_argument("-U", "--Umax", type=float, default=None,\
  help="Maximum value for wind magnitude in the profile. Used for scaling "
"the magnitude of the wind profile.")
parser.add_argument("-z", "--zmax", type=float, default=None,\
  help="Maximum height of the wind profile. Used for scaling the extent "
"of the wind profile..")
parser.add_argument("-f", "--filename", type=str, default="umag_profil.dat",
                    help="File with reference wind profile.")
parser.add_argument("-dp", "--pressuregrad", type=float, default=0.002,
                    help="Magnitude of the pressure gradient.")
args = parser.parse_args()
#==========================================================#

alpha   = (270. - args.winddir ) * np.pi/180.
# = = = = = = = = = = = = = = = 

z, Umag = np.loadtxt(args.filename, usecols=(0,1), unpack=True)
maxUmag = np.max(Umag)
maxz    = np.max(z)

if args.Umax==None:
    Uscale = maxUmag
else:
    Uscale  = args.Umax
if args.zmax==None:
    zscale = maxz
else:
    zscale  = args.zmax 

# Normalise profiles. If already normalised, nothing happens
# If no scaling is specified, use input profiles directly.
z = z/maxz
Umag = Umag/maxUmag

dpmag   = args.pressuregrad

u = Uscale * Umag * np.cos( alpha )
v = Uscale * Umag * np.sin( alpha )

pdx = -dpmag * np.cos( alpha )
pdy = -dpmag * np.sin( alpha )


ustr = ' '.join('{:.1f},'.format(ui) for ui in u[::6])
vstr = ' '.join('{:.1f},'.format(vi) for vi in v[::6])
zstr = ' '.join('{:.1f},'.format(zscale*zi) for zi in z[::6])

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
