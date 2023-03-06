#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import argparse
from scipy.interpolate import interp1d

#==========================================================#
parser = argparse.ArgumentParser(prog='flow_init_dpdxy.py',
                                 description="Generate PARIN entries "
                                 "for initial wind profiles and pressure "
                                 "gradient for a given wind direction. "
                                 "Velocity profile can be given as an "
                                 "external file which is then interpolated "
                                 "for desired intervals. The extent and "
                                 "magnitude of the profile is set using "
                                 "maximum height and maximum velocity "
                                 "magnitude. NB! Extrapolation of profiles "
                                 " can happen. Beware!")
parser.add_argument("-wd", "--winddir", type=float, default=270.,
  help="Wind direction in degrees relative to grid: 0 equals north (-y),"
                    " 180 south (+y), 270 west (+x). Default = 270.")
parser.add_argument("-um", "--Umax", type=float, default=None,\
  help="Maximum value for wind magnitude in the profile. Used for scaling "
"the magnitude of the wind profile.")
parser.add_argument("-zm", "--zmax", type=float, default=None,\
  help="Maximum height of the wind profile. Used for scaling the extent "
"of the wind profile.")
parser.add_argument("-f", "--filename", type=str, default=None,
                    help="File with reference wind profile.")
parser.add_argument("-dp", "--pressuregrad", type=float, default=0.002,
                    help="Magnitude of the pressure gradient.")
parser.add_argument("-z", "--uvheights", type=float, nargs='+',
                    help="Locations where velocity profiles are given. "
                    "Input a list of z coordinates separated by spaces. "
                    "Default: every 100 metres up to zmax.")
parser.add_argument("-i", "--interpolation", type=str, default="linear",
                    choices=["linear", "nearest", "nearest-up", "previous",
                             "next", "zero", "slinear", "quadratic", "cubic"],
                    nargs="?", help="Type of interpolation used for velocity "
                    "profiles.")
parser.add_argument("-o", "--outputPARIN", type=str, default=None,
                    help="Output the wind profile and pressure gradient "
                    "directly to an exsisting PARIN file. NB! Your old wind "
                    "profile and pressure gradient information will be "
                    "overwritten. Beware!")
args = parser.parse_args()
#==========================================================#

alpha   = (270. - args.winddir ) * np.pi/180.

# = = = = = = = = = = = = = = = 

if args.filename==None:
    z = np.arange(0,1000,100,dtype=float)
    Umag = np.array([0.0, 3.8, 5.8, 6.4, 6.9, 7.2, 7.4, 7.7, 7.9, 8.2])
else:
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

# Normalise and scale profiles. If already normalised, nothing should
# happen. If no scaling is specified, input profiles are used directly.
z = zscale*z/maxz
Umag = Uscale*Umag/maxUmag

dpmag   = args.pressuregrad

u = Umag * np.cos( alpha )
v = Umag * np.sin( alpha )

pdx = -dpmag * np.cos( alpha )
pdy = -dpmag * np.sin( alpha )


if args.uvheights==None:
    uvz = np.arange(0,zscale+0.1,100,dtype=float)
else:
    uvz = np.sort(np.asarray(args.uvheights,dtype=float))
    uvz = uvz[uvz>=0] # Remove negative heights

ui = interp1d(z,u,kind=args.interpolation,fill_value='extrapolate')
vi = interp1d(z,v,kind=args.interpolation,fill_value='extrapolate')

ustr = ' '.join('{:.1f},'.format(ui(zi)) for zi in uvz)
vstr = ' '.join('{:.1f},'.format(vi(zi)) for zi in uvz)
zstr = ' '.join('{:.1f},'.format(zi) for zi in uvz)

u_output = 'u_profile  = {}'.format(ustr)
v_output = 'v_profile  = {}'.format(vstr)
z_output = 'uv_heights = {}'.format(zstr)
dp_output = 'dpdxy =  {:.8f}, {:.8f},'.format( pdx, pdy )

if args.outputPARIN is not None:
    parinin  = open(args.outputPARIN,'r').readlines()
    parinout = open(args.outputPARIN,'w')
    wdoutlist = ['!- wd', 'dpdxy', 'u_profile', 'v_profile', 'uv_heights', 'wd -!']

    for j in parinin:
        if not any ( k in j for k in wdoutlist):
            parinout.write(j)
        if "&initialization_parameters" in j:
            parinout.write('\n        !- wd = {} deg'.format(args.winddir))
            parinout.write('\n        {}\n'.format(dp_output))
            parinout.write('\n        {}'.format(u_output))
            parinout.write('\n        {}'.format(v_output))
            parinout.write('\n        {}'.format(z_output))
            parinout.write('\n        ! - - - - - - wd -!\n')

    parinout.close()
    print("Wind profile and pressure gradient information written to "
          "file "+args.outputPARIN)
else:
    print('        !- wd = {} deg'.format(args.winddir))
    print('        {}\n'.format(dp_output))
    print('        {}'.format(u_output))
    print('        {}'.format(v_output))
    print('        {}'.format(z_output))
    print('        ! - - - - - - wd -!\n')
