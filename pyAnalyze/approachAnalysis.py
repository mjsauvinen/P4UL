#!/usr/bin/env python
import sys
import subprocess as sb
import numpy as np
import argparse
from utilities import filesFromList, writeLog
from plotTools import userLabels, extractFromCSV, addToPlot
import matplotlib.pyplot as plt
''' 
Description:


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''
#======== Function definitions =============================#

def p2pMaxMin( r ):
  # Peak to peak max and min evaluation routine.
  dr = (r[1:] - r[:-1])
  fpos = (dr>=0.).astype(int)
  fneg = (dr<0.).astype(int)

  rp_cum = 0.; rn_cum = 0.
  rp_max = 0.; rn_min = 0.
  i = 0
  for fp in fpos:
    if( fp == 0 ):
      if( rp_cum > rp_max ): rp_max = rp_cum
      rp_cum = 0.
    rp_cum += float(fp)*dr[i]; i+=1
    #print('rp_cum[{}] = {} '.format(i,rp_cum))

  i = 0
  for fn in fneg:
    if( fn == 0 ):
      if( rn_cum < rn_min ): rn_min = rn_cum
      rn_cum = 0.    
    rn_cum += float(fn)*dr[i]; i+=1
    #print('rn_cum[{}] = {} '.format(i,rn_cum))

  return rp_max, rn_min

#==========================================================#
parser = argparse.ArgumentParser(prog='approachAnalysis.py')
parser.add_argument("strKey", help="Search string for collecting files.",nargs='?',\
    default=".csv")
parser.add_argument("--magy", help="Magnitude of all variables.", action="store_true",\
    default=False)
parser.add_argument("--yx", help="Reverse axes: plot(x,y) --> plot(y,x)", action="store_true",\
    default=False)    
parser.add_argument("--labels", help="User specified labels.", action="store_true",\
    default=False)
parser.add_argument("--reuse", help="Reuse once specified variable selections.", action="store_true",\
    default=False)
parser.add_argument("-v", "--var", help="Variable Name in CSV-file", type=str, nargs='+',\
  default=['u','v','w'] )
parser.add_argument("-yl","--ylims", help="Y-axis limits: [min,max]. Default=[0,10]",\
  type=float,nargs=2,default=[0.,10.])
parser.add_argument("-fn","--figName", help="Name of the (temporary) figures. (default=tmp)",\
  type=str,default="tmp")
parser.add_argument("-fa","--fileAnim", help="Name of the animation file. (default=anim.gif)",\
  type=str,default="anim.gif")
parser.add_argument("-na", "--noAnim", help="Do not make an animation.",\
  action="store_true", default=False) 
args = parser.parse_args()    
writeLog( parser, args )
#==========================================================#

strKey      = args.strKey
figName     = args.figName
fileAnim    = args.fileAnim
noAnimation = args.noAnim
ylims       = args.ylims
varList     = args.var 

fileNos, fileList = filesFromList( "*"+strKey+"*" )

print(' The varList [-v, --var] option is over ridden at this point. ')
print(' Reading coordinate values from file {} ...'.format( fileList[0]) )
coordList = [ 'arc_length', 'Points:0', 'Points:1', 'Points:2']
xv = extractFromCSV( fileList[0]  , coordList )
s = xv[0].copy()  # arc_length
x = xv[1].copy(); y = xv[2].copy(); z = xv[3].copy()
xv = None
print(' Done.\n')

# -------------------------------------------------------- #
print(' Computing the mean velocity values ... ') 
varList = ['u', 'v', 'w']

Ux_mean = None; Uy_mean = None; Uz_mean = None

n = 0
for fn in fileNos:
  n += 1
  #pfig = pl.figure(num=1, figsize=(18.,9.))
  tv = extractFromCSV( fileList[fn]  , varList )
  u = tv[0].copy(); v = tv[1].copy(); w = tv[2].copy()
  tv = None
  
  # The velocity data may contain nan entries, which should be replaced by 0.
  u[np.isnan(u)] = 0.; v[np.isnan(v)] = 0.; w[np.isnan(w)] = 0.
  
  # Accumulate sums for mean values.
  if( Ux_mean == None ):
    Ux_mean = np.zeros( u.shape ) # Initialize
    Uy_mean = np.zeros( u.shape )
    Uz_mean = np.zeros( u.shape )
  Ux_mean += u; Uy_mean += v; Uz_mean += w
  
# Use the sums to compute the mean values.
Ux_mean /= float(n); Uy_mean /= float(n); Uz_mean /= float(n)
print(' Done.\n')

# -------------------------------------------------------- #

print(' Extract directional data from the approach line ... ')
#pfig = plotCSV( pfig, fileList[fn], args.yx, args.magy, args.reuse )

rad2deg = 180./np.pi
deg2rad = np.pi/180.

# Starting point: Rissala's approach line (from Paraview)
p1 = np.array([ x[0], y[0], z[0]  ])  # np.array([6800., 1250., 0.]) 
p2 = np.array([ x[-1],y[-1],z[-1] ])  # np.array([7700., 650., 72.]) 

da =  p2 - p1
da_mag = np.sqrt( np.sum( da**2 ) )
da_xy  = np.sqrt( np.sum( da[0:2]**2))

# Approach direction (normal vector)
na = da/da_mag

# Sharp angle between the runway and the mean wind
theta = np.arccos( da[0]/da_xy )
print(' Sharp angle between the runway and the mean wind: theta = {} deg'.format( theta*rad2deg ))
print(' Done.\n')


# -------------------------------------------------------- #
# Hornet's approach speed and velocity
Uappr_mag = 69. 
Ua = Uappr_mag*na


# Mean headwind
Uhw_mean = Ux_mean * np.cos( theta ) - Uy_mean * np.sin( theta )

# Speed relative to the ground ... perhaps not needed.
U_grd = Uappr_mag - Uhw_mean

# Approach angle
gamma = np.arctan( da[2]/da_xy )



# F18 Data:
rho = 1.2  # standard air
CL  = 1.2  # at 7deg angle of attack
CLa = 2.86 # 1/rad  (alpha in range [3deg, 10deg]) 
Aref=18.*3.
K = 0.5*rho*Aref

# Extract deviations in the headwind and compute the changes in AoA [alpha].
Lift = K*Uappr_mag**2*CL
n = 0
dL_max = 0.
dL_sum = 0.
dL_mxv = 0.   # Maximum variance.
dL_p2p_max = 0.
dL_p2p_min = 0.

for fn in fileNos:
  n += 1
  #pfig = pl.figure(num=1, figsize=(18.,9.))
  tv = extractFromCSV( fileList[fn]  , varList ) # NOTE: varList = ['u', 'v', 'w']
  du = tv[0]-Ux_mean
  dv = tv[1]-Uy_mean
  dw = tv[2]-Uz_mean   # Uz_mean could be replaced by 0.
  tv = None
  
  # The velocity data may contain nan entries, which should be replaced by 0.
  du[np.isnan(du)] = 0.; dv[np.isnan(dv)] = 0.; dw[np.isnan(dw)] = 0.
  
  dU_hw = du * np.cos( theta ) - dv * np.sin( theta )
  dalpha = np.arctan( dw/Uappr_mag)
  
  # Change in lift due to changes in AoA:
  dL_a = K*Uappr_mag**2*CLa*dalpha
  # Change in lift due to changes in head wind.
  dL_u = 2.*K*CL*Uappr_mag*dU_hw
  
  dLp_a = dL_a/Lift * 100.  # In percentage
  dLp_u = dL_u/Lift * 100.  
  dLp_mag= np.sqrt( (dLp_a+dLp_u)**2 )
  
  #fig = plt.figure(num=1, figsize=(18,9))
  fig, (ax1, ax2)  = plt.subplots(num=1, nrows=2, sharex=True, figsize=(18,11))
  lines11,=ax1.plot( s,dLp_a,'-o', linewidth=1.6 )
  lines12,=ax1.plot( s,dLp_u,'-o', linewidth=1.6 )
  ax1.legend( (lines11,lines12) , ('dL(alpha) [%]',' dL(u) [%]'), loc=1 )
  ax1.set_ylim([-8., 8.])
  ax1.set_xlim([ min(s) , 1.05*max(s)]) # s: arc_length
  ax1.set_title(' Changes in Lift due to turbulence ', fontsize=22)
  ax1.set_ylabel(' dL [%] ', fontsize=22); ax1.grid(True)
  
  
  lines2,=ax2.plot(s,dLp_mag,'-ro', linewidth=1.6 )
  ax2.legend( (lines2,) , (' ABS(SUM(dL)) [%]',), loc=1 )
  ax2.set_xlim([ min(s) , 1.05*max(s)]) # s: arc_length
  ax2.set_ylim([-1., 12.5]); ax2.set_xlim([ min(s) , max(s)])
  ax2.set_xlabel(' Distance along approach line [m] ', fontsize=22 )
  ax2.set_ylabel(' dL [%] ', fontsize=22 ); ax2.grid(True)
  
  # Maximum variance
  dL_ivar = np.var( dLp_mag[ du > 0 ] )  # Consider only nonzero values.
  if( dL_ivar > dL_mxv ): dL_mxv = dL_ivar
  
  # Mean variance
  dL_sum  += dL_ivar
  dL_var  = dL_sum/float(n)
  
  dL_imax = np.max(dLp_mag)
  if( dL_imax > dL_max ): dL_max = dL_imax
  
  dL_ip2p_mx, dL_ip2p_mn = p2pMaxMin( (dLp_a+dLp_u) )
  if( dL_ip2p_mx > dL_p2p_max ): dL_p2p_max = dL_ip2p_mx
  if( dL_ip2p_mn < dL_p2p_min ): dL_p2p_min = dL_ip2p_mn
  
  infoStr  =' Time = {:4d}s\n'.format((n-1)*2)
  infoStr +=' Current  P2P(dL) [max,min] = [{:4.1f}% , {:4.1f}%]\n'.format(dL_ip2p_mx, dL_ip2p_mn)
  infoStr +=' Running P2P(dL) [max,min] = [{:4.1f}% , {:4.1f}%]\n'.format(dL_p2p_max, dL_p2p_min)
  #infoStr +=' Max(dL) = {:4.1f}%\n'.format(dL_imax)
  infoStr +=' Running Max(dL) = {:4.1f}%\n'.format(dL_max)
  #infoStr +=' Var(dL) = {:4.1f}%\n'.format(dL_ivar)
  infoStr +=' Running Mean(Var(dL)) = {:4.1f}%\n'.format(dL_var)
  infoStr +=' Running Max(Var(dL))  = {:4.1f}%\n'.format(dL_mxv)
  
  
  plt.text( 1. , 5.5, infoStr , fontsize=20)
  figStr = '{}_{:04d}.jpg'.format(figName,n)
  print(' Saving figure {} '.format(figStr))
  fig.savefig(figStr)
  
  ax1.cla(); ax2.cla(); fig.clf()

if( not noAnimation ):
  cmd = 'convert {}_*  {} '.format(figName,fileAnim)
  print ' Executing command: ${}'.format(cmd)
  sb.call(cmd, shell=True)
  
print(' All Done! ')  
