#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse
from plotTools import addToPlot
from netcdfTools import read3dDataFromNetCDF
from utilities import filesFromList
try:
  import bootstrapped.bootstrap as bs
  import bootstrapped.stats_functions as bs_stats
except:
  msg= ''' The bootstrapped package 
    (https://pypi.org/project/bootstrapped) 
    is not installed on your system. Exiting ...
    '''
  sys.exit(msg)

#==========================================================#
#==========================================================#
hStr = '''Indecies in the plane perpendicular to the profile plane. 
Ex. [i,j] for z-profile, [i,k] for y-profile, [j,k] for x-profile. 
If not given, the middle Indecies are chosen.'''

parser = argparse.ArgumentParser()
parser.add_argument("fileKey", default=None,\
  help="Search string for collecting files.")
parser.add_argument("-v", "--varname",  type=str, default='u',\
  help="Name of the variable in NETCDF file. Default='u' ")
parser.add_argument("-m", "--mode", type=str, default='mean', choices=['mean', 'std', 'var'],\
  help="Mode: mean, std or var")
parser.add_argument("-px", "--paxis", type=str, default='z', choices=['x', 'y','z'],\
  help="Axis along which profile is constructed. Choices = 'x','y', and 'z'. Default='z'")
parser.add_argument("-i", "--ij", type=int, nargs=2, default=[None,None], help=hStr)
parser.add_argument("-n", "--niter",  type=int, default=5000,\
  help="Number of iterations taken by the bootstrap algorithm. Defaut=5000")
parser.add_argument("-a", "--alpha",  type=float, default=0.05,\
  help="Value defining the confidence interval. Ex. alpha=0.05 refers to 95th-CI. Default=0.05")
parser.add_argument("-wa", "--writeAscii", action="store_true", default=False,\
  help="Save profile data to an ascii file.")
parser.add_argument("-s", "--save", type=str, default=None, \
  help="Name of the saved figure. Default=None")
args = parser.parse_args()    

#==========================================================#
# Rename ... that's all.
fileKey    = args.fileKey
varname    = args.varname
mode       = args.mode
alpha      = args.alpha
paxis      = args.paxis
ij         = args.ij
niter      = args.niter
writeAscii = args.writeAscii
saveFig    = args.save
#==========================================================#

if( mode == 'mean' ): sfunc = bs_stats.mean
if( mode == 'std'  ): sfunc = bs_stats.std
if( mode == 'var'  ): sfunc = bs_stats.var

# Obtain a list of files to include.
fileNos, fileList = filesFromList( fileKey+'*' )

fig = plt.figure(num=1, figsize=(12,10))
ax = fig.add_axes( [0.1, 0.075 , 0.875 , 0.81] ) #[left, top, width, height]

for fn in fileNos:
  dataDict = read3dDataFromNetCDF( fileList[fn] , varname, 1 )
  vr = dataDict['v']
  x  = dataDict['x']; y = dataDict['y']; z = dataDict['z']
  time = dataDict['time']
  dataDict = None
  
  nt, nk, nj, ni = vr.shape
  
  if( paxis == 'z' ):
    if( None in ij ): ij = [ni/2,nj/2]
    ij[0] = min( ij[0], ni-1 ); ij[1] = min( ij[1], nj-1 )
    vs = vr[:,:,ij[1],ij[0]]
    N  = nk
    d  = z
  elif( paxis == 'y' ):
    if( None in ij ): ij = [ni/2,nk/2]
    ij[0] = min( ij[0], ni-1 ); ij[1] = min( ij[1], nk-1 )
    vs = vr[:,ij[1],:,ij[0]]
    N  = nj
    d  = y
  else:
    if( None in ij ): ij = [nj/2,nk/2]
    ij[0] = min( ij[0], nj-1 ); ij[1] = min( ij[1], nk-1 )
    vs = vr[:,ij[1],ij[0],:]
    N  = ni
    d  = x
  
  vr     = None
  vres   = np.zeros( (N) )
  vlower = np.zeros( (N) )
  vupper = np.zeros( (N) )
  
  for i in xrange(N):
    res = bs.bootstrap(vs[:,i], alpha=alpha, stat_func=sfunc, num_iterations=niter, num_threads=2)
    print(' Step {} of {} done ...'.format(i,N))
    vres[i]   = res.value
    vlower[i] = res.lower_bound
    vupper[i] = res.upper_bound
    res = None
  
  nameList = fileList[fn].split('_')
  nstr = '{}{}'.format(nameList[0].replace('/MASK','_'), nameList[-1].strip('.nc'))
  ax.plot(vres, d, lw=2, label='{}, {}'.format(mode, nstr))
  ax.fill_betweenx(d, vupper, vlower, facecolor='gray', alpha=0.25)
  ax.set_title('{} of {}, {} '.format(mode, varname, nstr))
  
  if( writeAscii ):
    fname = '{}_{}_{}.dat'.format(varname, mode, nstr)
    print(' Writing data to ascii file: {}'.format(fname))
    hStr = '{}({}) [{}, value, lower, upper] '.format(mode,varname,paxis)
    np.savetxt(fname, np.c_[d, vres, vlower, vupper], header=hStr)

plt.legend(loc=0)
plt.grid(True)

if( saveFig ):
  fig.savefig( saveFig, format='jpg', dpi=300)

plt.show()

