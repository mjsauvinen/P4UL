#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse
from plotTools import addToPlot
from netcdfTools import read3dDataFromNetCDF
from utilities import filesFromList, writeLog
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

parser = argparse.ArgumentParser(prog='bootstrapDataNetCdf.py')
parser.add_argument("fileKey", default=None,\
  help="Search string for collecting files.")
parser.add_argument("-v", "--varnames",  type=str, nargs='+', default=['u'],\
  help="Variable names in NETCDF file. Multiple entries means covariance. Default=['u']")
parser.add_argument("-vs", "--vstar",type=float, nargs='+', default=[1.],\
  help="Characteristic values v* (vs) used in (v+ =(v-v0)/v*). Default=[1].")
parser.add_argument("-v0", "--vref",type=float, nargs='+', default=[0.],\
  help="Reference value v0 (vref) used in (v+ =(v-v0)/v*). Default=[0].")
parser.add_argument("-xs", "--xscale",type=float, default=1.,\
  help="Coordinate scaling value (xs) used in (x+ =x/xs). Default=1.")
parser.add_argument("-m", "--mode", type=str, default='mean', choices=['mean', 'std', 'var'],\
  help="Mode: mean, std or var")
parser.add_argument("-px", "--paxis", type=str, default='z', choices=['x', 'y','z'],\
  help="Axis along which profile is constructed. Choices = 'x','y', and 'z'. Default='z'")
parser.add_argument("-i", "--ij", type=int, nargs=2, default=[None,None], help=hStr)
parser.add_argument("-n", "--niter",  type=int, default=5000,\
  help="Number of iterations taken by the bootstrap algorithm. Defaut=5000")
parser.add_argument("-a", "--alpha",  type=float, default=0.05,\
  help="Value defining the confidence interval. Ex. alpha=0.05 refers to 95th-CI. Default=0.05")
parser.add_argument("-fo", "--fileout",  type=str, default=None,\
  help="String for the output file name. Default: None (file is automatically named)")
parser.add_argument("-all", "--allfiles", help="Select all files automatically.",\
  action="store_true", default=False)
parser.add_argument("-wa", "--writeAscii", action="store_true", default=False,\
  help="Save profile data to an ascii file.")
parser.add_argument("-wao", "--writeAsciiOnly", action="store_true", default=False,\
  help="Save profile data to an ascii file.")
parser.add_argument("-s", "--save", type=str, default=None, \
  help="Name of the saved figure. Default=None")
args = parser.parse_args()    
writeLog( parser, args )
#==========================================================#
# Rename ... that's all.
fileKey    = args.fileKey
varnames   = args.varnames
v0         = np.array( args.vref  )  # Convert to numpy array
vs         = np.array( args.vstar )
xs         = args.xscale
mode       = args.mode
alpha      = args.alpha
paxis      = args.paxis
ij         = args.ij
niter      = args.niter
fileout    = args.fileout
saveFig    = args.save
allfiles   = args.allfiles
writeAscii     = args.writeAscii or args.writeAsciiOnly
writeAsciiOnly = args.writeAsciiOnly

#==========================================================#

if( mode == 'mean' ): sfunc = bs_stats.mean
if( mode == 'std'  ): sfunc = bs_stats.std
if( mode == 'var'  ): sfunc = bs_stats.var

# Obtain a list of files to include.
fileNos, fileList = filesFromList( fileKey+'*', allfiles )

fig = plt.figure(num=1, figsize=(12,10))
ax = fig.add_axes( [0.1, 0.075 , 0.875 , 0.81] ) #[left, top, width, height]

for fn in fileNos:
  vr = None 
  for i in xrange(len(varnames)):
    dataDict = read3dDataFromNetCDF( fileList[fn] , varnames[i], 1 )
    if( vr is None ):
      vr = dataDict['v']; vr -= v0[i]; vr /= vs[i]
      if( len(varnames) > 1): 
        # If we are computing Reynolds stresses, subtract the mean.
        vr -= np.mean( vr, axis=(0) )
    else:
      tr = dataDict['v']; tr -= v0[i]; tr /= vs[i]; tr -= np.mean( tr, axis=(0) )
      vr *= tr
      tr = None 
  
  
  x  = dataDict['x']; y = dataDict['y']; z = dataDict['z']
  time = dataDict['time']
  dataDict = None
  
  nt, nk, nj, ni = vr.shape
  
  if( paxis == 'z' ):
    if( None in ij ): ij = [ni/2,nj/2]
    ij[0] = min( ij[0], ni-1 ); ij[1] = min( ij[1], nj-1 )
    vb = vr[:,:,ij[1],ij[0]]
    N  = nk
    d  = z
  elif( paxis == 'y' ):
    if( None in ij ): ij = [ni/2,nk/2]
    ij[0] = min( ij[0], ni-1 ); ij[1] = min( ij[1], nk-1 )
    vb = vr[:,ij[1],:,ij[0]]
    N  = nj
    d  = y
  else:
    if( None in ij ): ij = [nj/2,nk/2]
    ij[0] = min( ij[0], nj-1 ); ij[1] = min( ij[1], nk-1 )
    vb = vr[:,ij[1],ij[0],:]
    N  = ni
    d  = x
  
  vr     = None
  vres   = np.zeros( (N) )
  vlower = np.zeros( (N) )
  vupper = np.zeros( (N) )
  d /= xs
  
  for i in xrange(N):
    res = bs.bootstrap(vb[:,i], alpha=alpha, stat_func=sfunc, num_iterations=niter, num_threads=2)
    print(' Step {} of {} done ...'.format(i,N))
    vres[i]   = res.value
    vlower[i] = res.lower_bound
    vupper[i] = res.upper_bound
    res = None
  
  if( fileout is None ):
    nameList = fileList[fn].split('_')
    nstr = '{}{}'.format(nameList[0].replace('/MASK','_'), nameList[-1].strip('.nc'))
  else:
    nameList = fileList[fn].split('UEX')
    nstr = 'UEX{}'.format(nameList[-1].strip('.nc'))

  
  if( writeAscii ):
    varnames = np.array( varnames )
    fname = '{}_{}_{}.dat'.format(varnames.tostring(), mode, nstr)
    print(' Writing data to ascii file: {}'.format(fname))
    hStr = '{}({}) [{}, value, lower, upper] '.format(mode,varnames.tostring(),paxis)
    np.savetxt(fname, np.c_[d, vres, vlower, vupper], header=hStr)

  if( not writeAsciiOnly ):
    ax.plot(vres, d, lw=2, label='{}, {}'.format(mode, nstr))
    ax.fill_betweenx(d, vupper, vlower, facecolor='gray', alpha=0.25)
    ax.set_title('{} of {}, {} '.format(mode, varnames.tostring(), nstr))

    plt.legend(loc=0); plt.grid(True)
    
    if( saveFig ):
      fig.savefig( saveFig, format='jpg', dpi=300)
    
    plt.show()

print(' Done! ')