#!/usr/bin/env python3

import netCDF4 as nc
import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse
import scipy.ndimage as sn # contains the filters
from plotTools import addImagePlot
from netcdfTools import read3dDataFromNetCDF
from utilities import selectFromList

#==========================================================#
def readVar( fn, vstr, cl=1 ):
  xDict = read3dDataFromNetCDF( fn , vstr , cl )
  v = xDict['v']; x = xDict['x']; y = xDict['y']; z = xDict['z']
  xDict = None
  return v, x, y, z
#==========================================================#
def U_hd( fn, cl=1, direction=False ):
  ut, xu, yu, zu = readVar( fn, 'u_xy', cl )
  vt, xv, yv, zv = readVar( fn, 'v_xy', cl )
  x = xv[:-1]; y = yu[:-1]; z = 0.5*(zu+zv)
  uc = 0.5*( ut[:,:,:-1,1:] + ut[:,:,:-1,:-1] )
  vc = 0.5*( vt[:,:,1:,:-1] + ut[:,:,:-1,:-1] )
  if( direction ):
    v = np.arctan( vc/(uc+1.E-5) ) * (180./np.pi)
  else:
    a = np.arctan( vc/(uc+1.E-5) )
    v = uc * np.cos(a) + vc * np.sin(a)
    
  return v, x, y, z

#==========================================================#
helpStr = '''
Diff mode: 
'd': root mean square diff (RMSD), 
'r': RMSD (relative delta), 
's': RMSD (scaled delta), 
'n': root normalized mean square diff.,
'f': fractional bias
'v': geometric variance
'''
methodList = ['d','r', 's','n','f','v','R']


parser = argparse.ArgumentParser(prog='compareNetCdf2D.py')
parser.add_argument("-f1", "--filename1",type=str, help="Name of the first (ref) input NETCDF file.")
parser.add_argument("-f2", "--filename2",type=str, help="Name of the second input NETCDF file.")
parser.add_argument("-v", "--varname",  type=str, default='u',\
  help="Name of the variable in NETCDF file. Default='u' ")
parser.add_argument("-v0", "--vref", type=float, nargs=2, default=[0.,0.],\
  help="Reference values 'v0' in v+ = (v - v0)/v* for -f1 and -f2. Default = [0,0]")
parser.add_argument("-vs", "--vstar", type=float, nargs=2, default=[1.,1.],\
  help="Characteristic value 'v*' in v+ = (v - v0)/v* for -f1 and -f2. Default = [1.,1.]")
parser.add_argument("-c", "--coarsen", type=int, nargs=2, default=[1,1],\
  help="Factor for coarsening the -f1 and -f2 data when read from file. Default = [1,1]")
parser.add_argument("-m","--mode", type=str, default='d', choices=methodList,\
  help=helpStr)
parser.add_argument("-w", "--writeFile", action="store_true", default=False,\
  help="Write the root-mean-square of the differences to a file.")
parser.add_argument("-nxx1", "--nexclx1", type=int, nargs=2, default=[None,None],\
  help="For -f1, exclude the [first,last] number of nodes from analysis in x-direction.")
parser.add_argument("-nxy1", "--nexcly1", type=int, nargs=2, default=[None,None],\
  help="For -f1, exclude the [first,last] number of nodes from analysis in y-direction.")
parser.add_argument("-nxx2", "--nexclx2", type=int, nargs=2, default=[None,None],\
  help="For -f2, exclude the [first,last] number of nodes from analysis in x-direction.")
parser.add_argument("-nxy2", "--nexcly2", type=int, nargs=2, default=[None,None],\
  help="For -f2, exclude the [first,last] number of nodes from analysis in y-direction.")
parser.add_argument("-xs", "--exclsmall", help="Exclude values below |0.01| from analysis.",\
  action="store_true", default=False)
parser.add_argument("-p", "--printOn", help="Print the numpy array data.",\
  action="store_true", default=False)
parser.add_argument("-s", "--save", action="store_true", default=False,\
  help="Save figures. Default=False")
parser.add_argument("--lims", help="User specified limits.", action="store_true", default=False)
parser.add_argument("--grid", help="Turn on grid.", action="store_true", default=False)
args = parser.parse_args()    

#==========================================================#
# Rename ... that's all.
f1       = args.filename1      # './DATA_2D_XY_AV_NETCDF_N02-1.nc'
f2       = args.filename2      # './DATA_2D_XY_AV_NETCDF_N02-2.nc'
varname  = args.varname
v0       = np.array(args.vref )
vs       = np.array(args.vstar)
cl       = np.array(args.coarsen)
mode     = args.mode
nxx1     = args.nexclx1
nxy1     = args.nexcly1
nxx2     = args.nexclx2
nxy2     = args.nexcly2
exclSmall= args.exclsmall 
writeFile= args.writeFile
printOn  = args.printOn
saveOn   = args.save
limsOn   = args.lims
gridOn   = args.grid

#----------------------------------------------------------#

Sdict = {'d':'RMSD','s':'RMSD (scaled)','r':'RMSD (rel)','n':'RNMSD','f':'FB',\
  'v':'VG','R':'R'}

# Shorter name
vn = varname.split('_')[0]
dirOn   = 'UD' in varname.upper()
horizOn = 'UH' in varname.upper()

# Default for excluded indices is [None,None]. If numerical values are given, 
# the latter needs to be made negative.
if( nxx1.count(None) == 0 ): nxx1[1]*=-1
if( nxy1.count(None) == 0 ): nxy1[1]*=-1
if( nxx2.count(None) == 0 ): nxx2[1]*=-1
if( nxy2.count(None) == 0 ): nxy2[1]*=-1

if( (not horizOn) and (not dirOn) ):
  #print('{}'.format(varname))
  v1, x1, y1, z1 = readVar( f1, varname, cl[0] )
  v2, x2, y2, z2 = readVar( f2, varname, cl[1] )
  
else:
  v1, x1, y1, z1 = U_hd( f1, cl[0], dirOn )  
  v2, x2, y2, z2 = U_hd( f2, cl[1], dirOn )

if( not dirOn ):
  v1 -= v0[0]; v1 /= vs[0]
  v2 -= v0[1]; v2 /= vs[1]

idk = selectFromList( z1 )

if( writeFile ):
  fout = open('{}_d{}.dat'.format(Sdict[mode],vn), 'w')
  fout.write('# file1 = {}, file2 = {}\n'.format(f1, f2))
  fout.write('# z_coord \t {}(d{})\n'.format(Sdict[mode],vn))
  
  #fout.write('{:.2f}\t{:.2e}'.format( z1[k1], dv ))

for k1 in idk:
  
  #k2 = np.where(z2==z1[k1])[0] # This outputs a list 
  k2 = np.where(np.abs(z2-z1[k1])==np.min(np.abs(z2-z1[k1])))[0]
  if( len(k2) == 0 ):
    print(' Coordinate {} not in file {}. Skipping.'.format(z1[k1],f2))
    continue
  else:
    k2 = k2[0]    # Take always the first term
  
  if( len(v1.shape) == 4): v1x  = np.mean(v1[:,k1,nxy1[0]:nxy1[1],nxx1[0]:nxx1[1]], axis=0)
  else:                    v1x  =         v1[  k1,nxy1[0]:nxy1[1],nxx1[0]:nxx1[1]]
    
  if( len(v2.shape) == 4): v2x =  np.mean(v2[:,k2,nxy2[0]:nxy2[1],nxx2[0]:nxx2[1]], axis=0)
  else:                    v2x =          v2[  k2,nxy2[0]:nxy2[1],nxx2[0]:nxx2[1]]

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  dims1  = np.array( v1x.shape )
  dims2  = np.array( v2x.shape )
  if( all( dims1 == dims2 ) ):
    print(' Dimensions of the two datasets match!: dims = {}'.format(dims1))  
  else:
    print(' Caution! Dataset dimensions do not match. dims1 = {} vs. dims2 = {}'.format(dims1, dims2))
    dx1 = (x1[2]-x1[1]); dy1 = (y1[2]-y1[1])
    dx2 = (x2[2]-x2[1]); dy2 = (y2[2]-y2[1])
    rr  = int(np.round(dx2/dx1, decimals=0)); rry = int(np.round(dy2/dy1, decimals=0))
    
    if( rr != rry ): sys.exit(' Resolution ratios are dissimilar. Exiting ...')
    v2f = np.zeros( dims1 ) # Fine resolution
    nc,ec = np.ogrid[ 0:dims1[0] , 0:dims1[1] ]  # northing, easting ... fine resolution
    nf,ef = np.ogrid[ 0:dims1[0] , 0:dims1[1] ]  # northing, easting ... fine resolution
    
    nc=nc//rr; ec=ec//rr  # coarse indices
    
    #nc = nc.astype(int);  ec = ec.astype(int)
    #nf = nf.astype(int);  ef = ef.astype(int)
    
    #np.savetxt('n.dat',np.c_[nf,nc], fmt='%.1f')
    #np.savetxt('e.dat',np.c_[ef.T,ec.T], fmt='%.1f')
    
    # Check bounds
    nf = np.minimum( nf , dims1[0]-1)
    ef = np.minimum( ef , dims1[1]-1)
    nc = np.minimum( nc , dims2[0]-1)
    ec = np.minimum( ec , dims2[1]-1)
    
    # Perform value placement
    v2f[nf, ef] += v2x[nc,ec]; v2x = None
    v2x = v2f
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - #

  if( not np.ma.isMaskedArray(v1x) and not np.ma.isMaskedArray(v2x) ):
    idm = (v1x == v2x)
    v1x = np.ma.masked_array( v1x, mask=idm ) 
    v2x = np.ma.masked_array( v2x, mask=idm )  
    idm = None 
  
  idm = np.ma.getmask(v1x); print(' Nm = {}'.format(np.count_nonzero(idm)))
  
  
  idz = (v2x == 0.0)
  idm += idz
  #idm = sn.binary_dilation(idm); print(' Nm = {}'.format(np.count_nonzero(idm)))
  v1x = np.ma.masked_array( v1x, mask=idm)
  v2x = np.ma.masked_array( v2x, mask=idm)
  #v2x  = np.ma.round( v2x, decimals=2 )
  #v1x  = np.ma.round( v1x, decimals=2 )
  
  
  
  if( exclSmall ):
    # Take values that are above 0.01 or below -0.01
    idx = np.array( (v1x < 5.E-2) )
    #idx = np.array( (v1x > -0.1) )
    m1x = np.ma.getmask(v1x)
    m1x += idx
  
    m2 = np.ma.getmask(v2x)
    m2 += idx
  
  '''
  id2x = np.array( np.abs(v2x) > 1E-2 ) 
  vm2  = np.ma.mean( v2x[id2x] )
  
  id1x = np.array( np.abs(v1x) > 1E-2 )
  vm1  = np.ma.mean( v1x[id1x] )
  dv   = (v2x[id1x] - v1x[id1x] )

  '''
  
  vm1 = np.mean( v1x )
  vm2 = np.mean( v2x )
  print('k={}: vm1 = {}, vm2 = {} '.format(k1,vm1,vm2))
  
  dv  = (v2x - v1x)
  
  
  # NOTE: We are using the desired indices obtained only from the reference (f1) data.
  
  idnn = ~(dv == np.nan )
  N = np.ma.count( np.ravel( dv[idnn] ) )
  print('k={}: Number of good points, N = {}'.format(k1,N))

  if( mode in ['r','s','d'] ):
    if( mode == 'r' ): dv /= np.abs( v1x + 1E-5 )
    if( mode == 's' ): dv /= ( vm1 + 1E-5 )
    #if( mode == 'd' ):  Keep as is: dv = (v2x - v1x)
    RES = np.sqrt(np.sum(dv**2)/N)
    SkewDiff = (1./N)*np.sum(dv**3) * ( 1./(N-1.)*np.sum(dv**2) )**(-1.5)
    print('{} (d{}) = {}, Sk(d{}) = {} '.format(Sdict[mode], vn , RES, vn, SkewDiff ))
    
  if( mode == 'n'):
    v_min_th = np.sqrt(1.e-5)
    vm1 = np.sign(vm1)*np.maximum( np.abs(vm1), v_min_th )
    vm2 = np.sign(vm2)*np.maximum( np.abs(vm2), v_min_th )
    denom = vm1*vm2
    enum  = np.sum(dv**2)/N
    RES = np.sqrt( enum/np.abs(denom) )
    #print(' enum = {}, denom = {} '.format(enum,denom))
    print('{} (d{}) = {}'.format(Sdict[mode], vn , RES))

  if( mode == 'f'):
    denom_min_th = 1.e-2
    dv = (vm2 - vm1)  # Note, mean values
    enum = np.maximum( dv, 1.e-3 )
    denom = 0.5*(np.abs(vm2)+np.abs(vm1))
    denom = np.sign(denom)*np.maximum( np.abs(denom), denom_min_th )
    RES = dv/denom
    #print(' enum = {}, denom = {} '.format(dv,denom))
    print('{} (d{}) = {}'.format(Sdict[mode], vn , RES))
    
  if( mode == 'v'):
    v_min_th = 1.e-1
    dv =  np.log( np.maximum( np.abs(v2x), v_min_th)/(np.maximum( np.abs(v1x), v_min_th )) )
    RES = np.exp( np.sum(dv**2)/N )
    print('{} (d{}) = {}'.format(Sdict[mode], vn , RES))
    
  if( mode == 'R'):
    dv = (v1x-vm1)*(v2x-vm2)
    RES = (np.sum(dv)/N)/(np.std(v1x)*np.std(v2x))
    print('{} (d{}) = {}'.format(Sdict[mode], vn , RES))
  
  if( writeFile ):
    fout.write('{:.2f}\t{:.2e}\n'.format( z1[k1], RES ))


  if( printOn ):
    dimsf  = np.array( np.shape( dv ) )
    xydims = dimsf
    figDims = 13.*(xydims[::-1].astype(float)/np.max(xydims))
    fig = plt.figure(num=1, figsize=figDims)
    labelStr = '({0}_2 - {0}_1)(z={1} m)'.format(vn, z1[k1])
    fig = addImagePlot( fig, dv[::-1,:], labelStr, gridOn, limsOn )
    
    
    fig2 = plt.figure(num=2, figsize=figDims)
    lbl = '(Ref {0})(z={1} m)'.format(vn, z1[k1])
    fig2 = addImagePlot( fig2, v1x[::-1,:], lbl, gridOn, limsOn )
    
    fig3 = plt.figure(num=3, figsize=figDims)
    lbl = '(f2 {0})(z={1} m)'.format(vn, z1[k1])
    fig3 = addImagePlot( fig3, v2x[::-1,:], lbl, gridOn, limsOn )
    
    #fig3 = plt.figure(num=3)
    #plt.hist( np.ravel(dv[idnn]), bins=25, \
    #  normed=True, log=True, histtype=u'bar', label=labelStr )
    
    if( saveOn ):
      figname = '{}_{}_z{}.jpg'.format(Sdict[mode],vn, int(z1[k1]))
      print(' Saving = {}'.format(figname))
      fig.savefig( figname, format='jpg', dpi=150)
      fig2.savefig( 'REF_'+figname, format='jpg', dpi=150)
      fig3.savefig( 'F2_'+figname, format='jpg', dpi=150)
      
      #fig3.savefig( figname.replace("RES","Hist"), format='jpg', dpi=150)
    plt.show()

if( writeFile ): fout.close()
