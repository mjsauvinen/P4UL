#!/usr/bin/python

import sys
import subprocess as sb
import pylab as pl
import argparse
from utilities import filesFromList
from plotTools import userLabels, plotCSV
''' 
Description:


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''

#==========================================================#
parser = argparse.ArgumentParser()
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
parser.add_argument("-yl","--ylims", help="Y-axis limits: [min,max]. Default=[0,10]",\
  type=float,nargs=2,default=[0.,10.])
parser.add_argument("-fn","--figName", help="Name of the (temporary) figures. (default=tmp)",\
  type=str,default="tmp")
parser.add_argument("-fa","--fileAnim", help="Name of the animation file. (default=anim.gif)",\
  type=str,default="anim.gif")
parser.add_argument("-na", "--noAnim", help="Do not make an animation.",\
  action="store_true", default=False) 
args = parser.parse_args()    
#==========================================================#
strKey      = args.strKey
figName     = args.figName
fileAnim    = args.fileAnim
noAnimation = args.noAnim
ylims       = args.ylims

n = 0
fileNos, fileList = filesFromList( "*"+strKey+"*" )

for fn in fileNos:
  n += 1
  pfig = pl.figure(num=1, figsize=(18.,9.))
  pfig = plotCSV( pfig, fileList[fn], args.yx, args.magy, args.reuse )
  pl.grid(True)
  pl.legend(loc=2)
  ax = pfig.get_axes()[0]
  ax.set_ylim( tuple(ylims) )
  figStr = figName+'_{:04d}.png'.format(n)
  print(' Saving fig: {}'.format(figStr))
  pl.savefig(figStr)
  ax.cla()
  pfig.clf()

if( not noAnimation ):
  cmd = 'convert {}_*  {} '.format(figName, fileAnim)
  print ' Executing command: ${}'.format(cmd)
  sb.call(cmd, shell=True)

