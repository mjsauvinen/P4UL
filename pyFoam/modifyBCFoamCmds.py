#!/usr/bin/python3

import sys 
import argparse
import subprocess as sb

# = = = = = = = = = = = = = = = = = = = = = = = = = = #
def goFwd():
  Question = "\n Continue? <Enter>=Yes , <a-z>=No :"
  q = raw_input(Question)
  if( q != ''): 
    sys.exit(1)

  return True

# = = = = = = = = = = = = = = = = = = = = = = = = = = #

sepStr = ' # = # = # = # = # = # = # = # = '
parser = argparse.ArgumentParser()
parser.add_argument("-t1", "--time", help="Current time step.",\
  default=None)
parser.add_argument("-t2", "--endTime", help="Next endTime in controlDict.",\
  default=None)
parser.add_argument("-v", "--values", help="U_BC values: Old, New", nargs=2,\
  default=[None, None])
parser.add_argument("-x", "--execute", help="Execute commands.",\
  action="store_true", default=False)
parser.add_argument("--noZip", help="Do not unzip and zip files.",\
  action="store_true", default=False)

  
args = parser.parse_args()    
t1   = args.time; t2 = args.endTime
vals = args.values
print args.values

cntrFile = "./system/controlDict"

cmds = []
if( not args.noZip ):
  cmds.append("gunzip -r ./processor*/{}/U.gz".format(t1))

cmds.append("grep -r uniform ./processor*/{}/U".format(t1))
cmds.append("sedFindAndReplace.py -f ./processor*/{0}/U  -s1 \"{1}\" -s2 \"{2}\" ".\
  format(t1,vals[0],vals[1]))

cmds.append("grep -r uniform ./processor*/{}/U".format(t1))

if( not args.noZip ):
  cmds.append("gzip -r ./processor*/{}/U".format(t1))


tStr1 = "endTime    {}".format(t1); tStr2 = "endTime    {}".format(t2)
cmds.append("sedFindAndReplace.py -f {0} -s1 \"{1}\" -s2 \"{2}\" ".\
  format(cntrFile,tStr1,tStr2))
  
for cmd in cmds:
  print cmd + "\n"

if( args.execute  ):
  goFwd()
  for cmd in cmds:
    print cmd + "\n"
    pipe = sb.Popen( cmd , shell=True, stdout=sb.PIPE).stdout
    output = pipe.read()
    print output
    goFwd()
    


