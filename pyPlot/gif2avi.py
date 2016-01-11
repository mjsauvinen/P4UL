#!/usr/bin/env python
import subprocess as sb
import argparse
''' 
Description:


Author: Mikko Auvinen
        mikko.auvinen@helsinki.fi 
        University of Helsinki &
        Finnish Meteorological Institute
'''

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filename", help="Name of input gif file", type=str)
parser.add_argument("-o", "--output", help="Name of output file", type=str, default='output')
parser.add_argument("-s", "--scale", help="Number of pixels in X and Y dir", type=int, nargs=2)
parser.add_argument("-fps", "--framerate", help="Frames per second", type=int, default=15)
parser.add_argument("--mp4", help="MP4 instead of AVI.", action="store_true", default=False)

args = parser.parse_args()

# Rename for convenience.
fps = args.framerate
sx  = args.scale[0] # 1024
sy  = args.scale[1] # 768
fileIn = args.filename
fileOut= args.output

if( args.mp4 ):
    # mencoder, MPG
    cmd = 'mencoder -fps {0} -vf scale={1}:{2} -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell:vbitrate=10000 -o {3}mp4 {4}gif'.format(\
    fps,sx,sy,fileOut.rstrip('mp4'), fileIn.rstrip('gif') )
else:
    # mencoder, AVI
    cmd = 'mencoder -fps {0} -vf scale={1}:{2} -ovc lavc -lavcopts vcodec=msmpeg4v2:autoaspect:vbitrate=10000 -o {3}avi {4}gif'.format(\
    fps,sx,sy,fileOut.rstrip('avi'), fileIn.rstrip('gif') )
    
print( cmd )
sb.call(cmd, shell=True)
print( 'Completed!' )

