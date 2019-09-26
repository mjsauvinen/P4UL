#!/usr/bin/env python3

import os, sys
import subprocess
import math
from pylab import *

def plotVelocityProfiles(timeDir, time0 ):
    time = float( timeDir ) 
    rpm = 1470.0                 # Specify this (rpm)
    omegaRad = rpm*(2.*math.pi)/60.
    omegaDeg = omegaRad * (180./math.pi)
    angle = int( round( omegaDeg * (time-time0) , 0) )
    if( angle <= 0 ):
        nRev = angle/359
        angle = -nRev*360 + angle
    elif( angle > 360 ):
        nRev = angle/360
        angle = angle - nRev*360

    print ' time = %f , angle = %03d '%(time, angle)
    fileCFD = timeDir+'/line_U.xy'
    fileEmp = '../mp12/profile_%03d.txt'%angle
    xc, Ux_c, Uy_c = mlab.load( fileCFD,skiprows=0,usecols=(0,1,2),unpack=True )
    xe, Uy_e, Ux_e = mlab.load( fileEmp,skiprows=0,usecols=(0,1,3),unpack=True )
    offSet = -0.0475
    xe[:] = -xe[:]/1000. + offSet

    fig=figure(num=1, figsize=(6.15,6.15));
    ax=fig.add_axes( [0.14, 0.075 , 0.85 , 0.81] ) #[left, up, width, height], fig.add_subplot(111)

	# mp11/mp12
    lines=ax.plot(Ux_c ,xc ,'b', -Uy_c, xc,'g', Ux_e , xe,'ko', Uy_e, xe,'ro', linewidth=4, markersize=11)
    ax.set_ylim([max(xc),min(xc)]); ax.set_xlim([-2.0 , 5.6 ])

	# mp31/mp32
    #lines=ax.plot(-Ux_c ,xc ,'b', -Uy_c, xc, 'g', Uy_e , xe,'ro', Ux_e , xe,'ko',linewidth=4, markersize=11)
    #ax.set_ylim([max(xc),min(xc)]); ax.set_xlim([-2.2 , 6.0 ])

    ax.grid(True)
    ax.set_xlabel('U (m/s)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Z-Coord. vs. Velocity, angle = %03d deg '%angle)

    fig.legend(lines,('OF-dev: Ux','OF-dev: Uy', 'LDA: Ux', 'LDA: Uy'), loc=(0.48,0.62))
	
    fig.suptitle(' OpenFOAM-dev (C) vs. Experiment (W1B) ')

    print ' Saving plot_U_%03d.png'%angle
    fname = 'plot_U_%03d.png'%angle
    fig.savefig(fname)
    ax.cla()# Clear axes  
    fig.clf() # Clear figure window

    if(True):
        fileStr = 'C-U_%03d.dat'%angle
        ff = open( fileStr , 'w'  )
    	for i in range(len(xc)):
    	    ff.write('%f \t %f \t %f \n' %(xc[i], Ux_c[i], Uy_c[i]) )
    	ff.close()
    	print ' file ', fileStr, ' written successfully. '

# - - - - - - - - - - - - - - - - - - - - -

try:
    timeZeroDegStr =  sys.argv[1]
    currentDir     =  sys.argv[2]
except:
    print "Usage:", sys.argv[0], " [time at 0deg]  [current directory] "
    sys.exit(1)

timeZeroDeg = float( timeZeroDegStr )

for f in os.listdir(currentDir):
    if( os.path.isdir(os.path.join(currentDir, f)) ):
        plotVelocityProfiles( f , timeZeroDeg)
	

subprocess.call(" convert plot_U* new_U.miff ", shell=True)
subprocess.call(" rm -f *.png ", shell=True)

#show()



