#!/usr/bin/bash
# #####################################################
# A simple bash script to clean up the run directory  #
# Author: Jukka-Pekka Keskinen,                       #
#         jukka-pekka.keskinen@fmi.fi                 #
#         Finnish Meteorological Institute            #
###################MMXXI###############################

k=f
getopts ":hi" o

case "${o}" in
    h) 
	echo "Usage: $0 [-ih]"
	echo "This scripts cleans a PALM run directory from the files created by an earlier run."
	echo " "
	echo "  -h     Show this help."
	echo "  -i     Prompt before every removal."
	exit;;
    i)
	k=i
	;;
esac

echo "Cleaning current directory"

rm -$k DEBUG_*
rm -$k DATA_*
rm -$k PROGRESS*
rm -$k RUN_CONTROL*
rm -$k HEADER*
rm -$k CPU_MEASURES*
rm -$k *.err
rm -$k *.out
rm -$k NO_COMBINE*

