# README #

### What is this repository for? ###

* This repository contains a wide range of python scripts for data pre- and post-processing.

### How do I get set up? ###

############################
Mandatory steps to make the python-scripts work in Linux:

1) Install the following packages: 

* python3-numpy

* python3-matplotlib

* python3-scipy

* python3-gdal 


2) Go to $HOME directory and clone the repository

cd $HOME

git clone https://github.com/mjsauvinen/P4UL.git


3) Add to $HOME/.bashrc

export PATH=$HOME/P4UL:$PATH

export PATH=$HOME/P4UL/pyPlot:$PATH

export PATH=$HOME/P4UL/pyFootprint:$PATH

export PATH=$HOME/P4UL/pyRaster:$PATH

export PATH=$HOME/P4UL/pyUtils:$PATH

export PATH=$HOME/P4UL/pyFoam:$PATH

export PATH=$HOME/P4UL/pyAnalyze:$PATH

export PATH=$HOME/P4UL/pyMisc:$PATH

export PATH=$HOME/P4UL/pyNetCDF:$PATH


export PYTHONPATH=$HOME/P4UL/pyLib/:$PYTHONPATH


4) While in $HOME/P4UL, run the following command to make the scripts executables:

chmod u+x py*/*.py


Now the scripts can be run everywhere in your system.


#########################
### Who do I talk to? ###

* Mikko Auvinen, Finnish Meteorological Institute.

#########################
##### Contributors ######

* Mikko Auvinen, Finnish Meteorological Institute

* Jukka-Pekka Keskinen, Finnish Meteorological Institute

* Sasu Karttunen, University of Helsinki

* Mona Kurppa, University of Helsinki

