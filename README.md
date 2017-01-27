# README #

### What is this repository for? ###

* This repository contains a wide range of python scripts for data pre- and post-processing.

### How do I get set up? ###

############################
Mandatory steps to make the python-scripts work in Linux:

0) Go to $HOME/bin and clone the repository

cd $HOME/bin

git clone https://github.com/mjsauvinen/P4US.git


1) Add to $HOME/.bashrc 

export PATH=$HOME/bin:$PATH

export PATH=$HOME/bin/pyPlot:$PATH

export PATH=$HOME/bin/pyFootprint:$PATH

export PATH=$HOME/bin/pyRaster:$PATH

export PATH=$HOME/bin/pyUtils:$PATH

export PATH=$HOME/bin/pyFoam:$PATH

export PATH=$HOME/bin/pyAnalyze:$PATH

export PATH=$HOME/bin/pyMisc:$PATH

export PATH=$HOME/bin/pyNetCDF:$PATH


export PYTHONPATH=$HOME/bin/pyLib/:$PYTHONPATH


2) While in $HOME/bin run the following command to make the programs executables:

chmod u+x py*/*.py


Now the scripts can be run everywhere in your system.


#########################
### Who do I talk to? ###

* Mikko Auvinen, University of Helsinki / Finnish Meteorological Institute.
