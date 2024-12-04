import numpy as np
import time
import datetime
import argparse

#==========================================================#
parser = argparse.ArgumentParser(prog='date2time.py')
parser.add_argument("-f", "--filename",type=str,\
  help="Name of the data input file.")
parser.add_argument("-fo", "--fileout",type=str,\
  help="Name of the output file.")  
parser.add_argument("-n", "--ncol", metavar='ncol',type=int,\
  help="Column number containing date as string. (First column is zero.)")
parser.add_argument("-cc", "--ccol", metavar='ccol',type=int, nargs='+',\
  help="Column numbers to copy into the same output file.")
parser.add_argument("-d", "--dateformat",type=str, default='%d/%m/%Y %H:%M:%S',\
  help="Dateformat. Default=%%d/%%m/%%Y %%H:%%M:%%S ")  # additional % is needed to escape
parser.add_argument("-sr", "--skiprows", metavar='skiprows',type=int, default=0,\
  help=" Number of rows to be skipped. Default=0.")

args = parser.parse_args()
#==========================================================#
filename   = args.filename
fileout    = args.fileout
dateformat = args.dateformat
ncol       = args.ncol 
ccol       = tuple( args.ccol ) # columns to copy with time data
sr         = args.skiprows  

try:    dat = np.loadtxt(filename, skiprows=sr, usecols=ncol, dtype=str, delimiter=',')
except: dat = np.loadtxt(filename, skiprows=sr, usecols=ncol, dtype=str )

N = len(dat)
t = np.zeros(N)

for i in range(N):
  t[i] = time.mktime(datetime.datetime.strptime(dat[i], dateformat ).timetuple())

t -= t[0]
dat = None 

try:    dc = np.loadtxt(filename, skiprows=sr, usecols=ccol, delimiter=',')
except: dc = np.loadtxt(filename, skiprows=sr, usecols=ccol )

np.savetxt(fileout, np.c_[ t , dc ], fmt='%e' )
