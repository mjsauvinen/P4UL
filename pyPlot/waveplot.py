#!/usr/bin/env python3
import sys
import pylab as pl
from utilities import filesFromList, selectFromList
from plotTools import userLabels, addToPlot
import argparse

# = # = # = # = # = Arguments and their parsing # = # = # = # = # = # 

parser = argparse.ArgumentParser()
parser.add_argument("strKey", type=str, help="Search string for collecting files.", default=None)
parser.add_argument("--select", help="Select the desired plots.", action="store_true", default=False)
args = parser.parse_args()


selections = [" Center Coords vs. Time ", " Umean vs. Time ", " Width vs. Time "]
xlb = [ "Time", "Time" ,  "Time" ]
ylb = [ "Xmean","Umean",  "Width"]


try: 
    args = parser.parse_args()
except: 
    args = parser.parse_args([0])   # Send strKey=0
    args.strKey = raw_input(" Enter search string: ")
    if( not args.strKey ): sys.exit(1)

while 1:
    
    fileNos, fileList = filesFromList( args.strKey ); print "fileList = {}".format(fileList)
    
    if(args.select):
        pId = selectFromList( selections )
    else:
        pId = range(len(selections))
    
    fig = []; lblList = []
    for j in range(len(pId)): 
        fig.append(pl.figure(num=j, figsize=(10.,8.)))
    
    for fn in fileNos:
        lblList.append(fileList[fn])
        
        for j in range(len(pId)):
            d = pl.loadtxt(fileList[fn])
            """
            Cols:   0       1      2      3       4       5     6
            Vars: Center  Umean   Umax   Umin    Width    Id   time
            """
            
            if( pId[j] == 0 ): x = d[:,6]; y = d[:,0]
            if( pId[j] == 1 ): x = d[:,6]; y = d[:,1]
            if( pId[j] == 2 ): x = d[:,6]; y = d[:,4]
            fig[j] = addToPlot( fig[j], x, y, fileList[fn], [selections[j],xlb[j],ylb[j]], 0)
    
    for j in range(len(pId)):
        ax, = fig[j].get_axes(); lns = ax.get_lines()
        ax.legend(lns , lblList , loc=1)
        
    pl.show()









  
