import paraview.simple as ps
import sys
import numpy as np
from plotTools import extractFromCSV

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def extractParaviewLineData( caseDir , varList, coords , calcDict,  csvOutputVar , time ):
    try:
        reader = ps.OpenFOAMReader( FileName=caseDir+'/case.foam')
    except:
        ps.Connect()
        reader = ps.OpenFOAMReader( FileName=caseDir+'/case.foam')
        
    try:    reader.CellArrays= varList
    except: print("Variables {} were not found in results. Exiting ...".format(varList))
        
    reader.MeshRegions = ['internalMesh']
    
    if ( time in reader.TimestepValues ):
        reader.SMProxy.UpdatePipeline(time)
        reader.UpdatePipelineInformation()
        #view = ps.CreateRenderView()
        #view.ViewTime = time
    else:
        print("Time-directory {} does not exist. Exiting ...".format(time))
        sys.exit(1)
       
    if( calcDict ): 
        calc = ps.Calculator(reader)
        calc.Function= calcDict["Function"]
        calc.ResultArrayName = calcDict["ResultArrayName"]

    plotLine = ps.PlotOverLine( Source="High Resolution Line Source" )
    plotLine.Source.Resolution = 450
    plotLine.Source.Point1 = coords[0]# [0.0, 0.0008, 0.0]
    plotLine.Source.Point2 = coords[1]# [0.59, 0.0008, 0.0]
    
    
    filePath = caseDir+"/tmp.csv"
    writer = ps.CreateWriter(filePath, plotLine)
    writer.WriteAllTimeSteps = False
    writer.UseScientificNotation = 1
    #print("Writing file tmp.csv ... ")
    writer.UpdatePipeline()
    #print("Done!")
    ps.Delete(reader); ps.Delete(calc)#; ps.Delete(writer)
    ps.Delete(plotLine)
    reader   = None; del reader
    calc     = None; del calc
    plotLine = None; del plotLine
    del writer
    
    if( np.mod( (100.*time) , 1 ) == 0 ):
        print("Disconnecting ... ")
        ps.Disconnect()  # Clear memory occasionally.

    data = extractFromCSV("tmp.csv", csvOutputVar )
    return data

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


