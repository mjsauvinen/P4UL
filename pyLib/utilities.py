import os, sys
import operator
import glob
import numpy as np

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def writeLog( parser , args, skip=False ):
  try:
    if( not skip ):
      import time
      td = time.localtime()
      progname = parser.prog
      fstr = '.'+progname.split('.py')[0]
      fstr += '_{}h{}min{}sec_{:02d}{:02d}{}.log'.format(\
        td.tm_hour,td.tm_min,td.tm_sec,td.tm_mday,\
        td.tm_mon,td.tm_year)
  
      fl = open(fstr, 'w')
      logStr = progname
      for arg, value in sorted(vars(args).items()):
        logStr+=" --{}  {} ".format(arg, value)
      fl.write( logStr )
      fl.write( '\n' )
      fl.close()
      logStr = fstr = None
  except:
    pass    

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def sortTimes(resultDir="./"):
  dirList = os.listdir(resultDir)         # Obtain the list of directories
  #dirList.sort()                   # Sort the list - this won't work	properly for floats (in string format)
  # We need a more robust sorting which uses numerical values
  # This is accomplished by using tuples.
  dirTuple = []
  for dirStr in dirList:
    try:
      dirTuple.append((dirStr, float(dirStr)) ) # (str,float)
    except ValueError:
      pass
  # Sort according to the numerical values: key=operator.itemgetter(1)
  sortedTuple = sorted(dirTuple, key=operator.itemgetter(1))
  sortedList = map(operator.itemgetter(0), sortedTuple) # Extract a list of sorted strings

  return sortedList

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def selectFromList( L ):
    n = 0
    for entry in L:
        print " # [{}]: {}".format(n,entry)
        n+=1
    
    #print "\n Enter the selection number(s): \n"
    Ids = []
    try:
        e = input(" Selection number(s): ")
        if(isinstance(e,tuple)): Ids.extend(e)
        elif(isinstance(e,int)): Ids.append(e)
        else: sys.exit("Invalid entry. Exiting ...")
    except:
        try: input(" Select All? [1-9]=> Yes, [Empty]=> No: ")
        except: sys.exit("Exiting Program ...")
    
    if( len(Ids) == 0 ): Ids.extend(range(len(L)))
    return Ids
    
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def filesFromList( searchStr ):
  print " Extracting files with search string (or path): %s"%searchStr
  fileList = []
  files = glob.glob(searchStr)   # obtain the list of files 
  fileList.extend(files)              # Lists are iterable
  fileList.sort()                     # Sort the list alphabetically

  n = 0
  for f in fileList:
    print " # ["+str(n)+"]: "+ str(f)
    n+=1

  print """
     Enter file number(s), use comma as separator:
     Example 1: File Numbers = 1
     Example 2: File Numbers = 0,2,3,
  """
  fileNos = []
  try:
    e = input(" File Numbers = ")
    if(isinstance(e,tuple)): fileNos.extend(e)
    elif(isinstance(e,int)): fileNos.append(e)
  except:
    try:
      select=input(" Select All? [1-9]=> Yes, [Empty]=> No: ")
    except:
      print ' Exiting program. '
      sys.exit(1)

  if( len(fileNos) == 0 ):
    fileNos.extend(range(len(fileList)))  

  return fileNos, fileList

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def removeEntriesFromList(L, delList):
  for entry in delList:
    try: L.remove(entry)
    except: pass
  
  return L

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
def reformatNumberInWord(word, separator):
  '''
  Search for number in the given string and change its format
  such that a basic sorting algorithm works on it.
  For example: 21 -> 0021.
  '''
  wrdList = word.split(separator)
  for i in xrange(len(wrdList)):
    try:
      iwrd = int(wrdList[i])
      wrdList[i]  = '{:04d}'.format(iwrd)
    except:
      pass
    
  # Put together the new filename 
  word_new = wrdList[0]
  for i in xrange(1,len(wrdList)):
    word_new += separator+wrdList[i]
    
  return word_new

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def extractMatchingTerms( trialList , targetList , verbose=False ):
    xList = []
    for x in trialList:
        if ( x in targetList ):
            xList.append(x)
        else:
            if( verbose ): print 'term = {} not present in the target list.'.format(x)
        
    if( len(xList) == 0 and verbose ): print "Returning an empty list"
    return xList

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def basicAnalysis( x, xStr, xRef, printOn ):
  x_mean = np.mean(x);  x_max = np.max(x)
  x_min  = np.min(x);   x_std = np.std(x)
  N  = float(len(x))
  dx = x-xRef
  x_rms  = np.sqrt( np.sum( x**2 )/N )
  dx_rms = np.sqrt( np.sum( dx**2 )/N )
  
  if( printOn ):
    pStr = '''
    mean({0})= {1}
    max({0}) = {2}
    min({0}) = {3}
    std({0}) = {4}
    rms({0}) = {5}
    drms({0}-{0}Ref) = {6}
    '''.format(xStr, x_mean, x_max, x_min, x_std, x_rms, dx_rms)
    print pStr
    
  return x_mean, x_max, x_min, x_std, x_rms

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
def openStlFile( solidName ):

  solidName = solidName.split('.')[0]  # Take away the possible .stl
  header = 'solid {}'.format(solidName)
  fx=file(solidName+'.stl' , 'w')
  fx.write(header)
  
  return fx

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
def closeStlFile( fs, solidName ):
  
  # Take away the possible .stl
  solidName = solidName.split('.')[0]  # Take away the possible .stl
  footer = '\nendsolid {}\n'.format(solidName)
  fs.write(footer)
  print ' Closing file {}.'.format(solidName)
  fs.close()
  print ' Done ! '

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
def writeStlFacet(fl, nv, v1, v2, v3 ):
  
  fstr= '''
facet normal {0} {1} {2}
outer loop
vertex {3} {4} {5}
vertex {6} {7} {8}
vertex {9} {10} {11}
endloop
endfacet'''.format(nv[0],nv[1],nv[2], v1[0],v1[1],v1[2], v2[0],v2[1],v2[2],\
             v3[0],v3[1],v3[2])
  
  fl.write(fstr)
  
  return fl

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def vtkWriteHeaderAndGridStructured2d( X, Y, Z, fileName, dataStr ):
  nPoints = X.size
  irows = len(X[:,0]); jcols = len(X[0,:])
  header = '# vtk DataFile Version 2.0\n'\
    +'{}\n'.format(dataStr)\
    +'ASCII\n'\
    +'DATASET STRUCTURED_GRID\n'\
    +'DIMENSIONS {}  {}  {}\n'.format(jcols, irows, 1)\
    +'POINTS {} float\n'.format( nPoints )
  
  print ' Writing header for file {} ...'.format( fileName )
  fileName = fileName.split('.vtk')[0]+'-2D.vtk'
  f = open(fileName, 'w')
  f.write( header )

  print ' Writing mesh data for file {} ...'.format( fileName )
  for i in xrange(irows):
    for j in xrange(jcols):
      s = '{0:.1f}\t{1:.1f}\t{2:.1f}\n'.format( X[i,j], Y[i,j], Z[i,j] )
      f.write(s)

  return f

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def vtkWritePointDataHeader( fx, V, nVars ):
  nPoints = V.size
  header ='POINT_DATA {}\n'.format(nPoints)\
    +'FIELD attributes {}\n'.format(nVars)
  
  fx.write(header)
  return fx

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def vtkWritePointDataStructured2D( fx, V, X, vStr ):
  # Check that dimensions agree
  irows  = len(X[:,0]); jcols  = len(X[0,:])
  icheck = len(V[:,0]); jcheck = len(V[0,:])
  if( irows != icheck or jcols != jcheck ):
    sys.exit("dim(V) /= dim(X). Exiting ...")
  

  try:
    nPoints = X.size
    FieldData ='{0} 1 {1} float\n'.format(vStr, nPoints)
    print(' Writing {} field data ...'.format(vStr))
    
    fx.write(FieldData)
    for i in xrange(irows):
      for j in xrange(jcols):
        s = '{0:12.4e} '.format(V[i,j])
        fx.write(s)
    print(' ... done!')
  except:
    pass
  
  return fx


# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def vtkWriteDataStructured2d( V, X, Y, Z, fileName, dataStr ):
  
  f_vtk = vtkWriteHeaderAndGridStructured2d( X, Y, Z, fileName, dataStr )
  
  try:
    f_vtk = vtkWritePointDataHeader( f_vtk, V, 1 )
    f_vtk = vtkWritePointDataStructured2D( f_vtk, V, X, dataStr )
  except:
    pass
  f_vtk.close()
  print(' Writing VTK-data complete!')

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


def vtkWriteUnsPointData( V, X, Y, Z, filename ):
  nPoints = X.size
  irows = len(X[:,0]); jcols = len(X[0,:])
  header = '# vtk DataFile Version 2.0\n'\
    +'Footprint Point Data\n'\
    +'ASCII\n'\
    +'DATASET UNSTRUCTURED_GRID\n'\
    +'POINTS {} float\n'.format( nPoints )
  
  pointdata ='POINT_DATA {}\n'.format(nPoints)\
    +'SCALARS fp float 1\n'\
    +'LOOKUP_TABLE fp\n'.format(nPoints)
  
  print ' Writing file {} ...'.format( filename )
  
  filename = filename.split('.vtk')[0]+'.vtk'
  f = open(filename, 'w')
  f.write( header )
  for i in xrange(irows):
    for j in xrange(jcols):
      s = '{0:.1f}\t{1:.1f}\t{2:.1f}\n'.format( X[i,j], Y[i,j], Z[i,j] )
      f.write(s)
  
  f.write(pointdata)
  for i in xrange(irows):
    for j in xrange(jcols):
      s = '{0:.2f} '.format(V[i,j])
      f.write(s)
  f.close()
  print ' ... done!'

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*