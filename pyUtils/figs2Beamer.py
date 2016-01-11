#!/usr/bin/env python
import sys
import glob
import argparse

#======= ARGUMENTS ========================================#
parser = argparse.ArgumentParser()
parser.add_argument("figKey", help="Search string for collecting pictures.", nargs='?', default="jpg")
parser.add_argument("-t", "--title", help="Title for the slides.", type=str, default="Title")
parser.add_argument("-st", "--subtitle", help="Subtitle for the slides.", type=str, default=" ")
args = parser.parse_args()    
#==========================================================#

searchStr = "*"+args.figKey+"*"
titleStr      = args.title
subtitleStr   = args.subtitle

cmd = "\\includegraphics[width=0.8\\textwidth]{{./Pics/{}}}"

btxt='''
% =*=*=*= FRAME =*=*=*==*=*=*=
\\begin{{frame}}{{ {} }}{{ {} }}

 %\\begin{{columns}}[t]
 %\\column{{0.5\\textwidth}}
 %\\begin{{itemize}}
 %  \\item text
 %\\end{{itemize}}
 %\\column{{0.5\\textwidth}}

 \\begin{{block}}{{ {} }}
   \\centering 
   \\includegraphics[width=0.8\\textwidth]{{./Pics/{}}}
 \\end{{block}}
 
 %\\end{{columns}}
\\end{{frame}}\n
'''

fileList = []
files = glob.glob(searchStr)        # obtain the list of files 
print(' files = {} '.format(files))
fileList.extend(files)              # Lists are iterable
fileList.sort()                     # Sort the list alphabetically

for f in fileList:
  print( btxt.format(titleStr, subtitleStr, f, f) )