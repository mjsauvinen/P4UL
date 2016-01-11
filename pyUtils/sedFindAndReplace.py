#!/usr/bin/env python
import subprocess as sb
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-s1", "--oldstr", help="Old string",  nargs='+')
parser.add_argument("-s2", "--newstr", help="New string",  nargs='+')
parser.add_argument("-f", "--filenames", help="Filenames", nargs='+' )
args = parser.parse_args()

fileList = args.filenames

ostr = args.oldstr
nstr = args.newstr

print 'old strings: {}'.format(ostr)
print 'new strings: {}'.format(nstr)

if(not len(nstr) == len(ostr)):
    sys.exit(' Number of Old and New strings do not match. Exiting ...')

for fn in fileList:
    for i in xrange(len(ostr)):
        sedCommand = "find {0} -type f -exec sed -i \'s/{1}/{2}/g\' {{}} \\;"\
        .format(fn, ostr[i], nstr[i])

        print sedCommand
        sb.call(sedCommand, shell=True)

print 'Done!'