"""This script will generated the required sample file for mothur make.contigs command. 
Usage: python make_file.py <fastq directory><file name>"""
import subprocess
import sys
import os
import re

raw_output = subprocess.check_output('find {} -type f -name "*.fastq.gz" | grep -v Apple'.format(sys.argv[1]), shell = True)
pathlist = sorted(raw_output.rstrip().split('\n'))

if len(pathlist) % 2 == 0:
    with open(sys.argv[2], 'w') as outfile:
        for i in range(0, len(pathlist), 2):
            filename = os.path.basename(pathlist[i])
            sample = re.sub("_S[0-9]*_L[0-9]*.*","", os.path.basename(pathlist[i]))
            outfile.write('\t'.join([str(sample), pathlist[i], pathlist[i + 1]]) + '\n')
else: 
    print "check: uneven number of files"