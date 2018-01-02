#!/usr/bin/env python
'''
Created on 1/03/2014

Stockholm format, secondary structure merger.
Secondary structures from RNAalifold.


@author: suu13
'''

from __future__ import print_function
import re
from docopt import docopt

__author__ = 'sium'

__licence__="""
MIT License

Copyright (c) 2017 Sinan Ugur Umu (SUU) sinanugur@gmail.com

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

__doc__="""Merge Stockholm alignments and RNAalifold structure to create a CM file.

Usage:
    stockholm-ss-merger.py <stockholm> <structure>
    stockholm-ss-merger.py (-h | --help)
    stockholm-ss-merger.py --version

Arguments:
    stockholm                           An alignment file in Stockholm format.
    structure                           A structure prediction by RNAalifold

Options:
    -h --help                          Show this screen.
    --version                          Show version.

"""


def RNA_Structure(RNAalifold): #get RNA structure from RNAaligment file
    R_Text=RNAalifold.read()
    return ((re.search(r"(.*) \(",R_Text)).group(1).strip())
    
    

def SS_Merger(Stockholm,RNAstructure):
    stockholm_lines=Stockholm.readlines()
    
    ID=[]
    W_Start=0
    W_Stop=0
    space_count=0
    for line in stockholm_lines:
        if(line[0]!="#" and len(line.split())==2):
            if(line.split()[0] in ID):
                print ("#=GC SS_cons" +' '*space_count + RNAstructure[W_Start:W_Stop] +"\n")
                W_Start=W_Start+W_Stop
                ID=[]
            else:
                ID.append(line.split()[0])
                space_count=len(line.strip())-len("#=GC SS_cons")-len(line.split()[1])
                W_Stop=len(line.split()[1])

        if (line[0:2]=="//"):
            print ("#=GC SS_cons" +' '*space_count + RNAstructure[W_Start:] +"\n")
            print (line.strip())
        else:
            if line.strip():
                print (line.strip())
            



def main():
    with open(arguments['<stockholm>']) as Stockholm, open(arguments['<structure>']) as RNAalifold:
        SS_Merger(Stockholm,RNA_Structure(RNAalifold))



if __name__ == '__main__':
    arguments = docopt(__doc__, version='0.95')
    main()
