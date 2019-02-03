#!/usr/bin/env python

'''
Created on 14/05/2013
Update on 01/02/2019

@author: suu13
'''
from __future__ import print_function



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




_doc_="""Extract selected sequences from a FASTA or FASTQ file

Usage:
    extract_from_FASTA.py --fasta <file> [--text <file>] [--reverse] [--nonexact]
    extract_from_FASTA.py (-h | --help)
    extract_from_FASTA.py --version


Arguments:
    -f <file>, --fasta <file>       A FASTA or FASTQ file of input sequences.
    -t <file>, --text <file>        A text file that contains sequence IDs per line.

Options:
    -h --help                   Show this screen.
    --version                   Show version.
    --reverse                   If this one is ON, remove the sequences with matched IDs and print the rest.
    --nonexact                  Do a non-exact match using find function, by default do an exact match.


"""



#prevent sigpipe error
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)
#########

from Bio import SeqIO
from docopt import docopt
from functools import reduce
import gzip
import sys


def function_read_text_file(txt_file):
    
    if txt_file == None:
        txt_lines = sys.stdin.readlines()
    else:
        with open(txt_file) as txt:
            txt_lines = txt.readlines()
    ids = list(filter(None, list(map(lambda x: x.strip(), txt_lines))))  # filter endlines and empty items from the text file
    return ids

def function_parse_sequence_file(fasta_file):
    file_sequences_dictionary={}
    if fasta_file.lower().endswith(".gz"):
        with gzip.open(fasta_file,"rt") as handle:
            
            if True:
                file_sequences = SeqIO.parse(handle, "fastq")
            else:
                file_sequences = SeqIO.parse(handle, "fasta")

            for i in file_sequences:
                file_sequences_dictionary[i.id] = i
    else:
        with open(fasta_file,"r") as handle:
            if True:
                file_sequences = SeqIO.parse(fasta_file, "fastq")
            else:
                file_sequences = SeqIO.parse(fasta_file, "fasta")
            
            for i in file_sequences:
                file_sequences_dictionary[i.id] = i

    return file_sequences_dictionary

def function_print_output(items_to_print):
    for item in items_to_print:
        #print (">%s\n%s" % (item.id,item.seq))
        SeqIO.write(item,sys.stdout,"fastq")

def keyword_finder_from_fasta_headers_nonexact(txt_file,fasta_file): #function to find nonexact match

    ids=function_read_text_file(txt_file)
    file_sequences_dictionary=function_parse_sequence_file(fasta_file)
    if arguments['--reverse']==False: #print the intersection
        union=list(filter(lambda x: reduce(lambda a,b: a or b,list(map(lambda y: x.find(y) >= 0,ids))),file_sequences_dictionary.keys()))
        items_to_print=map(lambda x: file_sequences_dictionary.pop(x), union) #pop out the items
        function_print_output(items_to_print)

    else: #print the non-intersection disjoint
        disjoint=list(filter(lambda x: not reduce(lambda a,b: a or b,list(map(lambda y: x.find(y) >= 0,ids))),file_sequences_dictionary.keys()))
        items_to_print=map(lambda x: file_sequences_dictionary.pop(x), disjoint) #pop out the intersection items
        function_print_output(items_to_print)

def keyword_finder_from_fasta_headers(txt_file,fasta_file): #exact match

    ids = function_read_text_file(txt_file)
    file_sequences_dictionary=function_parse_sequence_file(fasta_file)
    if arguments['--reverse']==False: #print the intersection
        union=list(filter(lambda x:x in ids,file_sequences_dictionary.keys()))
        items_to_print=map(lambda x: file_sequences_dictionary.pop(x), union) #pop out the items
        function_print_output(items_to_print)

    else: #print the non-intersection disjoint
        disjoint=list(filter(lambda x:x not in ids,file_sequences_dictionary.keys()))
        items_to_print=map(lambda x: file_sequences_dictionary.pop(x), disjoint)
        function_print_output(items_to_print)

def main():

    if(arguments['--nonexact'] is False):
        keyword_finder_from_fasta_headers(arguments['--text'],arguments['--fasta'])
    else:
        keyword_finder_from_fasta_headers_nonexact(arguments['--text'],arguments['--fasta'])


if __name__ == '__main__':
    arguments = docopt(_doc_, version='A sequence extraction tool from a FASTA file 1.5')
    main()

