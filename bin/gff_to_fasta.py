#!/usr/bin/env python
'''
Created on Apr 27, 2018

@author: sium
'''
# Sequence file parse using a GFF file
# Downloaded sequence file from NCBI

__author__ = 'sium'

__licence__="""
MIT License

Copyright (c) 2018 Sinan Ugur Umu (SUU) sinanugur@gmail.com

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

from docopt import docopt
from Bio import SeqIO
from BCBio import GFF

__doc__="""Convert GFF to FASTA.



Usage:
    gff_to_fasta.py <GFF> <FASTA>
    gff_to_fasta.py (-h | --help)
    gff_to_fasta.py --version

Arguments:
    GFF                                 A GFF file.
    FASTA                               A FASTA sequence file.

Options:
    -h --help                          Show this screen.
    --version                          Show version.

"""



def GFF_to_FASTA(gff_file,fasta_file):

    with open(gff_file) as gff, open(fasta_file) as fasta:

        seq_dict=SeqIO.to_dict(SeqIO.parse(fasta,"fasta"))

        for rec in GFF.parse(gff,base_dict=seq_dict):
            print (rec)
            #print (">{}\t{}\n".format(rec.seqid,rec.seq))

def main():
    GFF_to_FASTA(arguments["<GFF>"],arguments["<FASTA>"])


if __name__ == '__main__':
    arguments = docopt(__doc__, version='0.95')
    main()

