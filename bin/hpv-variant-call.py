#!/usr/bin/env python
'''
Created on 27/10/2016

@author: sium
'''
from __future__ import print_function


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

__doc__="""Variant caller for HPV project.

Usage:
    hpv-variant-call.py <BAM> (--chromosome <name> | --auto) [--reference <name>]
    hpv-variant-call.py <BAM> (--chromosome <name> | --auto) [--reference <name>] [--start=<number>] [--end=<number>]
    hpv-variant-call.py (-h | --help)
    hpv-variant-call.py --version

Arguments:
    BAM                                          BAM or SAM File name.
    -c <name>, --chromosome <name>               The name of the chromosome.
    -r <name>, --reference <name>                Reference FASTA file.
    -s <number>, --start <number>                Start position [default : 0]
    -e <number>, --end <number>                  End position

Options:
    -a --auto                          Autodetect chromosome name (with highest coverage) to be fetched. 
    -h --help                          Show this screen.
    --version                          Show version.


"""


#prevent sigpipe error
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)
#########


import pysam
from collections import Counter
from docopt import docopt

import sys
import functools


def hpv_variant_table_create(bam_file,chromosome,reference_filename,start,end):



    samfile = pysam.AlignmentFile(bam_file)




    if arguments['--auto']:

        hpv_chromosomes=list(filter(lambda x: x.find("HPV") >=0,samfile.references)) #find HPV chromosomes
        the_list_of_chromosome_counts=list(map(lambda chr: [chr, samfile.count(chr)],hpv_chromosomes)) #estimate HPV chromosome coverages
        autodetected_chromosome=functools.reduce(lambda x,y: x if x[1] > y[1] >= 0 else y,the_list_of_chromosome_counts) #find the highest coverage
        print ("The contig with the highest coverage is %s for the BAM file, %s " % (autodetected_chromosome[0],bam_file), file=sys.stderr)
        chromosome=autodetected_chromosome[0]

    if reference_filename is None:
        sequence = None
    else:
        fastafile= pysam.FastaFile(reference_filename)
        sequence = fastafile.fetch(chromosome)

    print("chr\tposition\treference\tcoverage\tA\tG\tC\tT\tdeletion")




    start= int(0 if start is None else start) #start position of the fetched location
    end=   int(samfile.lengths[samfile.references.index(chromosome)]) if end is None else end #calculate the end by using the chromosome name




    for position in range(start,end):
        position_counter = Counter()
        position_coverage= 0
        for pileupcolumn in samfile.pileup(chromosome, position, position+1, truncate=True,max_depth=1000000000):
            position_coverage=pileupcolumn.n
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_refskip:
                    if not pileupread.is_del:
                        position_counter[pileupread.alignment.query_sequence[pileupread.query_position]] +=1
                    else:
                        position_counter["deletion"]+=1

                else:
                    position_counter["deletion"] += 1

        print ("{chromosome}\t{position}\t{reference}\t{coverage}\t{A}\t{G}\t{C}\t{T}\t{deletion}".format(chromosome=chromosome,position=position+1,
                   reference= 'NA' if sequence is None else sequence[position],
                                                                       coverage=position_coverage-position_counter["deletion"],A=position_counter["A"],
                                                           G=position_counter["G"],C=position_counter["C"],
                                                           T=position_counter["T"],deletion=position_counter["deletion"]
                                                           ))


def main():

    hpv_variant_table_create(arguments['<BAM>'],arguments['--chromosome'],arguments['--reference'],arguments['--start'],arguments['--end'])


if __name__ == '__main__':
    arguments = docopt(__doc__, version='Variant caller for HPV project.')
    main()
