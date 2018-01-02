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
    hpv-variant-call.py <BAM> (--chromosome <name> | --auto) [--reference <name>] [--transformed]
    hpv-variant-call.py <BAM> (--chromosome <name> | --auto) [--reference <name>] [--start=<number>] [--end=<number>] [--transformed]
    hpv-variant-call.py <BAM> <FASTA> <TSV> [--chromosome <name> | --auto] [--reference <name>] [--start=<number>] [--end=<number>]
    hpv-variant-call.py (-h | --help)
    hpv-variant-call.py --version

Arguments:
    BAM                                          BAM or SAM File name.
    FASTA                                        Output FASTA file name for soft clipped sequences.
    TSV                                          Output Tab-seperated text file name for soft clipped sequences.
    -c <name>, --chromosome <name>               The name of the chromosome.
    -r <name>, --reference <name>                Reference FASTA file.
    -s <number>, --start <number>                Start position [default : 0]
    -e <number>, --end <number>                  End position

Options:
    -a --auto                          Autodetect chromosome name (with highest coverage) to be fetched. 
    -t --transformed                   Mapped HPV genomes are transformed.
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
from math import floor

from Bio import SeqIO
from Bio.Seq import Seq
from re import search
from re import match
from re import compile


def auto_detect_chromosome(samfile,bam_file):
    hpv_chromosomes = list(filter(lambda x: x.find("HPV") >= 0, samfile.references))  # find HPV chromosomes
    the_list_of_chromosome_counts = list(
        map(lambda chr: [chr, samfile.count(chr)], hpv_chromosomes))  # estimate HPV chromosome coverages
    autodetected_chromosome = functools.reduce(lambda x, y: x if x[1] > y[1] >= 0 else y,
                                               the_list_of_chromosome_counts)  # find the highest coverage
    print("The contig with the highest coverage is %s for the BAM file, %s " % (autodetected_chromosome[0], bam_file),
          file=sys.stderr)

    return(autodetected_chromosome[0])


def hpv_variant_table_create(bam_file,chromosome,reference_filename,start,end):

    samfile = pysam.AlignmentFile(bam_file)

    if arguments['--auto']:
        chromosome = auto_detect_chromosome(samfile,bam_file)

    if reference_filename is None:
        sequence = None


    else:
        #fastafile= pysam.FastaFile(reference_filename)
        #sequence = fastafile.fetch(chromosome)

        for record in SeqIO.parse(reference_filename,"fasta"):
            if record.id == chromosome:
                sequence=str(record.seq)
                break




    print("chr\tposition\treference\tcoverage\tA\tG\tC\tT\tdeletion\tskip")




    start= int(0 if start is None else start) #start position of the fetched location
    end=   int(samfile.lengths[samfile.references.index(chromosome)]) if end is None else int(end) #calculate the end by using the chromosome name
    length=int(samfile.lengths[samfile.references.index(chromosome)])

    second_half=length - floor(length/2) +1
    first_half=floor(length/2 -1)


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
                    position_counter["skip"] += 1


        function_position=lambda position: int(position + 1 + first_half) if position + 1 <= second_half else  int(position + 1 - second_half)

        pos=position+1 if not arguments['--transformed'] else function_position(position)


        print ("{chromosome}\t{position}\t{reference}\t{coverage}\t{A}\t{G}\t{C}\t{T}\t{deletion}\t{skip}".format(chromosome=chromosome,position=pos,
                   reference= 'NA' if sequence is None else sequence[position],
                                                                       coverage=position_counter["A"] + position_counter["G"] +position_counter["C"] +position_counter["T"] ,
                                                                        A=position_counter["A"],
                                                           G=position_counter["G"],
                                                        C=position_counter["C"],
                                                           T=position_counter["T"],deletion=position_counter["deletion"],skip=position_counter['skip']
                                                           ))


def fetch_soft_clipped(bam_file,chromosome,start,end,fasta_file,tsv_file):

    samfile = pysam.AlignmentFile(bam_file)

    if arguments['--auto']:
        chromosomes = list(auto_detect_chromosome(samfile,bam_file))

    elif chromosome is None:
        chromosomes = samfile.references

    else:
        chromosomes = list(chromosome)

    cigarsoft = compile("([1-9][0-9]+)S")



    with open(fasta_file,"w") as fasta,open(tsv_file,"w") as tsv:
        for chromosome in chromosomes:
            start = int(0 if start is None else start)  # start position of the fetched location
            end = int(samfile.lengths[samfile.references.index(chromosome)]) if end is None else int(
                    end)  # calculate the end by using the chromosome name

            for read in samfile.fetch(chromosome,start,end):
                if not read.is_unmapped and search(cigarsoft,read.cigarstring):
                    seq_position=0
                    for i in read.cigartuples:
                        if i[0] == 4 and i[1] >= 10: #detect soft clipped, 4 is for soft clip
                            if read.is_reverse:
                                sequence=str(Seq(read.seq[seq_position:seq_position + i[1]]).reverse_complement()) #take reverse complement if on opposite strand
                            else:
                                sequence=read.seq[seq_position:seq_position + i[1]]
                            print (">{read_id}\n{sequence}".format(read_id=read.query_name,sequence=sequence),file=fasta)
                            feat_start = read.reference_start if match(cigarsoft,read.cigarstring) else read.reference_end

                            print ("{ref_id}\t{feat_start}\t{feat_end}\t{name}\t{score}\t{strand}".format(ref_id=read.reference_name,
                                                                                                       feat_start=feat_start,
                                                                                                       feat_end=feat_start+i[1],
                                                                                        name=read.query_name,score=1,strand="."),file=tsv)

                            break
                        #elif i[0] != 3: #3 is for Ns
                        elif i[0] != 3:  # 3 is for Ns
                            seq_position=seq_position + i[1]



                else:
                    pass


def main():

    if arguments['<FASTA>']:
        fetch_soft_clipped(arguments['<BAM>'],arguments['--chromosome'],arguments['--start'],arguments['--end'],arguments['<FASTA>'],arguments['<TSV>'])
    else:
        hpv_variant_table_create(arguments['<BAM>'],arguments['--chromosome'],arguments['--reference'],arguments['--start'],arguments['--end'])


if __name__ == '__main__':
    arguments = docopt(__doc__, version='0.95')
    main()
