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
    hpv-variant-call.py <BAM> [<OUTCSV1>] (--chromosome <name> | --auto) [--reference <name>] [--start=<number>] [--end=<number>] [--transformed] [--discordant]
    hpv-variant-call.py <BAM> <OUTCSV1> <OUTCSV2> --discordant (--chromosome <name> | --auto) [--reference <name>] [--start=<number>] [--end=<number>] [--transformed]
    hpv-variant-call.py <BAM> <FASTA> <BED> [--chromosome <name> | --auto] [--reference <name>] [--start=<number>] [--end=<number>]
    hpv-variant-call.py (-h | --help)
    hpv-variant-call.py --version

Arguments:
    BAM                                          BAM or SAM File name.
    FASTA                                        Output FASTA file name for soft clipped sequences.
    BED                                          Output tab-seperated BED file name for soft clipped sequences.
    OUTCSV1                                      Write regular CSV output into a file, not STDOUT.
    OUTCSV2                                      If given, write discordant CSV into this file.
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


def auto_detect_chromosome_by_coverage(samfile,bam_file):
    hpv_chromosomes = list(filter(lambda x: x.find("HPV") >= 0, samfile.references))  # find HPV chromosomes
    the_list_of_chromosome_counts = list(
        map(lambda chr: [chr, samfile.count(chr)], hpv_chromosomes))  # estimate HPV chromosome coverages
    autodetected_chromosome = functools.reduce(lambda x, y: x if x[1] > y[1] >= 0 else y,
                                               the_list_of_chromosome_counts)  # find the highest coverage
    print("The contig with the highest coverage is %s for the BAM file, %s " % (autodetected_chromosome[0], bam_file),
          file=sys.stderr)

    return(autodetected_chromosome[0])



def auto_detect_hpv_type_from_file_name(samfile,bam_file):

    hpv_name=search('(HPV[0-9]+)',bam_file).group(1)
    hpv_regex = compile("\(" + hpv_name + "\)")

    autodetected_chromosome = list(filter(lambda x: search(hpv_regex,x), samfile.references))  # find HPV chromosome
    
    print("The contig with the highest coverage is %s for the BAM file, %s " % (autodetected_chromosome[0], bam_file),
                          file=sys.stderr)

    return (autodetected_chromosome[0])



def function_position_counter(pileupread,position_counter):
    if not pileupread.is_refskip:
        if not pileupread.is_del:
            position_counter[pileupread.alignment.query_sequence[pileupread.query_position]] += 1
        else:
            position_counter["deletion"] += 1

    else:
        position_counter["skip"] += 1



def hpv_variant_table_create(bam_file,chromosome,reference_filename,start,end,csv1,csv2):

    samfile = pysam.AlignmentFile(bam_file)

    if arguments['--auto']:
        chromosome = auto_detect_hpv_type_from_file_name(samfile,bam_file)

    if reference_filename is None:
        sequence = None


    else:
        #fastafile= pysam.FastaFile(reference_filename)
        #sequence = fastafile.fetch(chromosome)

        for record in SeqIO.parse(reference_filename,"fasta"):
            if record.id == chromosome:
                sequence=str(record.seq)
                break




    start= int(0 if start is None else start) #start position of the fetched location
    end=   int(samfile.lengths[samfile.references.index(chromosome)]) if end is None else int(end) #calculate the end by using the chromosome name
    length=int(samfile.lengths[samfile.references.index(chromosome)])

    second_half=length - floor(length/2) +1
    first_half=floor(length/2 -1)

    function_transformed_position = lambda position: int(
        position + 1 + first_half) if position + 1 <= second_half else  int(position + 1 - second_half)



    print("chr\tposition\treference\tcoverage\tA\tG\tC\tT\tdeletion\tskip",
          file= csv1 if csv1 else sys.stdout)


    if csv2:
        print("chr\tposition\treference\tcoverage\tA\tG\tC\tT\tdeletion\tskip",
          file=csv2)


    for position in range(start,end):
        position_counter = Counter()
        discordant_counter = Counter()
        position_coverage= 0
        for pileupcolumn in samfile.pileup(chromosome, position, position+1, truncate=True,max_depth=1000000000):
            position_coverage=pileupcolumn.n
            for pileupread in pileupcolumn.pileups:
                function_position_counter(pileupread,position_counter)
                if arguments['--discordant'] and (pileupread.alignment.reference_name != pileupread.alignment.next_reference_name):
                    function_position_counter(pileupread,discordant_counter)



        if not arguments['--transformed']:  # is this a shifted genome, no
            pos = position + 1
        else:
            pos = function_transformed_position(position)

        if csv2:
            print_variant_csv_files(position_counter,chromosome,sequence, position, pos, csv1)
            print_variant_csv_files(discordant_counter,chromosome,sequence,position, pos, csv2 )
        elif arguments['--discordant']:
            print_variant_csv_files(discordant_counter,chromosome,sequence,position,pos,csv1 if csv1 else sys.stdout)
        else:
            print_variant_csv_files(position_counter,chromosome,sequence,position,pos,csv1 if csv1 else sys.stdout)



def print_variant_csv_files(position_counter,chromosome,sequence,position,pos,where_to_print):


    print("{chromosome}\t{position}\t{reference}\t{coverage}\t{A}\t{G}\t{C}\t{T}\t{deletion}\t{skip}".format(
        chromosome=chromosome, position=pos,
        reference='NA' if sequence is None else sequence[position],
        coverage=position_counter["A"] + position_counter["G"] + position_counter["C"] + position_counter["T"],
        A=position_counter["A"],
        G=position_counter["G"],
        C=position_counter["C"],
        T=position_counter["T"], deletion=position_counter["deletion"], skip=position_counter['skip']
        ),file=where_to_print)


def fetch_soft_clipped(bam_file,chromosome,start,end,fasta_file,tsv_file):

    samfile = pysam.AlignmentFile(bam_file)

    if arguments['--auto']:
        chromosomes = list(auto_detect_hpv_type_from_file_name(samfile,bam_file))
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
                    #seq_position=0
                    #read_aligned_pairs=read.get_aligned_pairs()
                    #for i in read.cigartuples:
                        #if i[0] == 4 and i[1] >= 10: #detect soft clipped, 4 is for soft clip


                    if match(cigarsoft, read.cigarstring):
                        size=int(match(cigarsoft, read.cigarstring).group(1))
                        sequence=read.seq[0:size]
                    else:
                        size = int(search(cigarsoft, read.cigarstring).group(1))
                        sequence = read.seq[-size:]

                    if read.is_reverse:
                        sequence=str(Seq(sequence).reverse_complement()) #take reverse complement if on opposite strand


                    print (">{read_id}\n{sequence}".format(read_id=read.query_name,sequence=sequence),file=fasta)
                    feat_start = read.reference_start if match(cigarsoft,read.cigarstring) else read.reference_end

                    print ("{ref_id}\t{feat_start}\t{feat_end}\t{name}\t{score}\t{strand}".format(ref_id=read.reference_name,
                                                                                               feat_start=feat_start,
                                                                                               feat_end=feat_start+size,
                                                                                name=read.query_name,score=1,strand="."),file=tsv)

                        #break
                        #elif i[0] != 3: #3 is for Ns
                        #elif i[0] != 3:  # 3 is for Ns
                        #    seq_position=seq_position + i[1]



                else:
                    pass




def main():

    if arguments['<FASTA>']:
        fetch_soft_clipped(arguments['<BAM>'],arguments['--chromosome'],arguments['--start'],arguments['--end'],arguments['<FASTA>'],arguments['<BED>'])
    else:

        if arguments['<OUTCSV2>']:
            with open(arguments["<OUTCSV1>"], "w") as csv1, open(arguments["<OUTCSV2>"], "w") as csv2:
                hpv_variant_table_create(arguments['<BAM>'], arguments['--chromosome'], arguments['--reference'],
                                     arguments['--start'], arguments['--end'], csv1, csv2)

        elif arguments['<OUTCSV1>']:
            with open(arguments["<OUTCSV1>"], "w") as csv1:
                hpv_variant_table_create(arguments['<BAM>'], arguments['--chromosome'], arguments['--reference'],
                                         arguments['--start'], arguments['--end'], csv1, csv2=None)

        else:
            hpv_variant_table_create(arguments['<BAM>'], arguments['--chromosome'], arguments['--reference'],
                                     arguments['--start'], arguments['--end'], csv1=None, csv2=None)


if __name__ == '__main__':
    arguments = docopt(__doc__, version='0.95')
    main()
