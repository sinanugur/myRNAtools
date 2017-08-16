#!/usr/bin/env python
'''
Created on 13/01/2017

@author: sium
'''
from __future__ import print_function


__author__ = 'sium'


__doc__="""RNA sequencing reads count tool based on HTSeq.

Usage:
    rnaseq_tool.py <BAM> --gff <file> [--filter] [--feature_type=<gene>] [--identity_attribute=<ID>]
    rnaseq_tool.py <BAM> [--isomircount | --hairpincount]
    rnaseq_tool.py <BAM> --gff <file> [--filter] [--coverage] [--alignment <file>] [--feature_type=<gene>] [--identity_attribute=<ID>]
    rnaseq_tool.py (-h | --help)
    rnaseq_tool.py --version

Arguments:
    BAM                                          BAM or SAM File name.
    -g <file>, --gff <file>                      A GFF File.
    -f <feature>, --feature_type <feature>       Which feature to count [default: miRNA]
    -i <ID>, --identity_attribute <ID>           Which identifier to compile the results [default: ID]
    -a <file>, --alignment <file>                An alignment file.

Options:
    -h --help                          Show this screen.
    --version                          Show version.
    --collapsed                        The reads are collapsed.
    --isomircount                      Count isomiRs.
    --hairpincount                     Count hairpins.
    --filter                           Filter multimapped reads.
    --coverage                         Calculate coverage of features.




"""


#prevent sigpipe error
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)
#########


import HTSeq
import collections
from docopt import docopt

from Bio import SeqIO
from difflib import ndiff
from re import search

from Bio.SeqRecord import SeqRecord
from Bio import AlignIO,SeqIO

import pandas as pd



#chromosome mapping counts
def count_in_chromosomes(bam_file): # a very classical counter for BAM files.
    bam_reader = bam_or_sam_reader(bam_file)
    chrome_counts = collections.defaultdict( lambda: 0 ) #add zero count to a new chromosome
    for a_read in bam_reader:
        if a_read.aligned:
            chrome_counts[ a_read.iv.chrom ] += 1

    for chromename in sorted(chrome_counts.keys()):
        print ("%s\t%d" % (chromename,chrome_counts[chromename]))


def bam_or_sam_reader(bam_file):
    try:
        bam_reader = HTSeq.BAM_Reader(bam_file)
        return bam_reader
    except:
        bam_reader = HTSeq.SAM_Reader(bam_file)
        return bam_reader



def type_of_difference(sequence1, sequence2): #detect isomir type
    classification=[]
    difference=''.join([i.split()[0] for i in list(ndiff(sequence1, sequence2))])


    if bool(search('^\++|^\-+',difference)): #5' end isomir
        classification.append('5end')

    if bool(search('\++$|\-+$',difference)): #3' end isomir
        classification.append('3end')

    if bool(search('\++|\-+',difference.strip('+|-'))) : #SNP isomir
        classification.append('SNP')

    return classification


def count_hairpins(bam_file,mirna_hairpin_sequences_file):
    bam_reader = bam_or_sam_reader(bam_file)
    mirna_hairpin_sequences = list(SeqIO.parse(mirna_hairpin_sequences_file, "fasta"))

    mirna_hairpin_dictionary={}
    mirna_hairpin_counts = collections.defaultdict(lambda: 0)

    for i in mirna_hairpin_sequences:
        mirna_hairpin_dictionary[i.id] = i


    for a_read in bam_reader:
        if a_read.aligned:
            a_read_id = a_read.read.name
            mirna_hairpin_counts[a_read.iv.chrom] += int(a_read_id.split("-")[1])  # mirna hairpin

    for i in sorted(mirna_hairpin_counts.keys()):
        print ("%s\t%d\t%s\thairpin") % (i, mirna_hairpin_counts[i],mirna_hairpin_dictionary[i].seq)



def count_mirnas_isomirs(bam_file,mirna_sequences_file):
    bam_reader = bam_or_sam_reader(bam_file)
    mirna_sequences = list(SeqIO.parse(mirna_sequences_file, "fasta"))

    mirna_dictionary = {}
    isomir_dictionary = {}
    mirna_exact_counts = collections.defaultdict(lambda: 0)
    mirna_isomir_counts = collections.defaultdict(lambda: 0)

    for i in mirna_sequences:
        mirna_dictionary[i.id] = i



    for a_read in bam_reader:
        if a_read.aligned:

            #which_mirna=process.extractOne(a_read.iv.chrom, mirna_dictionary.keys(), scorer=fuzz.partial_ratio)[0]
            a_read_id=a_read.read.name
            which_mirna=a_read.iv.chrom

            which_mirna_sequence=mirna_dictionary[which_mirna].seq
            if which_mirna_sequence == a_read.read.seq: #mirna exact
                mirna_exact_counts[which_mirna] += int(a_read_id.split("-")[1])
            else:

                mirna_isomir_counts[a_read.read.seq]+=int(a_read_id.split("-")[1])
                isomir_dictionary[a_read.read.seq] = SeqRecord(seq=a_read.read.seq,
                    id=which_mirna,
                    name=which_mirna,
                    description=';'.join(type_of_difference(which_mirna_sequence,a_read.read.seq)))

    for i in sorted(mirna_exact_counts.keys()):
        print ("%s\t%d\t%s\tmature.exact") % (i, mirna_exact_counts[i],mirna_dictionary[i].seq)
    for i in sorted(mirna_isomir_counts.keys()):
        print ("%s\t%d\t%s\tisomir.%s") % \
              (isomir_dictionary[i].id, mirna_isomir_counts[i],isomir_dictionary[i].seq,isomir_dictionary[i].description)

    return




#An isoform report tool from GFF annotations. Only for collapsed reads.
def htseq_feature_count_gff_to_report_reads(bam_file,gff_file):
    bam_reader = bam_or_sam_reader(bam_file)
    gff_reader = HTSeq.GFF_Reader(gff_file)
    #coverage = HTSeq.GenomicArray("auto", stranded=True, typecode="i")
    features_array = HTSeq.GenomicArrayOfSets("auto",stranded=True)


    feature_to_count = arguments['--feature_type']
    feature_gff_id_to_group =  arguments['--identity_attribute']



    for feature in gff_reader:
        if feature.type == feature_to_count: #for example miRNA, exon etc.
            features_array[feature.iv] += feature.attr[feature_gff_id_to_group] # for example ID, Alias etc.

    #this part is important because ambigous reads have to be treated carefully. for now featurecounts method. check oneone entry
    for a_read in bam_reader:
        if a_read.aligned and (int(a_read.optional_field('XN')) == 0):  #remove out the reads with ambigous bases
            iset = None
            for iv2, step_set in features_array[a_read.iv].steps():
                if iset is None:
                    iset = step_set.copy()
                else:
                    iset.update(step_set)

            if len(iset) >= 1: # I think by default this is equal to one only, which means only one annotation is allowed for a single read.
                for i,f in enumerate(iset):
                    if not arguments['--filter']: #no need for filtering, print all alignments
                        print ("%s\t%s\t%d" % (str(a_read.read.seq,"utf-8"), f, int(a_read.read.name.split("-")[1])))
                    else:
                        try:
                            XS = int(a_read.optional_field('XS'))
                            AS = int(a_read.optional_field('AS'))
                            if arguments['--filter'] and AS >= XS: #count only equally good alignments, not necessarily unique
                                print ("%s\t%s\t%d" % (str(a_read.read.seq,"utf-8"),f,int(a_read.read.name.split("-")[1])))
                            else:
                                pass
                        except: #the alignment is unique, report it
                            print ("%s\t%s\t%d" % (str(a_read.read.seq,"utf-8"), f, int(a_read.read.name.split("-")[1])))

    return


# Report coverages of regions in a GFF file and if there is an alignment, use it.
def htseq_feature_count_gff_to_report_coverage(bam_file, gff_file):
    bam_reader = bam_or_sam_reader(bam_file)
    gff_reader = HTSeq.GFF_Reader(gff_file)
    coverage = HTSeq.GenomicArray("auto", stranded=True, typecode="i")
    features_array = HTSeq.GenomicArrayOfSets("auto", stranded=True)

    alignment_record = {}

    feature_to_count = arguments['--feature_type']
    feature_gff_id_to_group = arguments['--identity_attribute']

    for feature in gff_reader:
        if feature.type == feature_to_count:  # for example miRNA, exon etc.
            features_array[feature.iv] += feature.attr[feature_gff_id_to_group]  # for example ID, Alias etc.

    # this part is important because ambigous reads have to be treated carefully. for now featurecounts method. check oneone entry
    for a_read in bam_reader:
        if a_read.aligned and (int(a_read.optional_field('XN')) == 0):  #remove out the reads with ambigous bases
            iset = None
            for iv2, step_set in features_array[a_read.iv].steps():
                if iset is None:
                    iset = step_set.copy()
                else:
                    iset.update(step_set)

            if len(iset) >= 1:  # I think by default this is equal to one only, which means only one annotation is allowed for a single read.
                for i, f in enumerate(iset):
                    if not arguments['--filter']:  # no need for filtering, print all alignments
                        coverage[a_read.iv] += int(a_read.read.name.split("-")[1])
                    else:
                        try:
                            XS = int(a_read.optional_field('XS'))
                            AS = int(a_read.optional_field('AS'))
                            if arguments['--filter'] and AS >= XS:  # count only equally good alignments, not necessarily unique
                                coverage[a_read.iv] += int(a_read.read.name.split("-")[1])
                            else:
                                pass
                        except:  # the alignment is unique, report it
                            coverage[a_read.iv] += int(a_read.read.name.split("-")[1])




    if ['--alignment']:
        for record in AlignIO.read(arguments['--alignment'], "stockholm"):
            alignment_record[record.id]=record



        for feature in gff_reader:
            if feature.type == feature_to_count:  # for example miRNA, exon etc.
                id="%s:%s-%s(%s)" % (feature.iv.chrom,feature.iv.start,feature.iv.end,feature.iv.strand)

                feature_coverage=list(coverage[feature.iv])
                print (id,end="")
                for nucleotide in str(alignment_record[id].seq):
                    if nucleotide is not '-' and nucleotide is not '.':
                        print ("\t%s" % str(feature_coverage.pop(0)),end="")
                    else:
                        print("\t%s" % "NA", end="")

                print ("")

    else:
        for feature in gff_reader:
            if feature.type == feature_to_count:  # for example miRNA, exon etc.
                print ("%s:%s-%s(%s)\t%s") % (str(feature.iv.chrom,"utf-8"), str(feature.iv.start,"utf-8"), str(feature.iv.end,"utf-8"), str(feature.iv.strand,"utf-8"),
                                    list(coverage[feature.iv]))  # example od chr19:4724634-4724707(-))

    return



def main():


    if arguments['--gff'] and not arguments['--coverage']: #if it is not None
        htseq_feature_count_gff_to_report_reads(arguments['<BAM>'],arguments['--gff'])
    elif arguments['--gff'] and arguments['--coverage']:
        htseq_feature_count_gff_to_report_coverage(arguments['<BAM>'], arguments['--gff'])
    else:
        pass


if __name__ == '__main__':
    arguments = docopt(__doc__, version='RNA sequencing isomiR and isoform extracting tool.')
    main()