#!/usr/bin/env python
'''
Created on 27/10/2016

@author: sium
'''


__author__ = 'sium'


#prevent sigpipe error
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)
#########


import HTSeq
import collections
import argparse

from Bio import SeqIO
from difflib import ndiff
from re import search

from Bio.SeqRecord import SeqRecord

#chromosome mapping counts
def count_in_chromosomes(bam_file): # a very classical counter for BAM files.
    bam_reader = bam_or_sam_reader(bam_file)
    chrome_counts = collections.defaultdict( lambda: 0 ) #add zero count to a new chromosome
    for a_read in bam_reader:
        if a_read.aligned:
            chrome_counts[ a_read.iv.chrom ] += 1

    for chromename in sorted(chrome_counts.keys()):
        print "%s\t%d" % (chromename,chrome_counts[chromename])


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



def bam_or_sam_reader(bam_file):
    try:
        bam_reader = HTSeq.BAM_Reader(bam_file)
        return bam_reader
    except:
        bam_reader = HTSeq.SAM_Reader(bam_file)
        return bam_reader


#A simple GFF read counter, for non-collapsed reads
def htseq_feature_count_gff(bam_file,gff_file):
    bam_reader = bam_or_sam_reader(bam_file)
    gff_reader = HTSeq.GFF_Reader(gff_file)
    #coverage = HTSeq.GenomicArray("auto", stranded=True, typecode="i")
    features_array = HTSeq.GenomicArrayOfSets("auto",stranded=True)


    #feature_to_count = 'miRNA' if args.feature_type is None else args.feature_type
    #feature_gff_id_to_group =  'ID' if args.identity_attribute is None else args.identity_attribute

    feature_to_count='miRNA'
    feature_gff_id_to_group ='Alias'

    for feature in gff_reader:
        if feature.type == feature_to_count: #for example miRNA
            features_array[feature.iv] += feature.attr[feature_gff_id_to_group] # for example ID, Alias etc


    counts = {} #a dictionary to hold counts
    for feature in gff_reader:
        if feature.type == feature_to_count:
            counts[feature.attr[feature_gff_id_to_group]] = 0



    #this part is important because ambigous reads have to be treated carefully. for now featurecounts method.
    for a_read in bam_reader:
        if a_read.aligned:
            iset = None
            for iv2, step_set in features_array[a_read.iv].steps():
                if iset is None:
                    iset = step_set.copy()
                else:
                    iset.update(step_set)

            if len(iset) >= 1:
                for i,f in enumerate(iset):
                    counts[f] +=1

    for id in sorted(counts.keys()):
        print id, counts[id]

    return


def htseq_read_count_collapsed_reads(bam_file):
    bam_reader = bam_or_sam_reader(bam_file)

    for a_read in bam_reader:
        if a_read.aligned:
            try: #if the read is multi aligned

                if(a_read.option_field('XS') < a_read.option_field('AS')):
                    print ("%s\t%i" % (a_read.read.seq,int(a_read.read.name.split("-")[1])))

            except: #this read is single mapped
                    print ("%s\t%i" % (a_read.read.seq,int(a_read.read.name.split("-")[1])))



def main():

    if args.gff is not None:
        htseq_feature_count_gff(args.bam,args.gff) #featurecount style with non-collapsed reads
    if args.collapsedreadcount is True:
        htseq_read_count_collapsed_reads(args.bam)
    elif args.isomircount is True and args.fasta is not None:
        count_mirnas_isomirs(args.bam,args.fasta)
    elif args.hairpincount is True and args.fasta is not None:
        count_hairpins(args.bam, args.fasta)
    else:
        count_in_chromosomes(args.bam)


    """"
    program_function = {
        'pairfold': lambda: Pairfold_Execute(args.sRNA, args.targetRNA) if args.gff is not None else Pairfold_Execute_Parallel(
            args.sRNA, args.targetRNA),
        'bifold': lambda: bifold_Execute(args.sRNA, args.targetRNA) if args.cpu is None else bifold_Execute_Parallel(
            args.sRNA, args.targetRNA),
        'DuplexFold': lambda: DuplexFold_Execute(args.sRNA, args.targetRNA),
        'RNAup': lambda: RNAup_Execute_Parallel(args.sRNA, args.targetRNA, args.window),
        'AccessFold': lambda: AccessFold_Execute(args.sRNA, args.targetRNA),
        }

    program_function[args.program]()
    """

if __name__ == '__main__':
    Argument_Parser=argparse.ArgumentParser(prog="rnaseq_tool.py")
    Argument_Parser.add_argument('bam',type=str,help="BAM or SAM file.")
    Argument_Parser.add_argument('-gff',type=str,help="GFF file.")
    Argument_Parser.add_argument('-fasta', type=str, help="miRNAs FASTA file; either mature sequnces "
                                                          "for isomir counts or hairpins for hairpin counts")
    Argument_Parser.add_argument('--feature_type', type=str, help="Default = miRNA")
    Argument_Parser.add_argument('--identity_attribute', type=str, help="Default = ID")
    Argument_Parser.add_argument('--collapsedreadcount',action='store_true', help="Count and report all single mapped or higher quality collapsed reads"
                                                                           "(requires collapsed read mapped BAM file)")
    Argument_Parser.add_argument('--isomircount',action='store_true', help="Count and report isomir reads "
                                                                           "(requires mature miRNA mapped BAM file)")
    Argument_Parser.add_argument('--hairpincount', action='store_true', help="Count hairpin "
                                                                             "reads (requires hairpin mapped BAM file)")
    args=Argument_Parser.parse_args()

    main()