#!/usr/bin/env python
'''
Created on 14/05/2013

@author: suu13
'''

#prevent sigpipe error
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)
#########

from __future__ import print_function
import argparse
from Bio import SeqIO



def keyword_finder_from_fasta_headers_nonexact(txt_file,fasta_file): #function to find nonexact match
    with open(txt_file) as txt:
        txt_lines = txt.readlines()
    ids = filter(None, map(lambda x: x.strip(), txt_lines))  # filter endlines and empty items from the text file
    file_sequences = SeqIO.parse(fasta_file, "fasta")
    file_sequences_dictionary={}

    for i in file_sequences:
        file_sequences_dictionary[i.id] = i

    if args.reverse==False: #print the intersection

        disjoint=filter(lambda x: reduce(lambda a,b: a or b,map(lambda y: x.find(y) >= 0,ids)),file_sequences_dictionary.keys())
        map(lambda x: file_sequences_dictionary.pop(x), disjoint) #pop out the non intersection items
        map(lambda i: print (">%s\n%s" % (i,file_sequences_dictionary[i].seq)),file_sequences_dictionary.keys())
                
    else: #print the non-intersection disjoint

        union=filter(lambda x: not reduce(lambda a,b: a or b,map(lambda y: x.find(y) >= 0,ids)),file_sequences_dictionary.keys())
        map(lambda x: file_sequences_dictionary.pop(x), union) #pop out the intersection items
        map(lambda i: print (">%s\n%s" % (i,file_sequences_dictionary[i].seq)),file_sequences_dictionary.keys())


def keyword_finder_from_fasta_headers(txt_file,fasta_file): #exact match

    with open(txt_file) as txt:
        txt_lines = txt.readlines()
    ids = filter(None, map(lambda x: x.strip(), txt_lines))  # filter endlines and empty items from the text file

    file_sequences = SeqIO.parse(fasta_file, "fasta")
    file_sequences_dictionary={}
    for i in file_sequences:
        file_sequences_dictionary[i.id] = i


    if args.reverse==False: #print the intersection

        disjoint=filter(lambda x:x not in ids,file_sequences_dictionary.keys())
        map(lambda x: file_sequences_dictionary.pop(x), disjoint) #pop out the items
        map(lambda i: print (">%s\n%s" % (i,file_sequences_dictionary[i].seq)),file_sequences_dictionary.keys())

    else: #print the non-intersection disjoint

        union=filter(lambda x:x in ids,file_sequences_dictionary.keys())
        map(lambda x: file_sequences_dictionary.pop(x), union) #pop out the intersection items
        map(lambda i: print (">%s\n%s" % (i,file_sequences_dictionary[i].seq)),file_sequences_dictionary.keys())



def main():
    if args.nonexact is None:
        keyword_finder_from_fasta_headers(args.text,args.fasta)
    elif args.nonexact is not None:
        keyword_finder_from_fasta_headers_nonexact(args.text,args.fasta)
    else:
        pass

    
        
  

if __name__ == '__main__':
    Argument_Parser=argparse.ArgumentParser(prog="extract_from_FASTA.py")
    Argument_Parser.add_argument('-fasta',type=str,help="FASTA file to extract sequences.",required=True)
    Argument_Parser.add_argument('-text',type=str,help="A text file with one tag or keyword or ID per line.",required=True)
    Argument_Parser.add_argument('-reverse',action='store_true',help="If this one is ON, remove the sequences with matched IDs.")
    Argument_Parser.add_argument('--nonexact',type=int,help="Non-Exact match, by default exact match.",choices=[1,2],default=None) #1 original seq id, 2 assign the new seq id
    args=Argument_Parser.parse_args()
    
    main()
    