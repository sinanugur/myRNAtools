#!/usr/bin/env python

'''
Created on 18/04/2013

EMBL2GFF parser, it reads all sequence features and write to attributes

# embl2gff.py -embl file_name

# Unutma EMBL dosyasinin kendisi de ilk sirada annotation olarak cikiyor, dikkat

@author: suu13
'''
import argparse
from Bio import SeqIO
import signal
import sys



def EMBL2GFF(EMBL):
    output_format="%s\tEMBL\t%s\t%d\t%d\t.\t%s\t.\t"
    signal.signal(signal.SIGPIPE, signal.SIG_DFL) #broken pipe hatasindan kurtulmak icin yaptim
    with open(EMBL) as EMBL_file:

        try:
            records=SeqIO.parse(EMBL_file,"gb")
        except:
            records=SeqIO.read(EMBL_file,"embl")

        for record in records:
            for EMBL_features in record.features:
                #EMBL_qualifiers=str(EMBL_features.qualifiers)[1:-1]
                EMBL_qualifiers=EMBL_features.qualifiers
                EMBL_feature_type=EMBL_features.type
                Seq_Start=int(EMBL_features.location.start) + 1 #embl okuyucusu bir sebeple hep 1 eksik okuyor, bunu duzeltmek icin koydum, dikkat et
                Seq_End=int(EMBL_features.location.end)
                if EMBL_features.location.strand == -1:
                    Seq_Strand='-'
                else:
                    Seq_Strand='+'
                try:
                    print (output_format % (record.id,EMBL_feature_type,Seq_Start,Seq_End,Seq_Strand),end="")
                    for k in EMBL_qualifiers:
                        print ("%s \"%s\"; " % (k,EMBL_qualifiers[k][0]),end="")
                    print ('') 
                except IOError:
                    pass
                 
def EMBL2_close_gene_features(EMBL,geneloc):
    output_format="%s\tEMBL\t%s\t%d\t%d\t.\t%s\t.\t%s"
    signal.signal(signal.SIGPIPE, signal.SIG_DFL) #broken pipe hatasindan kurtulmak icin yaptim
    Upstream=[]
    Downstream=[]
    with open(EMBL) as EMBL_file:
        record=SeqIO.read(EMBL_file,"embl")
        for EMBL_features in record.features:
            EMBL_qualifiers=str(EMBL_features.qualifiers)[1:-1]

            EMBL_feature_type=EMBL_features.type
            Seq_Start=int(EMBL_features.location.start) + 1 #embl okuyucusu bir sebeple hep 1 eksik okuyor, bunu duzeltmek icin koydum, dikkat et
            Seq_End=int(EMBL_features.location.end)
            if EMBL_features.location.strand == -1:
                Seq_Strand='-'
                if abs(geneloc[1]-Seq_Start) <=100 and EMBL_feature_type == 'gene':
                    try:
                        #print output_format % (record.id,EMBL_feature_type,Seq_Start,Seq_End,Seq_Strand,EMBL_qualifiers)
                        #sys.stdout.write("Upstream\t%s\tReverse\t%d\t" % (EMBL_qualifiers,geneloc[1]-Seq_Start))
                        Upstream=[EMBL_qualifiers,"Reverse",geneloc[1]-Seq_Start]
                    except IOError:
                        pass
                elif abs(geneloc[0]-Seq_End) <=100 and EMBL_feature_type == 'gene':
                    try:
                        #print output_format % (record.id,EMBL_feature_type,Seq_Start,Seq_End,Seq_Strand,EMBL_qualifiers)
                        #sys.stdout.write("Downstream\t%s\tReverse\t%d\t" % (EMBL_qualifiers,geneloc[0]-Seq_End))
                        Downstream=[EMBL_qualifiers,"Reverse",geneloc[0]-Seq_End]
                    except IOError:
                        pass

            else:
                Seq_Strand='+'
                if abs(geneloc[0]-Seq_End) <=100 and EMBL_feature_type == 'gene':
                    try:
                        #print output_format % (record.id,EMBL_feature_type,Seq_Start,Seq_End,Seq_Strand,EMBL_qualifiers)
                        #sys.stdout.write("Upstream\t%s\tForward\t%d\t" % (EMBL_qualifiers,geneloc[0]-Seq_End))
                        Upstream=[EMBL_qualifiers,"Forward",geneloc[0]-Seq_End]
                    except IOError:
                        pass
                elif abs(geneloc[1]-Seq_Start) <=100 and EMBL_feature_type == 'gene':
                    try:
                        #print output_format % (record.id,EMBL_feature_type,Seq_Start,Seq_End,Seq_Strand,EMBL_qualifiers)
                        #sys.stdout.write("Downstream\t%s\tForward\t%d\t" % (EMBL_qualifiers,geneloc[1]-Seq_Start))
                        Downstream=[EMBL_qualifiers,"Forward",geneloc[1]-Seq_Start]
                    except IOError:
                        pass
    if(len(Upstream) is not 3):
        Upstream=["NA","NA","NA"]
    if(len(Downstream) is not 3):
        Downstream=["NA","NA","NA"]

    sys.stdout.write("%s\t%s\n" % ('\t'.join(str(x) for x in Upstream),'\t'.join(str(x) for x in Downstream)))
    return
def main():
    if args.genelocation is not None:
        EMBL2_close_gene_features(args.embl,args.genelocation)
    else:
        EMBL2GFF(args.embl)



if __name__ == '__main__':
    Argument_Parser=argparse.ArgumentParser(prog="embl2gff.py")
    Argument_Parser.add_argument('-embl',type=str,help="Name of the EMBL file",required=True)
    Argument_Parser.add_argument('-genelocation',type=int,nargs=2,help="Window to find close genes")
    args=Argument_Parser.parse_args()
    main()
    pass
