import os,sys,csv,random,re
from pathlib import Path
## Usage: python <scriptname> <txt filewith one queryid per line> <gff3 file>
## Usage: python <scriptname> <txt filewith one genomic coordinate per line> <gff3 file>
## for second usage comment PrintMappingQueryinGff usage and uncomment line below.

## For example:
## python QueryGff3.py <queryfile-with-oneid-perline> <gff3file>
## python QueryGff3.py querycoord_test.txt ~/MAGIC16/annotation/genome8/genes.gff
## python QueryGff3.py queryids.txt ~/MAGIC16/annotation/genome8/genes.gff

## To get bed file for input list of genes from their annot gff3 file
##thimmamp@KW61043:~/PythonPrograms/utilities$ python QueryGff3.py ~/MAGIC16/indels_results/genome33_Specific_genes.txt ~/MAGIC16/annotation/genome33/genes.gff genome33_Specific_genes.bed

def Queryids2Dict(infile):
    outDict = {}
    with open(infile, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            outDict[line]=1
    return outDict

def PrintMappingQueryinGff(indict, annotfile):
    with open(annotfile, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                tokens=line.rstrip('\n').split('\t')
                for k in indict.keys():
                    if k in tokens[-1]:
                        print(line)

def GetGenicregionBedForQueryId(indict, annotfile, outfile):
    # Chr12	cshl_gene	gene	24526624	24527466	.	-	.	ID=gene:OsIR64_12g0020410;biotype=protein_coding;logic_name=cshl_gene

    with open(annotfile, 'r') as f:
        with open(outfile, 'w') as of:
            for line in f:
                if not line.startswith('#'):
                    tokens=line.rstrip('\n').split('\t')
                    for k in indict.keys():
                        if k in tokens[-1]:
                            #print(line)
                            annotid = tokens[-1].split(';')[0].split(':')[1]
                            if k==annotid:
                                #print("Matching :", k, annotid)
                                of.write(tokens[0]+'\t'+tokens[3]+'\t'+tokens[4]+'\t'+annotid+'\t0\t'+tokens[6]+'\n')



def PrintMappingQueryRegioninGff(querydict, annotfile):
    with open(annotfile, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                tokens=line.rstrip('\n').split('\t')
                chromosome = tokens[0]
                start = int(tokens[3])
                end = int(tokens[4])
                for k in querydict.keys():
                    qchr = k.split(':')[0]
                    qstart = int(k.split(':')[1].split('-')[0])
                    qend = int(k.split(':')[1].split('-')[1])
                    if chromosome==qchr and start >= qstart  and end <= qend:
                        print(line)


## call functions
query = Queryids2Dict(sys.argv[1])
#PrintMappingQueryinGff(query, sys.argv[2])
#PrintMappingQueryRegioninGff(query, sys.argv[2])
GetGenicregionBedForQueryId(query, sys.argv[2], sys.argv[3])
