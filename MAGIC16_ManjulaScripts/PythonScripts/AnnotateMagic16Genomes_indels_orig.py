import os,sys,csv,random,re
from pathlib import Path

## Objective is to annotate indels for each genome on MAGIC16 group with IRGSP(genome1 as ref)
## Usage: python AnnotateMagic16Genomes_indels.py ~/MAGIC16/tmp/vcf/testing/genome1vsMagic16genomes.snpeff.SNPs.vcf2


def GetAnnotationtoDictfromGff(annotfile, field):
    outDict={}
    with open(annotfile, mode='r') as f:
        for line in f:
            if not line.startswith('#'):
                tokens=line.rstrip('\n').split('\t')
                if tokens[2]==field:
                    key=tokens[0]+':'+tokens[3]+'-'+tokens[4]
                    outDict[key]=tokens[-1]
                    #print(key,tokens[-1])
                    
    return outDict                

######
## features for parsing
CDS
chromosome
exon
five_prime_UTR
gene
intergenic
intron
mRNA
three_prime_UTR
########

def GetGffFeaturestoDict(gff_file):
    outDict={}
    with open(gff_file, mode='r') as f:
        for line in f:
            if not line.startswith("#"):
                tokens = line.rstrip('\n').split('\t')
                if tokens[2]!="chromosome":
                    key = tokens[0]+':'+tokens[3]+'-'+tokens[4]+'_'+tokens[2]
    return outDict               


def GetGffAnnotationToDict(gff_file):
    outDict={}
    with open(gff_file, mode='r') as f:
        for line in f:
            if not line.startswith("#"):
                tokens = line.rstrip('\n').split('\t')
                if tokens[2]=="gene":
                    gene = tokens[-1].split(';')[0].split(':')[1]
                    outDict[gene]={}
                if tokens[2]=='mRNA':
                    txptid = tokens[-1].split(';')[0].split(':')[1] 
                    # example ID=transcript:Osaus.01G000020_01;Parent=gene:Osaus.01G000020;biotype=protein_coding;transcript_id=Osaus.01G000020_01;canonical_transcript=1
                    geneid = tokens[-1].split(';')[1].split(':')[1]
                    txptloci = tokens[0]+':'+tokens[3]+'-'+tokens[4]
                    txptkey = txptid+';'+txptloci
                    if not txptid in outDict[geneid].keys():
                        outDict[geneid][txptid] = [txptkey]  # get list of components with their chromosomal location
                if tokens[2]=='CDS': #CDS, chromosome, exon, five_prime_UTR, three_prime_UTR
                    # example ID=CDS:Osaus.01G000010_01;Parent=transcript:Osaus.01G000010_01;protein_id=Osaus.01G000010_01
                    cdsid = tokens[-1].split(';')[0].split('=')[1]
                    txptid = tokens[-1].split(';')[1].split(':')[1]
                    geneid = txptid.split('_')[0]
                    cdskey = cdsid+';'+tokens[0]+':'+tokens[3]+'-'+tokens[4]
                    if geneid in outDict.keys() and txptid in outDict[geneid].keys():
                       components = outDict[geneid][txptid]
                       if not cdskey in components:
                          components.append(cdskey)
                if tokens[2]=='exon': #CDS, chromosome, exon, five_prime_UTR, three_prime_UTR
                    # example Parent=transcript:Osaus.01G000020_01;Name=Osaus.01G000020_01.exon.1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Osaus.01G000020_01.exon.1;rank=1
                    exonid = tokens[-1].split(';')[4].split('=')[1]
                    exonumber = tokens[-1].split(';')[-1].split('=')[1]
                    txptid = tokens[-1].split(';')[0].split(':')[1]
                    geneid = txptid.split('_')[0]
                    exonkey = exonid+';'+exonumber+';'+tokens[0]+':'+tokens[3]+'-'+tokens[4]
                    if geneid in outDict.keys() and txptid in outDict[geneid].keys():
                        components = outDict[geneid][txptid]
                        if not exonkey in components:
                            components.append(exonkey)
                if tokens[2]=='three_prime_UTR': #CDS, chromosome, exon, five_prime_UTR, three_prime_UTR
                    # example Chr01   NAM     three_prime_UTR 21657   21789   .       +       .       Parent=transcript:Osaus.01G000020_01
                    txptid = tokens[-1].split(':')[1]
                    geneid = txptid.split('_')[0]
                    threeprimekey = tokens[2]+';'+tokens[0]+':'+tokens[3]+'-'+tokens[4]
                    if geneid in outDict.keys() and txptid in outDict[geneid].keys():
                        components = outDict[geneid][txptid]
                        if not cdskey in components:
                            components.append(threeprimekey)
                if tokens[2]=='five_prime_UTR': #CDS, chromosome, exon, five_prime_UTR, three_prime_UTR
                    # example Chr01   NAM     five_prime_UTR  14336   14625   .       +       .       Parent=transcript:Osaus.01G000020_02
                    txptid = tokens[-1].split(':')[1]
                    geneid = txptid.split('_')[0]
                    fiveprimekey = tokens[2]+';'+tokens[0]+':'+tokens[3]+'-'+tokens[4]
                    if geneid in outDict.keys() and txptid in outDict[geneid].keys():
                        components = outDict[geneid][txptid]
                        if not cdskey in components:
                            components.append(fiveprimekey)
    return outDict                    

def GetAnnotationtoDict(annotfile):
    outDict = {}
    with open(annotfile, mode='r') as f:
        for line in f:
            if not line.startswith("##"):
                line = line.rstrip('\n')
                tokens=line.split('\t')
                #print(tokens[2], tokens[-2])
                tmp=tokens[-2].split(';')
                annset=set()
                for item in tmp:
                    if "ANN=" in item:
                        tos=item.split('|')
                        for annitem in tos:
                            if '_variant' in annitem:
                                annset.add(annitem)
                            if 'gene:' in annitem:
                                annset.add(annitem)
                #print(tokens[2], annset)
                outDict[tokens[2]]=annset
    return outDict          

#def AddAnnotation(infile, annotdict, outfile):
def AddAnnotation(infile, annotdict):
    outDict = {}
    with open(infile, mode='r') as f:
        #with open(outfile, mode='w') as of:
            for line in f:
                tokens=line.rstrip('\n').split('\t')
                for k, v in annotdict.items():
                    chromosome=k.split(':')[0]
                    start = int(k.split(':')[1].split('-')[0])
                    end = int(k.split(':')[1].split('-')[1])
                    #check if SNP's loci lies between annotation loci
                    if tokens[1]==chromosome and int(tokens[2]) >= start and int(tokens[2]) <= end:
                        valtokens = v.split(';')
                        gene=''
                        for item in valtokens:
                            #if "exon_id=" in item:
                            #    gene = item.split('=')[1] ## here is exon = gene
                            if "ID=" in item:
                                gene=item.split(':')[1]
                        locigene=k+'_'+gene        
                        if locigene in outDict.keys():
                            snpset = outDict[locigene]
                            snpset.add(tokens[0])
                            outDict[locigene]=snpset
                        else:
                            snpset = set()
                            snpset.add(tokens[0])
                            outDict[locigene]=snpset
                        #print(line, k, v)
                        #of.write(tokens[0]+'\t'+tokens[1]+'\t'+tokens[2]+'\t'+v+'\n')
    return outDict                   

def writeDict2File(indict, outfilename):
    with open(outfilename, 'w') as of:
        for k,v in indict.items():
            of.write(k+'\t'+str(len(v))+'\t'+str(v)+'\n')

def PrintDict(indict):
    for k, v in indict.items():
        for txptkey, components in indict[k].items():
            print(k, txptkey, components,'\n')
        #print(k,len(v), v, '\n')


def ReadVcfAnnotFiletoDict(infile):
    outDict = {}
    with open(infile,'r') as f:
        for line in f:
            tokens = line.rstrip('\n').split('\t')
            tokens[-1] = tokens[-1].replace('{','')
            tokens[-1] = tokens[-1].replace('}','')
            outDict[tokens[0]] = tokens[-1]
    return outDict


## call functions
genetxptcomponents = GetGffAnnotationToDict(sys.argv[1])
PrintDict(genetxptcomponents)



'''
annotype="gene"
annotdict = GetAnnotationtoDictfromGff(sys.argv[1],annotype) ## "gene", "mRNA", "exon"
outdir="/home/thimmamp/MAGIC16/indels_results/IndividualGenomes/AnnotatedVcfs/"

#print("SampleName\tInsertions\tDeletions\tBoth\tOthers")
for fn in os.listdir(sys.argv[2]):
#    if fn.endswith(".insertions.only.forcmplot.txt"):
    if fn.endswith(".deletions.only.forcmplot.txt"):
        print(os.path.join(sys.argv[2], fn))
        samplename = fn.split('modified')[0]
#        #samplename = fn.split('.')[0]
        outfilename=samplename+"_deletions_mappedtoAnnotation_"+annotype+".txt"
        print("Sample Name: ", samplename, " Output file name: ", outfilename,'\n')
#        #queryGenomeAnnotdict = AddAnnotation(os.path.join(sys.argv[2],fn), annotdict, outfilename)
        queryGenomeAnnotdict = AddAnnotation(os.path.join(sys.argv[2],fn), annotdict)
        writeDict2File(queryGenomeAnnotdict, outdir+outfilename)
        
#FilterVCFforMinCountAllele(sys.argv[1])
#queryGenomeAnnotdict = AddAnnotation(sys.argv[2], annotdict, sys.argv[3])
#for k, v in queryGenomeAnnotdict.items():
#    print(k,len(v),'\n')

## Get dict for each genome from MAGIC16 cohort
'''
## All 16 genomes
######
## replace _deletions with nothing for insertions annotated file

#genome1_dict = ReadVcfAnnotFiletoDict("/home/thimmamp/MAGIC16/indels_results/IndividualGenomes/AnnotatedVcfs/genome1_GJ_temp_IRGSP_deletions_mappedtoAnnotation_gene.txt")
#genome4_dict = ReadVcfAnnotFiletoDict("/home/thimmamp/MAGIC16/indels_results/IndividualGenomes/AnnotatedVcfs/genome4_cA1_N22_deletions_mappedtoAnnotation_gene.txt")
#genome5_dict = ReadVcfAnnotFiletoDict("/home/thimmamp/MAGIC16/indels_results/IndividualGenomes/AnnotatedVcfs/genome5_GJ_trop1_Azucena_deletions_mappedtoAnnotation_gene.txt")
#genome6_dict = ReadVcfAnnotFiletoDict("/home/thimmamp/MAGIC16/indels_results/IndividualGenomes/AnnotatedVcfs/genome6_XI_1B1_IR_64_deletions_mappedtoAnnotation_gene.txt")
#genome7_dict = ReadVcfAnnotFiletoDict("/home/thimmamp/MAGIC16/indels_results/IndividualGenomes/AnnotatedVcfs/genome7_cB_ARC10497_deletions_mappedtoAnnotation_gene.txt")
#genome8_dict = ReadVcfAnnotFiletoDict("/home/thimmamp/MAGIC16/indels_results/IndividualGenomes/AnnotatedVcfs/genome8_XI_2B_LarhaMugad_deletions_mappedtoAnnotation_gene.txt")
#genome9_dict = ReadVcfAnnotFiletoDict("/home/thimmamp/MAGIC16/indels_results/IndividualGenomes/AnnotatedVcfs/genome9_XI_3B2_LiuXu_deletions_mappedtoAnnotation_gene.txt")
#genome10_dict = ReadVcfAnnotFiletoDict("/home/thimmamp/MAGIC16/indels_results/IndividualGenomes/AnnotatedVcfs/genome10_XI_3B1_KhaoYaiGuang_deletions_mappedtoAnnotation_gene.txt")
#genome11_dict = ReadVcfAnnotFiletoDict("/home/thimmamp/MAGIC16/indels_results/IndividualGenomes/AnnotatedVcfs/genome11_XI_3A_LIMA_deletions_mappedtoAnnotation_gene.txt")
#genome12_dict = ReadVcfAnnotFiletoDict("/home/thimmamp/MAGIC16/indels_results/IndividualGenomes/AnnotatedVcfs/genome12_cA2_NATELBORO_deletions_mappedtoAnnotation_gene.txt")
#genome13_dict = ReadVcfAnnotFiletoDict("/home/thimmamp/MAGIC16/indels_results/IndividualGenomes/AnnotatedVcfs/genome13_XI_1B2_PR_106_deletions_mappedtoAnnotation_gene.txt")
#genome14_dict = ReadVcfAnnotFiletoDict("/home/thimmamp/MAGIC16/indels_results/IndividualGenomes/AnnotatedVcfs/genome14_GJ_trop2_KetanNangka_deletions_mappedtoAnnotation_gene.txt")
#genome15_dict = ReadVcfAnnotFiletoDict("/home/thimmamp/MAGIC16/indels_results/IndividualGenomes/AnnotatedVcfs/genome15_GJ_subtrp_ChaoMeo_deletions_mappedtoAnnotation_gene.txt")
#genome16_dict = ReadVcfAnnotFiletoDict("/home/thimmamp/MAGIC16/indels_results/IndividualGenomes/AnnotatedVcfs/genome16_XI_2A_GobolSail_deletions_mappedtoAnnotation_gene.txt")
#genome27_dict = ReadVcfAnnotFiletoDict("/home/thimmamp/MAGIC16/indels_results/IndividualGenomes/AnnotatedVcfs/genome27_XI_adm_MH63_deletions_mappedtoAnnotation_gene.txt")
#genome28_dict = ReadVcfAnnotFiletoDict("/home/thimmamp/MAGIC16/indels_results/IndividualGenomes/AnnotatedVcfs/genome28_XI_1A_ZhenShan97_deletions_mappedtoAnnotation_gene.txt")
#
### get the merged dict to overall view
#from collections import ChainMap
#alldict =[genome1_dict, genome4_dict, genome5_dict,genome6_dict,genome7_dict,genome8_dict,genome9_dict,genome10_dict,genome11_dict,genome12_dict,genome13_dict,genome14_dict,genome15_dict,genome16_dict,genome27_dict,genome28_dict]
#allkey = set().union(*alldict)
##print(str(allkey), '\n')
#with open(sys.argv[1],'w') as of:
#    of.write("Gene\tgenome1_GJ_temp_IRGSP\tgenome4_cA1_N22\tgenome5_GJ_trop1_Azucena\tgenome6_XI_1B1_IR_64\tgenome7_cB_ARC10497\tgenome8_XI_2B_LarhaMugad\tgenome9_XI_3B2_LiuXu\tgenome10_XI_3B1_KhaoYaiGuang\tgenome11_XI_3A_LIMA\tgenome12_cA2_NATELBORO\tgenome13_XI_1B2_PR_106\tgenome14_GJ_trop2_KetanNangka\tgenome15_GJ_subtrp_ChaoMeo\tgenome16_XI_2A_GobolSail\tgenome27_XI_adm_MH63\tgenome28_XI_1A_ZhenShan97\n")
#    for gene in allkey:
#        outstr=gene+'\t'
#        if gene in genome1_dict.keys():
#            outstr+=genome1_dict[gene]+'('+str( len(genome1_dict[gene].split(', ')))+')\t'
#        else:
#            outstr+='NA\t'
#        if gene in genome4_dict.keys():
#            outstr+=genome4_dict[gene]+'('+str( len(genome4_dict[gene].split(', ')) )+')\t'
#        else:
#            outstr+='NA\t'
#        if gene in genome5_dict.keys():
#            outstr+=genome5_dict[gene]+'('+str( len(genome5_dict[gene].split(', ')) )+')\t'
#        else:
#            outstr+='NA\t'
#        if gene in genome6_dict.keys():
#            outstr+=genome6_dict[gene]+'('+str( len(genome6_dict[gene].split(', ')) )+')\t'
#        else:
#            outstr+='NA\t'
#        if gene in genome7_dict.keys():
#            outstr+=genome7_dict[gene]+'('+str( len(genome7_dict[gene].split(', ')) )+')\t'
#        else:
#            outstr+='NA\t'
#        if gene in genome8_dict.keys():
#            outstr+=genome8_dict[gene]+'('+str( len(genome8_dict[gene].split(', ')) )+')\t'
#        else:
#            outstr+='NA\t'
#        if gene in genome9_dict.keys():
#            outstr+=genome9_dict[gene]+'('+str( len(genome9_dict[gene].split(', ')) )+')\t'
#        else:
#            outstr+='NA\t'
#        if gene in genome10_dict.keys():
#            outstr+=genome10_dict[gene]+'('+str( len(genome10_dict[gene].split(', ')))+')\t'
#        else:
#            outstr+='NA\t'
#        if gene in genome11_dict.keys():
#            outstr+=genome11_dict[gene]+'('+str( len(genome11_dict[gene].split(', ')))+')\t'
#        else:
#            outstr+='NA\t'
#        if gene in genome12_dict.keys():
#            outstr+=genome12_dict[gene]+'('+str( len(genome12_dict[gene].split(', ')))+')\t'
#        else:
#            outstr+='NA\t'
#        if gene in genome13_dict.keys():
#            outstr+=genome13_dict[gene]+'('+str( len(genome13_dict[gene].split(', ')))+')\t'
#        else:
#            outstr+='NA\t'
#        if gene in genome14_dict.keys():
#            outstr+=genome14_dict[gene]+'('+str( len(genome14_dict[gene].split(', ')))+')\t'
#        else:
#            outstr+='NA\t'
#        if gene in genome15_dict.keys():
#            outstr+=genome15_dict[gene]+'('+str( len(genome15_dict[gene].split(', ')))+')\t'
#        else:
#            outstr+='NA\t'
#        if gene in genome16_dict.keys():
#            outstr+=genome16_dict[gene]+'('+str( len(genome16_dict[gene].split(', ')))+')\t'
#        else:
#            outstr+='NA\t'
#        if gene in genome27_dict.keys():
#            outstr+=genome27_dict[gene]+'('+str( len(genome27_dict[gene].split(', ')))+')\t'
#        else:
#            outstr+='NA\t'
#        if gene in genome28_dict.keys():
#            outstr+=genome28_dict[gene]+'('+str( len(genome28_dict[gene].split(', ')))+')\n'
#        else:
#            outstr+='NA\n'
#        of.write(outstr)    
#        
