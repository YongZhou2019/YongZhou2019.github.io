import os,sys,csv,random,re,subprocess
from pathlib import Path
from subprocess import Popen, PIPE
from itertools import islice

## Objective is to annotate indels for each genome on MAGIC16 group with IRGSP(genome1 as ref)
## Usage: python AnnotateMagic16Genomes_indels.py ~/MAGIC16/tmp/vcf/testing/genome1vsMagic16genomes.snpeff.SNPs.vcf2
## usage: python AnnotateMagic16Genomes_indels.py ~/MAGIC16/annotation/genome1_longestisofm_introns_intergenic.gff ~/MAGIC16/tmp/vcf/genome1.filterCount5.genomewide.INDELs.withID.PlinkFormat.vcf

## for randomisation of indels
## python AnnotateMagic16Genomes_indels_test.py ~/MAGIC16/tmp/vcf/genome1.filterCount5.genomewide.INDELs.withID.PlinkFormat_withnocomments.vcf ~/MAGIC16/tmp/vcf/genome6_sample_names_from_header.txt ~/MAGIC16/indels_results/genome1Vs3K_indels_10randomisation_output.txt

## python AnnotateMagic16Genomes_indels_test.py ~/MAGIC16/indels_results/mapman_annotation_forparsing.csv ~/MAGIC16/annotation/IRGSPgramene_IRGSPRAP.one2one ~/MAGIC16/annotation/ReportedGeneFunrice.bed ~/MAGIC16/annotation/Oryza_sativa.IRGSP-1.0.48.gff3 ~/MAGIC16/indels_results/genome1_start_lost_details.txt ~/MAGIC16/indels_results/genome1_start_lost_details_MAPmanandRAPannotation.txt

## For getting dispensable genes, indelcounts and annotation of IRGSP genes
## python AnnotateMagic16Genomes_indels_dev.py ~/MAGIC16/annotation/Nipponbare_IRGSPRAP.one2one ~/MAGIC16/annotation/Oryza_sativa.IRGSP-1.0.48.gff3 ~/MAGIC16/annotation/mapman_annotation_forparsing.csv ~/MAGIC16/indels_results/MAGIC16Orthologs_dispensable_genes_indelcounts.txt > MAGIC16Orthologs_dispensable_genes_indelcounts_withannotation.txt


## for getting indel counts for specific genes
#python AnnotateMagic16Genomes_indels_dev.py ~/MAGIC16/annotation/genome6_genes.gff /home/thimmamp/MAGIC16/tmp/vcf/genome6_processed/genome6.nostar.bimultiallele.50bpindel_annotated_col1_9.vcf ~/MAGIC16/magic16_homologs_20230125/specific/genome6.IR64.gset.specific
## dispensable orthologs header
## MH63RS3	ZS97RS3	117425	125619	125827	127518	127564	127652	127742	128077	132278	132424	Azucena	IR64RS2	Nipponbare	117534
header2genome = {'MH63RS3' : 'genome27',
                 'ZS97RS3' : 'genome28',
                 '117425' : 'genome7',
                 '125619' : 'genome8',
                 '125827' : 'genome9',
                 '127518' : 'genome10',
                 '127564' : 'genome11',
                 '127652' : 'genome12',
                 '127742' : 'genome13',
                 '128077' : 'genome14',
                 '132278' : 'genome15',
                 '132424': 'genome16',
                 'Azucena' : 'genome5',
                 'IR64RS2' : 'genome33',
                 'Nipponbare' : 'genome1',
                 '117534' : 'genome4',
                 }

gff2shortname = {'Osativa.' : 'OsNip_',
                 'Osazucena_' : 'OsAzu_',
                 'Os132424_' : 'OsGoSa_',
                 'Os132278_' : 'OsCMeo_',
                 'Os128077_' : 'OsKeNa_',
                 'Os127742_' : 'OsPr106_',
                 'Os127652_' : 'OsNaBo_',
                 'Os127564_' : 'OsLima_',
                 'Os127518_' : 'OsKYG_',
                 'Os125827_' : 'OsLiXu_',
                 'Os125619_' : 'OsLaMu_',
                 'Os117425_' : 'OsARC_',
                 'Osaus.' : 'OsN22_',
                 'OsIR64_' : 'OsIR64_',
                 'Osmh63.' : 'Osmh63.',
                 'Oszs97.' : 'Oszs97.',
        }

## infile format is
## 13020,OsGoSa_05g0013030,OsGoSa_05g0012710	OsAzu_05g0012740	OsIR64_05g0012520,OsIR64_05g0012500,OsIR64_05g0012820	OsNip_05G013270,OsNip_05G013260	OsN22_05G012430,OsN22_05G012760,OsN22_05G012750	
#Osmh63.12G021880	Oszs97.12G021310,Oszs97.12G021320	OsARC_12g0021710,OsARC_12g0021700	OsLaMu_12g0021840,OsLaMu_12g0021850	OsLiXu_12g0021630	OsKYG_12g0021660,OsKYG_12g0021670	OsLima_12g0021450,OsLima_12g0021440	OsNaBo_12g0021080	OsPr106_12g0021320	OsKeNa_12g0022000,OsKeNa_12g0022010	OsCMeo_12g0022370,OsCMeo_12g0022380	OsGoSa_12g0021550	OsAzu_12g0021920,OsAzu_12g0021930	OsIR64_12g0021110,OsIR64_12g0021120	OsNip_12G022960,OsNip_12G022970	OsN22_12G021340,OsN22_12G021330	
## 15th column has Nipponbare 
## for core analysis we are going to use Nipponbare genes, one entry per ortholog row ie if there are multiple orthologs in Nipponbare, we use anyone of this
def GetCoreGenesForM16Orthologs(infile):
    outDict={}
    with open(infile, 'r') as f:
        next(f) # skip header line
        for line in f:
            token=line.rstrip('\n').split('\t')
            if ',' in token[14]: #get first entry if multiple orthologs
                outDict[token[14].split(',')[0]] = 1
            else:    
                outDict[token[14]]=1
    return outDict        

def GetDispensableGenesForM16Orthologs(infile):
    outDict={}
    header=["MH63RS3","ZS97RS3","117425","125619","125827","127518","127564","127652","127742","128077","132278","132424","Azucena","IR64RS2","Nipponbare","117534"]
    with open(infile, 'r') as f:
        next(f) # skip header line
        for line in f:
            #print(line)
            tokens=line.rstrip('\n').split('\t')
            if tokens[14]!=".": ## First check in genome1 entry for this table
                if ',' in tokens[14]:
                    k = tokens[14].split(',')[0]
                    #print(k, header2genome[header[14]])
                    outDict[k]=header2genome[header[14]]
                else:
                    #print(tokens[14], header2genome[header[14]])
                    outDict[tokens[14]]=header2genome[header[14]]
            else:
                for i in range(16): ## look for genes from first 14 , skip 15th genomes in the table    
                    if i!= 14 and tokens[i]!=".":
                        if ',' in tokens[i]:
                            k = tokens[i].split(',')[0]
                            #print(k, header2genome[header[i]])
                            outDict[k]=header2genome[header[i]]
                            break
                        else:
                            #print(tokens[i], header2genome[header[i]])
                            outDict[tokens[i]]=header2genome[header[i]]
                            break
    return outDict 

def GetFile2Dict(infile):
    outDict = {}
    with open(infile, 'r') as f:
        for line in f:
            token=line.rstrip('\n')
            outDict[token]=1
    return outDict        

def GetNipId2RAPid(infile):
    outDict={}
    with open(infile, 'r') as f:
        for line in f:
            tokens=line.rstrip('\n').split('\t')
            outDict[tokens[0]]=tokens[1]
    return outDict        

def GetAnnotationtoDictfromGff(annotfile, field):
    outDict={}
    with open(annotfile, mode='r') as f:
        for line in f:
            if not line.startswith('#'):
                tokens=line.rstrip('\n').split('\t')
                if tokens[2]==field:
                    ## get genes and their coords
                    temp = tokens[-1].split(';')[0] ## ID
                    keylist = re.split(':|=', temp) ## gene with delimiter either : or =
                    key = keylist[-1]
                    if "UnG" not in key:
                    #key = tokens[-1].split(';')[0].split('=')[1] ## gene
                        val = tokens[0]+':'+tokens[3]+'-'+tokens[4]
                    #key=tokens[0]+':'+tokens[3]+'-'+tokens[4]
                    #outDict[key]=tokens[-1]
                        key1=''
                        for k,v in gff2shortname.items():
                            if k in key:
                                key1=key.replace(k,v)
                        if key1!="":        
                            outDict[key1]=val
    return outDict

def GetAnnotationtoDictfromGffForQuery(annotfile, field, querygenes):
    outDict={}
    with open(annotfile, mode='r') as f:
        for line in f:
            if not line.startswith('#'):
                tokens=line.rstrip('\n').split('\t')
                if tokens[2]==field:
                    ## get genes and their coords
                    temp = tokens[-1].split(';')[0] ## ID
                    keylist = re.split(':|=', temp) ## gene with delimiter either : or =
                    key = keylist[-1]
                    #key = tokens[-1].split(';')[0].split('=')[1] ## gene
                    val = tokens[0]+':'+tokens[3]+'-'+tokens[4]
                    #key=tokens[0]+':'+tokens[3]+'-'+tokens[4]
                    #outDict[key]=tokens[-1]
                    key1=''
                    for k,v in gff2shortname.items():
                        if k in key:
                            temp=key.replace(k,v)
                            if temp in querygenes.keys(): ## coords for only query genes
                                key1=temp
                    if key1!="":        
                        outDict[key1]=val
    return outDict

## input file format
## thimmamp@kw61043:~/MAGIC16/indels_results$ head ~/MAGIC16/annotation/cloned.gene.list
##1	RAP2018-11-26	gene	27143	28644	.	+	.	ID=gene:Os01g0100700;Name=RPS5;biotype=protein_coding;description=Similar to 40S ribosomal protein S5-1;gene_id=Os01g0100700;logic_name=rapdb_genes
##1	RAP2018-11-26	gene	35623	41136	.	+	.	ID=gene:Os01g0100900;Name=SPL1;biotype=protein_coding;description=Sphingosine-1-phosphate lyase%2C Disease resistance respons;gene_id=Os01g0100900;logic_name=rapdb_genes
##1	RAP2018-11-26	gene	82426	84095	.	+	.	ID=gene:Os01g0101700;Name=OsDjC1;biotype=protein_coding;description=Similar to chaperone protein dnaJ 20;gene_id=Os01g0101700;logic_name=rapdb_genes
##1	RAP2018-11-26	gene	89763	91465	.	-	.	ID=gene:Os01g0102000;Name=OsNPC6;biotype=protein_coding;description=Phosphoesterase family protein;gene_id=Os01g0102000;logic_name=rapdb_genes
def GetCloneidForRAPid(clonegff):
    outDict={}
    with open(clonegff, mode='r') as f:
        for line in f:
            if not line.startswith('#'):
                det=line.rstrip('\n').split('\t')[8].split(';')
                rapid=''
                cloneid=''
                description=''
                for item in det:
                    if "ID=gene" in item:
                        rapid=item.split(':')[1]
                    if "Name" in item:
                        cloneid=item.split('=')[1]
                    else:
                        cloneid='NA'
                    if "description=" in item:
                        description=item.split('=')[1]
                    #else:
                    #    description='NA'
                outDict[rapid]=cloneid+';'+description        
    return outDict



## function to read text file with one sample name per line to a python set
## create a random set of samples with specified number
def GetRandomSampleset(infile, percentage):
    samples = set()
    with open(infile, mode='r') as f:
        for line in f:
            token=line.rstrip('\n')
            samples.add(token)
    #print(len(samples), ' :: ' , samples)
    elements = round(len(samples)*float(percentage))
    selection = set(random.sample(samples,elements))
    #print(percentage, str(elements), selection)
    return selection


def GetOrthologsofQueryGenome(querydict, querygenome,orthogenome, orthologfile):
    ortho_dict={}
    with open(orthologfile, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        #colnames = reader.fieldnames
        for row in reader:
            for k in querydict.keys():
                if k==row[querygenome]:
                    if row[orthogenome]!=".":
                        if row[orthogenome] in ortho_dict.keys() and k not in ortho_dict[row[orthogenome]]:
                            ortho_dict[row[orthogenome]].append(k)
                        else:
                            ortho_dict[row[orthogenome]]=[k]
    return ortho_dict                            

## Get number of ins/dels per Mb of Chr of 


## input file format
#GeneName        GeneId  TranscriptId    BioType variants_impact_HIGH    variants_impact_LOW     variants_impact_MODERATE        variants_impact_MODIFIER        variants_effect_3_prime_UTR_truncation  variants_effect_3_prime_UTR_variant     variants_effect_5_prime_UTR_truncation  variants_effect_5_prime_UTR_variant     variants_effect_conservative_inframe_deletion   variants_effect_conservative_inframe_insertion  variants_effect_disruptive_inframe_deletion     variants_effect_disruptive_inframe_insertion    variants_effect_downstream_gene_variant variants_effect_exon_loss_variant       variants_effect_frameshift_variant      variants_effect_intron_variant  variants_effect_non_coding_transcript_variant   variants_effect_splice_acceptor_variant variants_effect_splice_donor_variant    variants_effect_splice_region_variant   variants_effect_start_lost      variants_effect_start_retained_variant  variants_effect_stop_gained     variants_effect_stop_lost       variants_effect_stop_retained_variant   variants_effect_upstream_gene_variant
#Osir64_01g0000010.01.exon.1     GENE_Osir64_01g0000010.01.exon.1        transcript:Osir64_01g0000010.01 protein_coding  0       0       0       96      0       1       0       0       0       0       0       0       33      0       0       21      0       0       0       0       0       0       0       0       0       41
#Osir64_01g0000020.01.exon.1     GENE_Osir64_01g0000020.01.exon.1        transcript:Osir64_01g0000020.01 protein_coding  7       0       5       38      0       0       0       5       0       4       0       1       17      0       7       0       0       0       0       0       0       0       1       0       0       16
#Osir64_01g0000030.01.exon.1     GENE_Osir64_01g0000030.01.exon.1        transcript:Osir64_01g0000030.01 protein_coding  10      0       7       33      0       0       0       0       1       5       1       0       16      0       10      0       0       0       0       0       0       0       0       0       0       17
## function parses this files, looks for variants_impact_HIGH > 1 and non-zero value for each variants_effects and write one file per variants_effect
## Usage : thimmamp@kw61043:~/PythonPrograms/utilities$ python AnnotateMagic16Genomes_indels_dev.py ~/MAGIC16/tmp/vcf/genome6_processed/genome6.nostar.bimultiallele.50bpindel_annotated_col1_9.csv2.genes.txt

def GetVariantEffectwithHighImpact(inputgenesfile, variant_effect_type, genome):

    ofname=genome+"_"+variant_effect_type+"_genes.txt"
    with open(inputgenesfile, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        print("printing variants effect for ", effect)
        with open(ofname, 'w') as of:
            for row in reader:
                #print(row['variants_impact_HIGH'], row[effect])
                if int(row['variants_impact_HIGH']) > 1 and int(row[variant_effect_type]) > 0:
                    #print(row['GeneName'], row['GeneId'], row['TranscriptId'], row['BioType'], row['variants_impact_HIGH'],row[effect])
                    of.write(row['GeneName']+'\t'+row['GeneId']+'\t'+row['TranscriptId']+'\t'+row['BioType']+'\t'+row['variants_impact_HIGH']+'\t'+row[effect]+'\n')
        print("done printing variants effect for ", effect)


def ProcessSubsetSamplesinVcfSingleallele(invcffile,samplesfile, percentage):
    selectedsamples = GetRandomSampleset(samplesfile, percentage)
    inscount=0
    delscount=0
    with open(invcffile, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if ',' in row['ALT']: # ignore multi-allele
                next;
            else:
                permutedlist = list()
                for sample in selectedsamples:
                    permutedlist.append(row[sample]) ## store the value of homo/heterozygous for permuted samples in a list
                if len(row['ALT']) <= 50:
                        if len(row['REF']) > len(row['ALT']):
                                localdelcount=0
                                for val in permutedlist:
                                    if val=="0|1" or val=="1|1" or val=="0/1" or val=="1/1":
                                        localdelcount+=1
                                if localdelcount > 1:
                                    delscount+=1
                        if len(row['REF']) < len(row['ALT']):
                                localinscount=0
                                for val in permutedlist:
                                    if val=="0|1" or val=="1|1" or val=="0/1" or val=="1/1":
                                        localinscount+=1
                                if localinscount > 1:
                                       inscount+=1
    #print(delscount, inscount)                        
    if int(inscount) > 0 or int(delscount) > 0: 
       return(len(selectedsamples), inscount, delscount)

def ProcessSubsetSamplesinVcf(invcffile,samplesfile, percentage):
    selectedsamples = GetRandomSampleset(samplesfile, percentage)
    #print(len(selectedsamples))
    #variantcols = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
    inscount=0
    delscount=0
    with open(invcffile, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if ',*' in row['ALT']:
                next;
            else:    
                #first9tokens = row[variantcols]
                permutedlist = list()
                for sample in selectedsamples:
                    permutedlist.append(row[sample])
                #print(row['#CHROM'],row['POS'],row['ID'],row['REF'],row['ALT'],row['QUAL'],row['FILTER'],row['INFO'],row['FORMAT'],permutedlist)    
                if ',' in row['ALT']: ## indicates multiple INS and DELS
                    temp = row['ALT'].split(',')
                    iteminscount=[]
                    itemdelscount=[]
                    for item in temp:
                        if len(item) <= 50:
                            if len(row['REF']) > len(item):
                                    localdelcount=0
                                    for val in permutedlist:
                                        if val=="0|1" or val=="1|1" or val=="0/1" or val=="1/1":
                                            localdelcount+=1
                                    if localdelcount > 1:
                                       itemdelscount.append(localdelcount)
                            if len(row['REF']) < len(item):
                                    localinscount=0
                                    for val in permutedlist:
                                        if val=="0|1" or val=="1|1" or val=="0/1" or val=="1/1":
                                           localinscount+=1
                                    if localinscount > 1:
                                       iteminscount.append(localinscount)
                    if 0 not in iteminscount: 
                       inscount+=1
                    if 0 not in itemdelscount:
                       delscount+=1
                else:
                    if len(row['ALT']) <= 50:
                        if len(row['REF']) > len(row['ALT']):
                                localdelcount=0
                                for val in permutedlist:
                                    if val=="0|1" or val=="1|1" or val=="0/1" or val=="1/1":
                                        localdelcount+=1
                                if localdelcount > 1:
                                    delscount+=1
                        if len(row['REF']) < len(row['ALT']):
                                localinscount=0
                                for val in permutedlist:
                                    if val=="0|1" or val=="1|1" or val=="0/1" or val=="1/1":
                                        localinscount+=1
                                if localinscount > 1:
                                       inscount+=1
    return(len(selectedsamples), inscount, delscount)                                 



def PrintDict(indict):
    for k, v in indict.items():
        #print(k, " has ", len(v), " number of values")
        print(k, v)

def WriteDict(indict, outfile):
    with open(outfile, 'w') as outf:
        #outf.write('Key\tValues\n')
        outf.write('Geneid\tListofVariants\tNumberofINS\tNumberofDELS\n')
        for k, v in indict.items():
            outf.write(k+'\t'+v+'\n')
            #if type(v) is list:
            #    outf.write(k+'\t'+','.join(v)+'\n')
            #else:    
            #    outf.write(k+'\t'+v+'\n')

def GetGffCoordForFeature(infile, gff_file):
    ## input file with feature id ie gene, mRNA, exon, CDS, etc
    ## output would be feature's coord as ChrNN:start-end from a gff file
    ## output dict with feature id as key and set of coords as key
    featureDict = GetAnnotationtoDictfromGff(gff_file, 'exon')
    outDict={}
    with open(infile, 'r') as inf:
        for line in inf:
            idtokens = line.rstrip('\n').split('\t')
            feature_id = idtokens[0]
            feature_loci=''
            for k,v in featureDict.items():
                if feature_id in v: # check if feature id in the 9th details column of gff
                    feature_loci = k ## get feature loci
            if feature_loci:        
                if feature_id in outDict.keys():
                    if feature_loci not in outDict.keys():
                        outDict[feature_id].add(feature_loci)
                else:
                    outDict[feature_id]={feature_loci}
    return outDict
    #PrintDict(outDict)
            
    
def GetGffFeaturestoDict(gff_file):
    outDict={}
    with open(gff_file, mode='r') as f:
        for line in f:
            if not line.startswith("#"):
                tokens = line.rstrip('\n').split('\t')
                if tokens[2]!="chromosome": 
                    ## value = ID=five_prime_utr-1;Parent=transcript:Osativa.01G000010_01
                    key = tokens[0]+':'+tokens[3]+'-'+tokens[4]+';'+tokens[2]
                    outDict[key] = tokens[-1]
    return outDict               


## This function takes first 9 columns of a vcf file, classfies the variant column as insertion, deletion and SNP
## input is standard vcf file
## output is first 9 columns of vcf and 10th column as variant classified
def WriteVariantClassification(invcfile, outvcfile):
    with open(invcfile, 'r') as f:
        with open(outvcfile, 'w') as f1:
            for line in f:
                if not line.startswith("#"):
                    tokens = line.rstrip('\n').split('\t')[0:9] # get first 9 columns
                    if ',*' in tokens[4]:
                        next;
                    else:
                        #print(line)
                        indel_class=set()
                        ref = tokens[3]
                        if ',' in tokens[4]:
                            temp = tokens[4].split(',')
                            for item in temp:
                                if len(item) <= 50:
                                    #print(line, "\n","mutiple SV, shorter than 50bp")
                                    if len(ref) > len(item):
                                        indel_class.add("deletion")
                                    elif len(ref) < len(item):
                                        indel_class.add("insertion")
                                    elif len(ref) == len(item):
                                        indel_class.add("SNP")
                        else:
                            if len(tokens[4]) <= 50:
                            #print(line, "\n","single SV, shorter than 50bp")
                                if len(ref) > len(tokens[4]):
                                    indel_class.add("deletion")
                                elif len(ref) < len(tokens[4]):
                                    indel_class.add("insertion")
                                elif len(ref) == len(tokens[4]):
                                    indel_class.add("SNP")
                        #f1.write(line.rstrip('\n')+'\t'+','.join(indel_class)+'\n')
                        f1.write('\t'.join(tokens)+'\t'+','.join(indel_class)+'\n')

def GetVariantClassifiction(invcfile):
    # input is a vcf file with 9 columns
    # format
    # #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
    #1       9033    Chr01_9033_T_TAACCCTA   T       TAACCCTA,TAACCCTAAACCCTA        2068.22 .       AC=21,4;AF=0.036,0.00678;AN=590;DP=1454;ExcessHet=0;FS=0;InbreedingCoeff=0.3336;MLEAC=256,52;MLEAF=0.434,0.088;MQ=37.85;QD=28.73;SOR=7.223;TYPE=ins,ins GT
    #1       9068    Chr01_9068_AAACCCTAAACCCTAAACCCTAAACCCTAACCCTAAACCCTAAACCCT_A   AAACCCTAAACCCTAAACCCTAAACCCTAACCCTAAACCCTAAACCCT        A       820.8   .       AC=15;AF=0.003439;AN=4362;BaseQRankSum=-0.083;DP=9413;ExcessHet=0.0528;FS=8.649;InbreedingCoeff=0.3757;MLEAC=22;MLEAF=0.005044;MQ=42.06;MQRankSum=0;QD=9.12;ReadPosRankSum=0.319;SOR=2.888;TYPE=del    GT
    # out dict will have k=Chr01:9033-9034, val=Chr01_9033_T_TAACCCTA;T;TAACCCTA,TAACCCTAAACCCTA;multi_allele;ins,ins
    outDict={}
    with open(invcfile, mode='r') as f:
        for line in f:
            if not line.startswith('#'):
                tokens=line.rstrip("\n").split("\t")[0:9] # in case if longer than 9 columns
                chromosome=tokens[2].split('_')[0]
                variant_position=int(tokens[2].split('_')[1])
                variant_type=''
                k=chromosome+':'+str(variant_position)+'-'+str(variant_position+1)
                if ',' in tokens[4]: 
                    #variant_type='multi_allele'
                    next; ## IGNORE multi_alle
                else:
                    variant_type='bi_allele'
                variant_classification = tokens[7].split(';')[-1].split('=')[1]    
                allelecount = tokens[7].split(';')[0].split('=')[1]
                val=tokens[2]+';'+tokens[3]+';'+tokens[4]+';'+variant_type+';'+variant_classification+';'+str(allelecount)
                outDict[k]=val
    return outDict            

'''
def GetVariationClassification(infile, outfile):
    #print("From GetVariationClassification function... Processing ", infile, "...")
    outDict = {}
    with open(infile, mode='r') as f:
         for line in f:
             if not line.startswith("#"): #ignore comments
                 tokens=line.rstrip("\n").split("\t")[0:9]
                 ## get 3rd column details
                 tmp = tokens[2].split('_')
                 chromosome = tmp[0]
                 variant_pos = tmp[1]
                 ref = tmp[2]
                 variant = tmp[3]
                 indel_class=set()
                 #indel_class=[]
                 if ',*' in tokens[4]:
                     next
                 else:
                     if ',' in tokens[4]:
                       temp = tokens[4].split(',')
                       for item in temp:
                           if len(item) <= 50:
                              #print(line, "\n","mutiple SV, shorter than 50bp")
                              if len(ref) > len(item):
                                  #and "deletion" not in indel_class:
                                indel_class.add("deletion")
                              elif len(ref) < len(item):
                                #and "insertion" not in indel_class:
                                indel_class.add("insertion")
                              elif len(ref) == len(item):
                                #and "SNP" not in indel_class:
                                indel_class.add("SNP")
                     else:
                        if len(tokens[4]) <= 50:
                            #print(line, "\n","single SV, shorter than 50bp")
                            if len(ref) > len(tokens[4]):
                                indel_class.add("deletion")
                            elif len(ref) < len(tokens[4]):
                                indel_class.add("insertion")
                            elif len(ref) == len(tokens[4]):
                                indel_class.add("SNP")
                 ## write the class to a file
                 #outDict[tokens[2]]=tokens[0]+'\t'+tokens[1]+'\t'+tokens[3]+'\t'+tokens[4]+'\t'+tokens[7]+'\t'+','.join(indel_class)
                     if len(indel_class)!=0:
                        outDict[tokens[2]]=chromosome+'\t'+variant_pos+'\t'+ref+'\t'+variant+'\t'+tokens[7]+'\t'+','.join(indel_class)
    return outDict


def GetVariantsClassifiedForFeatureCoord(featuredict, annotdict, vcfDict, outfile):
    with open(outfile, 'w') as of:
        of.write("featureid\tgene\ttranscript\tfeature_coord\tvariant\tChromosome\tvariant_position\tRef\tAlt\tvcf_details\tINDEL_classification\n")
        for k, v in featuredict.items():
            for k1,v1 in annotdict.items():
                if k in v1:
                    featurechr = k1.split(':')[0]
                    start = int(k1.split(':')[1].split('-')[0])
                    end = int(k1.split(':')[1].split('-')[1])
                    parent=v1.split(';')[0]
                    txpt=parent.split(':')[1]
                    gene=txpt.split('.')[0]
                    for vcfk, vcfv in vcfDict.items():
                        variantchr = vcfk.split(':')[0]
                        variantstr = int( vcfk.split(':')[1].split('-')[0] )
                        variantend = int( vcfk.split(':')[1].split('-')[1] )
                        if variantchr==featurechr and variantstr >= start and variantend <= end:
                            varianttokens=vcfv.split(';')
                            #above tokens have variantid;REF;ALT;alleletype;indel_classfn
                            of.write(k+'\t'+gene+'\t'+txpt+'\t'+k1+'\t'+varianttokens[0]+'\t'+featurechr+'\t'+varianttokens[0].split('_')[1]+'\t'+
                                    varianttokens[1]+'\t'+varianttokens[2]+'\t'+varianttokens[3]+'\t'+varianttokens[-1]+'\n')
                            #print(k+'\t'+gene+'\t'+txpt+'\t'+k1+'\t'+varianttokens[0]+'\t'+featurechr+'\t'+varianttokens[0].split('_')[1]+'\t'+
                            #        varianttokens[1]+'\t'+varianttokens[2]+'\t'+varianttokens[3]+'\t'+varianttokens[-1])


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
'''
def writeDict2File(indict, outfilename):
    with open(outfilename, 'w') as of:
        for k,v in indict.items():
            #of.write(k+'\t'+str(len(v))+'\t'+str(v)+'\n')
            of.write(k+'\t'+str(v)+'\n')

def ReadVcfAnnotFiletoDict(infile):
    outDict = {}
    with open(infile,'r') as f:
        for line in f:
            tokens = line.rstrip('\n').split('\t')
            tokens[-1] = tokens[-1].replace('{','')
            tokens[-1] = tokens[-1].replace('}','')
            outDict[tokens[0]] = tokens[-1]
    return outDict


def GetMapmanAnnot2Dict(annotfile):
    outDict = {}
    with open(annotfile, 'r') as f:
        for line in f:
            if not line.startswith("BINCODE"):
                tokens = line.rstrip('\n').split('\t')
                if tokens[2]!="":
                    function=''
                    cellcomponent=''
                    tokens[2]= tokens[2].replace('o','O')
                    tokens[2]= tokens[2].replace('g','G')
                    gene = tokens[2]
                    temp = tokens[3].split(':')
                    funcs= temp[3].split(' & ')
                    if ',' in funcs[0]:
                        function = funcs[0].split(',')[0]
                        cellcomponent = funcs[0].split(',')[1]
                    else:
                        function = funcs[0]
                        cellcomponent = 'NA'
                    outDict[gene]=function+';'+cellcomponent
    #PrintDict(outDict)            
    return outDict           

### Add RAP symbol and description to variant annoation
## input fileformat
##thimmamp@kw61043:~/MAGIC16/annotation$ head IRGSPgramene_IRGSPRAP.one2one
## Gramene_id   RAP_id  col3    Chromosome
##Osativa.01G000010	Os01g0100100	RBH	ChrN
##Osativa.01G000020	Os01g0100200	RBH	ChrN
##Osativa.01G000030	Os01g0100300	RBH	ChrN
##Osativa.01G000040	Os01g0100400	RBH	ChrN
##Osativa.01G000050	Os01g0100500	RBH	ChrN
def GetRAPidforGramene(infile):
    outDict = {}
    with open(infile, 'r') as f:
        for line in f:
            tokens = line.rstrip('\n').split('\t')
            #outDict[tokens[1]] = tokens[0]
            outDict[tokens[0]] = tokens[1] ## gramene id is unique hence key, rapid repeated so value.
    return outDict

def GetGeneFunction(infile):
    outDict = {}
    with open(infile, 'r') as f:
        for line in f:
            tokens = line.rstrip('\n').split('\t')
            outDict[tokens[1]] = tokens[5]
    return outDict        

# Get description from gff3 for gramene
## 1       RAP2018-11-26   gene    12808   13978   .       -       .       ID=gene:Os01g0100466;biotype=protein_coding;description=Hypothetical protein;gene_id=Os01g0100466;logic_name=rapdb_genes
#1       RAP2018-11-26   gene    16399   20144   .       +       .       ID=gene:Os01g0100500;biotype=protein_coding;description=Immunoglobulin-like domain containing protein;gene_id=Os01g0100500;logic_name=rapdb_genes

def GetGeneDescription(gff_file):
    outDict ={}
    with open(gff_file, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                tokens=line.rstrip("\n").split("\t")
                if tokens[2]=="gene":
                    gene=''
                    tmp = tokens[-1].split(';')
                    for item in tmp:
                        if "ID=gene" in item:
                            gene = item.split(':')[1]
                        if "description" in item:
                            outDict[gene] = item.split('=')[1]
    return outDict                        

def GetGeneSymbol(inbed):
    outDict = {}
    with open(inbed, 'r') as f:
        for line in f:
            #print(line)
            if not "Start" in line:
                tokens=line.rstrip('\n').split('\t')
                outDict[tokens[3]] = tokens[5]
    return outDict        


def AnnotateVariantUsingGff(variantdict, annotdict):
    featureDetailsDict={}
    featureDict = {}
    outDict = {}
    for k,v in variantdict.items():
        chromosome = k.split('_')[0]
        variantcoord = int(k.split('_')[1])
        variantclass = v.split('\t')[-1]
        for ak,av in annotdict.items():
            achr = ak.split(';')[0].split(':')[0]
            astart = int(ak.split(';')[0].split(':')[1].split('-')[0])
            aend = int(ak.split(';')[0].split(':')[1].split('-')[1])
            feature = ak.split(';')[1]
            if chromosome == achr and astart <= variantcoord and aend >= variantcoord:
               #print("Variant: ",k, "Annot: ", ak)
               if k in outDict.keys() and feature not in featureDict[k]:
                   featureDict[k].append(feature)
                   featureDetailsDict[k].append(ak+'_'+av)
               else:
                   featureDict[k] = [feature]
                   featureDetailsDict[k] = [ak+'_'+av]
                   outDict[k]=[feature+'_'+av]
        print(k, variantclass, ','.join(featureDict[k]), ','.join(featureDetailsDict[k]))             
        #outDict[k]=variantclass+'\t'+','.join(featureDict[k])+'\t'+','.join(featureDetailsDict[k])
    return outDict

def AddMapmanAnnot2File(infile, annotdict, grameneannot, outfile):
    with open(infile, 'r') as f:
        with open(outfile, 'w') as f1:
            for line in f:
                line= line.rstrip('\n')
                if line.startswith('featureid'):
                    f1.write(line+'\tFunction\tCellular_Component\tRAPGeneSymbol\tGeneDescription\n')
                else:
                    annotstr=''
                    gene = line.split('\t')[1]
                    if gene in annotdict.keys():
                        annotstr='\t'+'\t'.join(annotdict[gene].split(';'))
                        #line=line+'\t'+'\t'.join(annotdict[gene].split(';'))+'\n'
                        #f1.write(line)
                    else:
                        #line=line+'\tNA\tNA\n'
                        annotstr='\tNA\tNA'
                        #f1.write(line)
                    if gene in grameneannot.keys():
                        annotstr+='\t'+grameneannot[gene]
                    else:
                        annotstr+='\tNA\tNA'

                    #f1.write(line+'\n')
                    line+=annotstr
                    #print(line,'\n')
                    f1.write(line+'\n')
    print("Done adding annotation...Check output file\n")                    

## Add RAPid and clone id to RAPandMapmann annotated file
def AddRapCloneid(infile, rapdict, gramenedict, outfile):
    with open(infile, 'r') as f:
        with open(outfile, 'w') as f1:
            for line in f:
                line= line.rstrip('\n')
                if line.startswith('featureid'):
                    f1.write(line+'\tRAPid\tCloneid\n')
                else:
                    gene = line.split('\t')[1]
                    rapid='NA'
                    cloneid='NA'
                    if gene in gramenedict.keys():
                       rapid=gramenedict[gene]
                       if rapid in rapdict.keys():
                          cloneid=rapdict[rapid]
                    print(line+'\t'+rapid+'\t'+cloneid)       
                    f1.write(line+'\t'+rapid+'\t'+cloneid+'\n')
    print("Done adding RAP and Clone ids...Check output file\n")                    
                    

#def GetIndelCountsForGenes(gffdict, vcfdict, queryfile):
def GetIndelCountsForGenes(gffdict, vcfdict):
    outDict = {}
    totalgeneicins=0
    totalgeneicdels=0
    #for k, v in queryfile.items():
    #with open(queryfile, 'r') as f:
        #for line in f:
            #line = line.rstrip("\n")
    for k,v in gffdict.items():
        l = gffdict[k] # gene as key, loci as value 'chr:str-end'
        geneinscount=0
        genedelcount=0
        genevariantdetails=[]
        gc = l.split(":")[0]
        gstr = int(l.split(":")[1].split("-")[0])
        gend = int(l.split(":")[1].split("-")[1])
        for v, vl in vcfdict.items():
            vchr = v.split(":")[0]
            vstr = int(v.split(":")[1].split("-")[0])
            vend = int(v.split(":")[1].split("-")[1])
            if gc==vchr and vstr >= gstr and vend <= gend:
                genevariantdetails.append(vl)
                varclass = vl.split(";")[-2]
                if "," in varclass: # multiple ins/dels ## IGNORE FOR NOW
                    next;
                    #temp=varclass.split(",")
                    #for item in temp:
                    #    if item=="ins":
                    #        geneinscount+=1
                    #        totalgeneicins+=1
                    #    elif item=="del":
                    #        genedelcount+=1
                    #        totalgeneicdels+=1
                else:
                    if varclass=="ins":
                        geneinscount+=1
                        totalgeneicins+=1
                    elif varclass=="del":
                        genedelcount+=1
                        totalgeneicdels+=1
                #print(genevariantdetails, geneinscount, genedelcount)                
        if len(genevariantdetails)>=1:  
            outDict[k]=";;".join(genevariantdetails)+'\t'+str(geneinscount)+'\t'+str(genedelcount)
        else:    
            outDict[k]='NA\t'+str(geneinscount)+'\t'+str(genedelcount)
    print("total number of geneic INS: ", str(totalgeneicins),"Total number of geneic DELS: ", str(totalgeneicdels))        
    return outDict

def GetGenomeSpecificDispensableGenes(masterdict, genome):
    outDict={}
    for k,v in masterdict.items():
        #print(k,v)
        if v==genome:
            if v in outDict.keys():
                outDict[v].append(k)
            else:
                outDict[v]=[k]
    return outDict             

## call functions

## To get genomwise total genic INS and dels
#for g in ["genome1", "genome4", "genome5", "genome7", "genome8", "genome9", "genome10", "genome11", "genome12", "genome13", "genome14", "genome15", "genome16", "genome27", "genome28", "genome33"]:
#g = sys.argv[1]
#print("Working with ", g)
#for g in ["genome4"]:
#annotfile = "/home/thimmamp/MAGIC16/annotation/"+g+"/genes.gff"
#vcfile = "/home/thimmamp/MAGIC16/tmp/vcf/"+g+"_processed/"+g+".nostar.bimultiallele.50bpindel_annotated_col1_9.vcf"
#outfile = "/home/thimmamp/MAGIC16/indels_results/"+g+"_genic_InDelsStats_Details.txt"
#featuredict = GetAnnotationtoDictfromGff(annotfile, 'gene')
    #print(dict(islice(featuredict.items(), 0,4)))
#vcfdict = GetVariantClassifiction(vcfile)
    #print(dict(islice(vcfdict.items(), 0,4)))
    #PrintDict(vcfdict)
#feature_indels = GetIndelCountsForGenes(featuredict, vcfdict)
    #print(dict(islice(feature_indels.items(), 0,4)))
    #PrintDict(feature_indels)
#WriteDict(feature_indels, outfile)
#print(" Done Working with ", g)


################################

## Add description, cloneid, mapman annotation to MAGIC16 core, orthologs and specific
nip2rap = GetNipId2RAPid(sys.argv[1]) ## Nipponbare_IRGSPRAP.one2one file as input
rap2clone = GetCloneidForRAPid(sys.argv[2]) ## RAPid 2 cloneid and description sep by ';'
mapmandict = GetMapmanAnnot2Dict(sys.argv[3]) ## Osativa.01GXX to celllar function details from mampman file
## Now read input ortholog file and add annotation
with open(sys.argv[4], 'r') as f:
    for line in f:
        #if not line.startswith("Geneid"):
        #print(line)
        tokens=line.rstrip('\n').split('\t')
        temp=tokens[0].split(' ')
        nipid=temp[0]
        #print(temp[0], temp[1], tokens[1], tokens[2])

        if 'OsNip' in nipid:
            rapid=''
            outstr=''
            #print("NIP id: ", temp[0], "Variant details: ", temp[1], "NumofINS: ",tokens[1], "NumofDELS: ", tokens[2])
            if nipid in nip2rap.keys():
                rapid=nip2rap[nipid]
                outstr+=rapid+'\t'
            else:
                outstr+='NA\t'
            if rapid!='':
                if rapid in rap2clone.keys():
                    clone,desc=rap2clone[rapid].split(';')
                    outstr+=clone+'\t'+desc+'\t'
                else:
                    outstr+='NA\tNA\t'
            else:
                outstr+='NA\tNA\t'
            nipchangeid=nipid.replace('OsNip_', 'Osativa.')
            if nipchangeid in mapmandict:
                outstr+=mapmandict[nipchangeid]
            else:
                outstr+='NA'
            print(nipid+'\t'+ temp[1]+'\t'+ tokens[1]+ '\t'+ tokens[2]+'\t'+outstr)
        else:
            print(temp[0]+'\t'+temp[1]+'\t'+tokens[1]+'\t'+tokens[2]+'\tNA\tNA\tNA\tNA')


'''
########################
## To get indel counts for core genes in MAGIC16 orthologs
### To get count of INS and DELS for each MAGIC16 genomes' genes and their vcf
##featuredict = GetAnnotationtoDictfromGff(sys.argv[1], 'gene')
#PrintDict(featuredict)
#print("Number of genes in the genome :", len(featuredict.keys()))

##vcfdict = GetVariantClassifiction(sys.argv[2]) # input annotated vcf file
#print("Number of INDEL variants in the genome :", len(vcfdict.keys()))
#PrintDict(vcfdict)

## get core/dispensable genes dict
##querygenesdict = GetCoreGenesForM16Orthologs(sys.argv[3])
## get indel counts for each gene in the querygenesdict
##gene2svcount=GetIndelCountsForGenes(featuredict, vcfdict, querygenesdict) # genes in a dict
#gene2svcount=GetIndelCountsForGenes(featuredict, vcfdict, sys.argv[3]) # genes in a file
#print("Number of genes with variants and indel counts for specific genes in the genome :", len(gene2svcount.keys()))
##PrintDict(gene2svcount)
#########################


##### To get indelcounts for dispensable ortho genes
## step1 get dispensable genes and their genomes
dispgenes2genome=GetDispensableGenesForM16Orthologs(sys.argv[1])
print("Done reading Dispensable orthologs...")
print("Number of Orthologous genomes: ", len(dispgenes2genome.keys()))
#PrintDict(dispgenes2genome)

## step2 gather genes from same genome
#g = "genome1"
for g in ["genome1",  "genome5", "genome7", "genome8", "genome9", "genome10", "genome11", "genome12", "genome13", "genome14", "genome15", "genome16", "genome27", "genome28", "genome33"]:
    genome_dispensable=GetGenomeSpecificDispensableGenes(dispgenes2genome, g)
    print("Number of dispensable genes from genome ", g, " are ", len(genome_dispensable[g]))
    #PrintDict(genome_dispensable)
    genes_dict={}
    for k, v in genome_dispensable.items():
        genes=v
        for gene in genes:
            genes_dict[gene]=1
    
    ## step3 get gene loci and indels within them for each of the dispensable genome's genes
    annotfile = "/home/thimmamp/MAGIC16/annotation/"+g+"/genes.gff"
    print("reading annotation file ", annotfile)
    ## step3 get gene loci and indels within them for each of the dispensable genome's genes
    featuredict = GetAnnotationtoDictfromGffForQuery(annotfile, 'gene', genes_dict)
    print("Number of genes from gff reading: ", len(featuredict.keys()))
    #PrintDict(featuredict)
    vcfile = "/home/thimmamp/MAGIC16/tmp/vcf/"+g+"_processed/"+g+".nostar.bimultiallele.50bpindel_annotated_col1_9.vcf"
    #print("reading vcf file ", vcfile)
    vcfdict = GetVariantClassifiction(vcfile) # input annotated vcf file
    #PrintDict(vcfdict)
    ## get indel counts for each gene in the querygenesdict
    #gene2svcount=GetIndelCountsForGenes(featuredict, vcfdict, genes_dict) # genes in a dict
    gene2svcount=GetIndelCountsForGenes(featuredict, vcfdict) # genes in a dict
    outfile = "/home/thimmamp/MAGIC16/indels_results/"+g+"_Dispensable_genic_InDelsStats_Details_testing.txt"
    print("Number of genes from genes2svcount reading: ", len(gene2svcount.keys()))
    #PrintDict(gene2svcount)
    WriteDict(gene2svcount, outfile)
    print("Done writing to file ", outfile, "...")

#outdir="/home/thimmamp/MAGIC16/indels_results/"
#fprefix = os.path.basename(sys.argv[3])
#fname=outdir+fprefix+"_orthologs_indels_count_onlybiallele.txt"
#with open(fname, 'w') as f:
#    f.write("Gene\tVariantDetails\tNumofINS\tNumofDEL\n")
#    for k,v in gene2svcount.items():
#        f.write(k+'\t'+v+'\n')

## get feature file to dict
#features = GetFile2Dict(sys.argv[1])
## Get feature loci and details from gff
#featuredict = GetAnnotationtoDictfromGff(sys.argv[2], 'exon') # input gff file 
#variantid2details = GetVariantClassifiction(sys.argv[3]) # input annotated vcf file
## Call variant and feature annotation, and write to file
#GetVariantsClassifiedForFeatureCoord(features, featuredict, variantid2details, sys.argv[4])

### Get orthologs
#querygenomedict = GetFile2Dict(sys.argv[1])
#query2targetorthodict = GetOrthologsofQueryGenome(querygenomedict, 'genome6', 'genome1', sys.argv[2])
#WriteDict(query2targetorthodict, sys.argv[3])


## append MAPManRAPannot with rapid and cloneid
## Usage: thimmamp@kw61043:~/MAGIC16/indels_results$ python ~/PythonPrograms/utilities/AnnotateMagic16Genomes_indels_dev.py ~/MAGIC16/annotation/cloned.gene.list ~/MAGIC16/annotation/IRGSPgramene_IRGSPRAP.one2one genome1_start_lost_details_MAPmanandRAPannotation.txt genome1_start_lost_details_MAPmanandRAPannotation_withrap_cloneid.txt
#rapid2cloneid = GetCloneidForRAPid(sys.argv[1])
#gramene2rap = GetRAPidforGramene(sys.argv[2])
#AddRapCloneid(sys.argv[3], rapid2cloneid, gramene2rap, sys.argv[4])


## working with snpeff output genes file for variant effect
#variant_effects = ['variants_effect_3_prime_UTR_truncation', 'variants_effect_3_prime_UTR_variant', 'variants_effect_5_prime_UTR_truncation', 'variants_effect_5_prime_UTR_variant', 'variants_effect_conservative_inframe_deletion', 'variants_effect_conservative_inframe_insertion', 'variants_effect_disruptive_inframe_deletion', 'variants_effect_disruptive_inframe_insertion', 'variants_effect_downstream_gene_variant', 'variants_effect_exon_loss_variant', 'variants_effect_frameshift_variant', 'variants_effect_intron_variant', 'variants_effect_non_coding_transcript_variant', 'variants_effect_splice_acceptor_variant', 'variants_effect_splice_donor_variant', 'variants_effect_splice_region_variant', 'variants_effect_start_lost', 'variants_effect_start_retained_variant', 'variants_effect_stop_gained', 'variants_effect_stop_lost', 'variants_effect_stop_retained_variant', 'variants_effect_upstream_gene_variant']

#for effect in variant_effects:
#    genome=sys.argv[1].split('\.')[0]
#    GetVariantEffectwithHighImpact(sys.argv[1], effect, genome)

#PrintDict(gramene2rap)
#annotDict = GetGffFeaturestoDict(sys.argv[1])
#vcfDict = GetVariationClassification(sys.argv[2])
#PrintDict(vcfDict)
#annotatedVariantDict = AnnotateVariantUsingGff(vcfDict, annotDict)
#writeDict2File(annotatedVariantDict, sys.argv[3])

## first file, genome1_stop_lost.txt & second file gff file
## Usage: thimmamp@kw61043:~/PythonPrograms/utilities$ python AnnotateMagic16Genomes_indels_test.py ~/MAGIC16/indels_results/genome1_stop_lost.txt ~/MAGIC16/annotation/genome1_longestisofm_introns_intergenic.gff
## Usage with output file:
## python AnnotateMagic16Genomes_indels_test.py ~/MAGIC16/indels_results/genome1_stop_gained.txt ~/MAGIC16/annotation/genome1_longestisofm_introns_intergenic.gff ~/MAGIC16/tmp/vcf/genome1.filterCount5.genomewide.INDELs.withID.PlinkFormat.col1_9.vcf genome1_stop_gained_details.txt

#feature2lociDict = GetGffCoordForFeature(sys.argv[1], sys.argv[2])
#vcfDict = GetVariationClassification(sys.argv[3])
#GetVariantsClassifiedForFeatureCoord(feature2lociDict, vcfDict, sys.argv[4])

#WriteVariantClassification(sys.argv[1], sys.argv[2])

## Usage for adding mapman dict to indel classified snpeff output
## python AnnotateMagic16Genomes_indels_test.py ~/MAGIC16/indels_results/mapman_annotation_forparsing.csv ~/MAGIC16/indels_results/genome1_start_lost_details.txt ~/MAGIC16/indels_results/genome1_start_lost_details_function_CC_test.txt

## Write variant classfified vcf file
#WriteVariantClassification(sys.argv[1], sys.argv[2])

## Get number of INS and DELS for random subset of vcf samples using bcftools and vcftools
## bcftools view -s ERS467753,ERS467754,ERS467755 genome1.filterCount5.genomewide.INDELs.withID.PlinkFormat.vcf | bcftools view -i '(ILEN >= -50 && ILEN <= 50)' -| vcf-annotate --fill-type  | grep -oP "TYPE=\w+" | sort | uniq -c | grep del | cut -d' ' -f1
#1558785
##(base) thimmamp@kw61043:~/MAGIC16/tmp/vcf$ bcftools view -s ERS467753,ERS467754,ERS467755 genome1.filterCount5.genomewide.INDELs.withID.PlinkFormat.vcf | bcftools view -i '(ILEN >= -50 && ILEN <= 50)' -| vcf-annotate --fill-type  | grep -oP "TYPE=\w+" | sort | uniq -c | grep ins | cut -d' ' -f1
#1073249
## random shuffling of a file
## shuf -n 300 genome6_sample_names_from_header.txt


#Test Random function
#GetRandomSampleset(sys.argv[1], 0.1)
#with open(sys.argv[3],'w') as f:
#    for i in range(100): #change to 10
#        for percentage in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
#            #numsamples, ins, dels = ProcessSubsetSamplesinVcf(sys.argv[1], sys.argv[2], percentage)
#            numsamples, ins, dels = ProcessSubsetSamplesinVcfSingleallele(sys.argv[1], sys.argv[2], percentage)
            #print(percentage, numsamples, ins, dels)
#            f.write(str(percentage*10)+'\t'+str(i+1)+'\t'+str(numsamples)+'\t'+str(ins)+'\t'+str(dels)+'\n')

##########
##mapmandict = GetMapmanAnnot2Dict(sys.argv[1])
#RAP geneid to Gramene
##rap2gramene = GetRAPidforGramene(sys.argv[2])
# RAP geneid to Symbol
##rap2genesymbol = GetGeneSymbol(sys.argv[3])
##genedesc = GetGeneDescription(sys.argv[4])
##GramenetoRAP = dict()
##for k, v in rap2gramene.items():
##    symboldescri=''
##    if k in rap2genesymbol.keys():
##        symboldescri+=rap2genesymbol[k]+'\t'
##    else:
##        symboldescri+='NA\t'
##    if k in genedesc.keys():
##        symboldescri+=genedesc[k]
##    else:
##        symboldescri+='NA'
##    GramenetoRAP[v] = symboldescri    
    #print(k, v, symboldescri,'\n')

##AddMapmanAnnot2File(sys.argv[5], mapmandict, GramenetoRAP, sys.argv[6])
'''
