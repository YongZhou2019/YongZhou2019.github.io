import os,sys,csv,random,re,subprocess
from pathlib import Path
from subprocess import Popen, PIPE

## Objective is to annotate indels for each genome on MAGIC16 group with IRGSP(genome1 as ref)
## Usage: python AnnotateMagic16Genomes_indels.py ~/MAGIC16/tmp/vcf/testing/genome1vsMagic16genomes.snpeff.SNPs.vcf2
## usage: python AnnotateMagic16Genomes_indels.py ~/MAGIC16/annotation/genome1_longestisofm_introns_intergenic.gff ~/MAGIC16/tmp/vcf/genome1.filterCount5.genomewide.INDELs.withID.PlinkFormat.vcf

## for randomisation of indels
## python AnnotateMagic16Genomes_indels_test.py ~/MAGIC16/tmp/vcf/genome1.filterCount5.genomewide.INDELs.withID.PlinkFormat_withnocomments.vcf ~/MAGIC16/tmp/vcf/genome6_sample_names_from_header.txt ~/MAGIC16/indels_results/genome1Vs3K_indels_10randomisation_output.txt

## python AnnotateMagic16Genomes_indels_test.py ~/MAGIC16/indels_results/mapman_annotation_forparsing.csv ~/MAGIC16/annotation/IRGSPgramene_IRGSPRAP.one2one ~/MAGIC16/annotation/ReportedGeneFunrice.bed ~/MAGIC16/annotation/Oryza_sativa.IRGSP-1.0.48.gff3 ~/MAGIC16/indels_results/genome1_start_lost_details.txt ~/MAGIC16/indels_results/genome1_start_lost_details_MAPmanandRAPannotation.txt

## python AnnotateMagic16Genomes_indels_dev.py /home/thimmamp/MAGIC16/tmp/vcf/genome1.filterCount5.genomewide.INDELs.withID.PlinkFormat_withnocomments_singleAllele<chr>.vcf ~/MAGIC16/tmp/vcf/genome6_sample_names_from_header.txt

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


def ProcessSubsetSamplesinVcfBiMultiAllele(invcffile, samplesfile, percentage):
    selectedsamples = GetRandomSampleset(samplesfile, percentage)
    inscount=0
    delscount=0
    with open(invcffile, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            permutedlist = list()
            for sample in selectedsamples:
                permutedlist.append(row[sample]) ## store the value of homo/heterozygous for permuted samples in a list
                #if len(row['ALT']) <= 50: ## already filtered for length
                #if len(row['REF']) > len(row['ALT']):
                for val in permutedlist:
                    if val=="0|1" or val=="1|1" or val=="0/1" or val=="1/1":
                       typeval =row['INFO'].split(';')[-1]
                       if ',' in typeval:
                           insdel=typeval.split(',')
                           for item in insdel:
                               if item=="ins":
                                   inscount+=1
                               elif item=="del":
                                   delscount+=1
                       else:
                           if typeval=="ins":
                               inscount+=1
                           elif typeval=="del":
                               delscount+=1
    #print(len(selectedsamples), delscount, inscount)                        
    if int(inscount) > 0 or int(delscount) > 0: 
        return(len(selectedsamples), inscount, delscount)



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
        #for txptkey, components in indict[k].items():
        #    print(k, txptkey, components,'\n')
        print(k, v, '\n')

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

def GetVariationClassification(infile, outfile):
    print("From GetVariationClassification function... Processing ", infile, "...")
    outDict = {}
    with open(infile, mode='r') as f:
         for line in f:
             #print(line)
             #if not line.startswith("##") or not line.startswith("#"): #ignore comments
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


def GetVariantsClassifiedForFeatureCoord(indict, vcfDict, outfile):
    with open(outfile, 'w') as of:
        of.write("featureid\tgene\ttranscript\tfeature_coord\tvariant\tChromosome\tvariant_position\tRef\tAlt\tvcf_details\tINDEL_classification\n")
        for k,v in indict.items():
            featureval = list(v)[0]
            chromosome = featureval.split(':')[0]
            start = int(featureval.split(':')[1].split('-')[0])
            end = int(featureval.split(':')[1].split('-')[1])
            for vcfk,vcfv in vcfDict.items():
                tokens =vcfk.split('_')
                if tokens[0] == chromosome:
                    if int(tokens[1]) >= start and int(tokens[1]) <= end:
                        txpt = k.split('.exon')[0]
                        gene = k.split('_')[0]
                        of.write(k+'\t'+gene+'\t'+txpt+'\t'+featureval+'\t'+vcfk+'\t'+vcfv+'\n')
                        #print(k+'\t'+gene+'\t'+txpt+'\t'+featureval+'\t'+vcfk+'\t'+vcfv+'\n')


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
def GetRAPidforGramene(infile):
    outDict = {}
    with open(infile, 'r') as f:
        for line in f:
            tokens = line.rstrip('\n').split('\t')
            outDict[tokens[1]] = tokens[0]
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

## call functions
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

#for i in range(3):
#    for percentage in [0.1, 0.2]:
#        selectedsamples = GetRandomSampleset(sys.argv[1], percentage)
#        samples = ','.join(selectedsamples)
        #print(str(i), str(percentage), samples)
#        output = subprocess.Popen(['wc -l'], stdin=open(sys.argv[1]), stdout=PIPE)
        #output = Popen(["wc -l"] + sys.argv[1], stdout = PIPE)
#        response = output.communicate()[0]
#        print(response)
        #rv, output = commands.getstatusoutput('wc -l sys.argv[1]')
        #print(output)


#Test Random function
#GetRandomSampleset(sys.argv[1], 0.1)
chromosome=sys.argv[1].split('_')[3]
#print(chromosome)

numrand=100
percentages = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

outprefix="/home/thimmamp/MAGIC16/indels_results/genome6Vs3K_indels_"
outfile=outprefix+str(numrand)+"rand_"+chromosome+"_output_test_orig.txt"
#print(outfile)

with open(outfile,'w') as f:
    for i in range(numrand): #change to 10
        for percentage in percentages:
            numsamples, ins, dels = ProcessSubsetSamplesinVcf(sys.argv[1], sys.argv[2], percentage)
            #numsamples, ins, dels = ProcessSubsetSamplesinVcfSingleallele(sys.argv[1], sys.argv[2], percentage)
            #numsamples, ins, dels = ProcessSubsetSamplesinVcfBiMultiAllele(sys.argv[1], sys.argv[2], percentage)
            #print(percentage, numsamples, ins, dels)
            f.write(str(percentage*10)+'\t'+str(i+1)+'\t'+str(numsamples)+'\t'+chromosome+'\t'+str(ins)+'\t'+str(dels)+'\n')

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
###########
