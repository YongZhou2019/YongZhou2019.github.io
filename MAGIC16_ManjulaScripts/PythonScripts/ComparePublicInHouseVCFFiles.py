import os,sys,csv,random,re,gzip
from pathlib import Path

## Usage : thimmamp@kw61043:~/PythonPrograms/utilities$ python ClassifyVar2indels.py ~/MAGIC16/tmp/
## Usage : python ClassifyVar2indels.py /home/thimmamp/MAGIC16/tmp/vcf/testing/ > genome1VsMagic16_lessthan50bp_deletions.txtfor getting data for cmplot to look insertion/deletion density
## the input file was already filtered for indels < 50bp and have coverage in atleast 1/16 accessions.

## usage:
## python ~/PythonPrograms/utilities/ComparePublicInHouseVCFFiles.py ../public/SAP.raw.vcf.gz SAP_only_SNP.vcf 

## Original genome1 vs 3K samples SNP file ~/MAGIC16/tmp/vcf/genome1.filterCount5.genomewide.INDELs.withID.PlinkFormat.vcf
## To get individual genomes' vcf modified, filtered for length <= 50 and allele count
## python ClassifyVar2indels_Samplewise.py ~/MAGIC16/tmp/notes.txt ~/MAGIC16/tmp/vcf/genome1.filterCount5.genomewide.INDELs.withID.PlinkFormat.vcf

def NcbidtoChr(infile):
    outDict = {}
    with open(infile, mode='r') as f:
        for line in f:
            line=line.rstrip('\n')
            tokens=line.split("\t")
            tokens[0]=tokens[0].rstrip()
            outDict[tokens[0]]=tokens[1]
    print(outDict)            
    return outDict           

def PidtoSampleid(infile):
    outDict = {}
    with open(infile, mode='r') as f:
        for line in f:
            line=line.rstrip('\n')
            tokens=line.split("\t")
            outDict[tokens[-1]]=tokens[0]
    print(outDict)            
    return outDict           

    
def GetSampleDict(infile):
    outDict = {}
    with open(infile, mode='r') as f:
        for line in f:
            if line.startswith("genome"):
                line=line.rstrip('\n')
                tokens=line.split("\t")
                k = tokens[-1]
                v = tokens[0]+'_'+tokens[1]
                outDict[k]=v
    print(outDict)            
    return outDict           


def WriteOnlySNPfromVCF(infile, outfilename):
    with gzip.open(infile, 'rb') as f:
        with open(outfilename, 'w') as of:
            for line in f:
                line = line.decode('utf-8')
                #print(line)
                if line.startswith("#CHROM"):
                    of.write(line)
                elif not line.startswith("#"):
                    tokens=line.split('\t', 6)
                    if len(tokens[3])==1 and len(tokens[4])==1: ## filter only SNPs
                        of.write(line)
                        





def FilterVCFforMinCountAllele(infile):
    with open(infile, mode='r') as f:
         for line in f:
             line = line.rstrip('\n')
             allele_count=line.split("\t")[9:] # look for atleast one 0|0 in the last 16 columns
             count=0
             for item in allele_count:
                 #print(item)
                 temp=[]
                 if "|" in item:
                    temp=item.split("|")
                 elif "/" in item:
                    temp=item.split("/")
                 #temp=filter(None, re.split(""[|\]", item))
                 #temp=item.split("[|/]")
                 if temp[0]!="." and temp[1]!=".":
                    if int(temp[0])>=1 and int(temp[1])>=1:
                       #print("Allele count valid: ",line)
                       count+=1
             if count>=1:
               print(line) 


def GetVariationClassification(infile):
    print("From GetVariationClassification function... Processing ", infile, "...")
    outDict = {}
    with open(infile, mode='r') as f:
         for line in f:
             #print(line)
             if not line.startswith("##") or not line.startswith("#"): #ignore comments
                 tokens=line.rstrip("\n").split("\t")[0:5]
                 #print(tokens)
                 indel_class=[]
                 if ','in tokens[-1]:
                     temp = tokens[-1].split(',')
                     for item in temp:
                         if len(item) <= 50:
                            #print(line, "\n","mutiple SV, shorter than 50bp")
                            if len(tokens[-2]) > len(item) and "deletion" not in indel_class:
                               indel_class.append("deletion")
                            elif len(tokens[-2]) < len(item) and "insertion" not in indel_class:
                               indel_class.append("insertion")
                            elif len(tokens[-2]) == len(item) and "SNP" not in indel_class:
                               indel_class.append("SNP")
                 else:
                     if len(tokens[-1]) <= 50: 
                        #print(line, "\n","single SV, shorter than 50bp")
                        if len(tokens[-2]) > len(tokens[-1]):
                           indel_class.append("deletion")
                        elif len(tokens[-2]) < len(tokens[-1]):
                           indel_class.append("insertion")
                        elif len(tokens[-2]) == len(tokens[-1]):
                           indel_class.append("SNP")
                 #print(tokens[2], tokens[0]+':'+tokens[1]+':'+tokens[3]+':'+tokens[4]+':'+','.join(indel_class))         
                 outDict[tokens[2]]=tokens[0]+'\t'+tokens[1]+'\t'+tokens[3]+'\t'+tokens[4]+'\t'+','.join(indel_class)
    return outDict


def GetVariationStats(inputdict, samplename):
    #print("Processing summary stats for", samplename, "...")
    inscount = 0
    delscount = 0
    bothcount = 0
    otherscount = 0
    #ins_file=samplename+"modified.filtered.insertions.only.vcf"
    #del_file=samplename+"modified.filtered.deletions.only.vcf"
    ins_file=samplename+"modified.filtered.insertions.only.forcmplot.txt"
    del_file=samplename+"modified.filtered.deletions.only.forcmplot.txt"
    with open(ins_file, "w") as ins:
        with open(del_file, "w") as dels:
            for k,v in inputdict.items():
                #print(k,'\t',v)
                svclass = v.split('\t')[-1]
                if "," in svclass:
                    #print(k, "insertion_and_deletion")
                    bothcount+=1
                elif svclass == "insertion":
                    chromosome, pos, other = inputdict[k].split('\t', 2) 
                    ins.write(k+'\t'+chromosome+'\t'+pos+'\n')
                    #print(k,'\t',v)
                    #print(k, "insertion only")
                    inscount+=1
                elif svclass == "deletion":
                    chromosome, pos, other = inputdict[k].split('\t', 2) 
                    dels.write(k+'\t'+chromosome+'\t'+pos+'\n')
                    #print(k,'\t',v)
                    #print(k, "deletion only")
                    delscount+=1
                else:
                    #print(k, "others")
                    otherscount+=1
 
    ## print Summary
    print(samplename+'\t'+str(inscount)+'\t'+str(delscount)+'\t'+str(bothcount)+'\t'+str(otherscount))


## call functions
## python ClassifyVar2indels_Samplewise.py ~/MAGIC16/tmp/notes.txt ~/MAGIC16/tmp/vcf/genome1.filterCount5.genomewide.INDELs.withID.PlinkFormat.vcf
#samplesdict = GetSampleDict(sys.argv[1])
#WriteSampleSpecificVcf(sys.argv[2], samplesdict)


WriteOnlySNPfromVCF(sys.argv[1], sys.argv[2])
#ncbi2chr = NcbidtoChr(sys.argv[1])
#pid2ers = PidtoSampleid(sys.argv[2])



## commented for testing
##print("SampleName\tInsertions\tDeletions\tBoth\tOthers")
##for fn in os.listdir(sys.argv[1]):
##    if fn.endswith(".vcf"):
##        print(os.path.join(sys.argv[1], fn))
##        samplename = '_'.join(fn.split('.')[0:2])
##        #samplename = fn.split('.')[0]
##        #print("Sample Name: ", samplename,'\n')
##        my_dict = GetVariationClassification(os.path.join(sys.argv[1], fn))
##        GetVariationStats(my_dict, samplename)
        
## call function
#FilterVCFforMinCountAllele(sys.argv[1])
