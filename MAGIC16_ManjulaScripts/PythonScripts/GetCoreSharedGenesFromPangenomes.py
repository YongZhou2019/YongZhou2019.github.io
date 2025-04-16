import os,sys,csv,random,re,subprocess
from pathlib import Path
from subprocess import Popen, PIPE


## Usage:
## thimmamp@kw61043:~/PythonPrograms/utilities$ python GetCoreSharedGenesFromPangenomes.py ~/MAGIC16/annotation/MAGIC16_genomes_orthologs.txt
## Usage: python GetCoreSharedGenesFromPangenomes.py ~/MAGIC16/annotation/MAGIC16_genomes_orthologs.txt ~/MAGIC16/annotation/cloned.gene.list
## above includes clone id from IRGSP and numberof variants for each gene from MAGIC16 genomes

## python GetCoreSharedGenesFromPangenomes.py /home/thimmamp/MAGIC16/tmp/vcf/genome1_processed/genome1.nostar.bimultiallele.50bpindel_annotated_col1_9.vcf ~/MAGIC16/tmp/annotation/genome1_RAP_chr_modified.gff 

## testing 
## python GetCoreSharedGenesFromPangenomes.py /home/thimmamp/MAGIC16/tmp/vcf/genome1_processed/genome1.nostar.bimultiallele.50bpindel_annotated_col1_9.vcf ~/MAGIC16/tmp/annotation/genome1_RAP_chr_modified.gff ~/MAGIC16/annotation/cloned.gene.list ~/MAGIC16/annotation/MAGIC16_genomes_orthologs.txt


## Get annotated vcf to dict
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
#1	828	Chr01_828_ACTATTGATTTTGTC_A	ACTATTGATTTTGTC	A	17227.1	.	AC=38;AF=0.006293;AN=6038;BaseQRankSum=1.75;DP=46079;ExcessHet=0;FS=0;InbreedingCoeff=0.6095;MLEAC=36;MLEAF=0.005962;MQ=59.99;MQRankSum=0;QD=29.11;ReadPosRankSum=1.28;SOR=0.15;TYPE=del	GT
#1	959	Chr01_959_G_GA	G	GA	243898	.	AC=1367;AF=0.227;AN=6032;BaseQRankSum=-0.108;DP=39869;ExcessHet=0;FS=0;InbreedingCoeff=0.8967;MLEAC=1349;MLEAF=0.224;MQ=60;MQRankSum=0;QD=27.67;ReadPosRankSum=0.005;SOR=0.63;TYPE=ins	GT
#1	1581	Chr01_1581_GAAGCGGTGGAAGGATA_G	GAAGCGGTGGAAGGATA	G	1100.92	.	AC=8;AF=0.001324;AN=6044;DP=44983;ExcessHet=0;FS=0;InbreedingCoeff=0.2462;MLEAC=5;MLEAF=0.0008273;MQ=60;QD=29.52;SOR=0.826;TYPE=del	GT

#def GetVcf2Dict(invcf, chromosome):
def GetVcf2Dict(invcf):
    outDict={} 
    # outDict[Chr01_1581_GAAGCGGTGGAAGGATA_G]='bi_allele';'del'';'deletion'
    with open(invcf, 'r') as f:
        for line in f:
            if not line.startswith("#CHROM"):
                tokens=line.rstrip('\n').split('\t')
                alleletype=''
                indeltype=''
                #if chromosome+"_" in tokens[2]: # get chromosome wise variant
                if ',' in tokens[4]:
                    alleletype="multi_allele"
                else:
                    alleletype="bi_allele"
                detail=tokens[-2].split(';')[-1].split('=')[1]
                if 'ins' in detail and not 'del' in detail:
                    indeltype='insertion'
                elif not 'ins' in detail and 'del' in detail:
                    indeltype='deletion'
                else:
                    indeltype='both'
                value=alleletype+';'+detail+';'+indeltype
                #print(tokens[2], value)
                outDict[tokens[2]]=value
    return outDict           

#def GetGff2Dict(gffile, chromosome, feature):
def GetGff2Dict(gffile, feature):
    outDict={}
    with open(gffile, 'r') as f:
        for line in f:
            #if not line.startswith("#") and line.startswith(chromosome):
            if not line.startswith("#"):
                tokens=line.rstrip('\n').split('\t')
                if tokens[2]==feature:
                    temp=tokens[8].split(';')
                    key=''
                    if feature=="gene":
                        value=tokens[0]+":"+tokens[3]+"-"+tokens[4] #get loci of feature
                        for item in temp:
                            if "ID" in item:
                                key=item.split(":")[-1]
                        outDict[key]=value
                    elif feature=="exon":
                        exonid=''
                        for item in temp:
                            if "Parent" in item:
                                txpt=item.split(":")[-1]
                                key=txpt[0:len(txpt)-2]
                            if "Name=" in item:
                                exonid=item.split("=")[-1]
                        value=exonid+";"+tokens[0]+":"+tokens[3]+"-"+tokens[4]
                        if key in outDict.keys():
                             outDict[key].append(value)
                        else:
                             outDict[key]=[value]
    return outDict                       

## this function goes through a feature annotation with genomic coords and gets a list of all variants
## overlapping those coords
def GetVariantsForFeature(gffdict, vcfdict, feature):
    outDict = {}
    for gk, gv in gffdict.items():
        if feature=="gene":
            #print(gk,gv)
            chromosome=gv.split(":")[0]
            start = int(gv.split(":")[1].split('-')[0])
            end = int(gv.split(":")[1].split('-')[1])
            inscount=0
            delscount=0
            bothcount=0
            vardetails=[]
            for vk, vv in vcfdict.items():
                vcfchr=vk.split('_')[0]
                position = int(vk.split('_')[1])
                if vcfchr==chromosome and position >= start and position <= end:
                    indeltype=vv.split(";")[-1]
                    if indeltype=="insertion":
                        inscount+=1
                    elif indeltype=="deletion":
                        delscount+=1
                    elif indeltype=="both":
                        bothcount+=1
                    vardetails.append(vk+';'+vv)
            #print(gk, str(inscount)+';'+str(delscount)+';'+str(bothcount), vardetails)
            outDict[gk]=str(inscount)+';'+str(delscount)+';'+str(bothcount)
            #if gk not in outDict.keys():
            #    outDict[gk]={}
                #print(gk, str(inscount)+';'+str(delscount)+';'+str(bothcount), vardetails)
                #outDict[gk]["counts"]=str(inscount)+';'+str(delscount)+';'+str(bothcount)
                #outDict[gk]["detail"]=vardetails
    return outDict        

#1	RAP2018-11-26	gene	27143	28644	.	+	.	ID=gene:Os01g0100700;Name=RPS5;biotype=protein_coding;description=Similar to 40S ribosomal protein S5-1;gene_id=Os01g0100700;logic_name=rapdb_genes
#1	RAP2018-11-26	gene	35623	41136	.	+	.	ID=gene:Os01g0100900;Name=SPL1;biotype=protein_coding;description=Sphingosine-1-phosphate lyase%2C Disease resistance respons;gene_id=Os01g0100900;logic_name=rapdb_genes
#1	RAP2018-11-26	gene	82426	84095	.	+	.	ID=gene:Os01g0101700;Name=OsDjC1;biotype=protein_coding;description=Similar to chaperone protein dnaJ 20;gene_id=Os01g0101700;logic_name=rapdb_genes

def GetIRGSPGenewithCloneid(infile):
    outDict={}
    with open(infile, 'r') as f:
        for line in f:
            tokens=line.rstrip('\n').split('\t')[-1].split(';')
            rapgene=''
            clone=''
            for det in tokens:
                if "ID=" in det:
                    rapgene=det.split(':')[1]
                if "Name=" in det:
                    clone=det.split('=')[1]
            outDict[rapgene]=clone
    return outDict            
            
def GetNumberofVariantsperGene(infile):
    ## key is part of gene name from gff
    ## value is part of equivalent gene name from ortho table
    gff2Orthodict={"Os117425_":"OsARC_",
                    "Os125619_" : "OsLaMu",
                    "Os125827_" : "OsLiXu_",
                    "Os127518_" : "OsKYG_",
                    "Os127564_" : "OsLima_",
                    "Os127652_" : "OsNaBo_",
                    "Os127742_" : "OsPr106_",
                    "Os128077_" : "OsKeNa_",
                    "Os132278_" : "OsCMeo_",
                    "Os132424_" : "OsGoSa_",
            }
    outDict = {}
    with open(infile, 'r') as f:
        for line in f:
            #print(line)
            tokens=line.rstrip('\n').split()
            gene = tokens[-1][0:len(tokens[-1])-3]
            for key in gff2Orthodict.keys():
                if key in gene:
                    gene=gene.replace(key, gff2Orthodict[key])
            #print(gene, tokens[0])
            #t = re.split('\.01|\_01', tokens[-1])# 01 removed to match with annotation
            outDict[gene]=tokens[0]
    return outDict         

def PrintSimpleDict(indict):
    for k, v in indict.items():
        print(k, v, '\n')

def PrintDict(indict, g):
    corecount=0
    sharecount=0
    for k, v in indict.items():
        if len(v)==0:
            #print(k,"core")
            corecount+=1

        else:
            #print(k,v,"shared")
            sharecount+=1
    print(g, str(corecount), str(sharecount))         
    #for txptkey, components in indict[k].items():
    #    print(k, txptkey, components,'\n')
    #print(k, v, '\n')


def CreateReportofCoreSharedGenes(indict, g, clonedict, querynumvariants, outfile):
    all_genomes=['genome27','genome28',	'genome7','genome8','genome9','genome10','genome11','genome12',	'genome13', 'genome14',	'genome15','genome16','genome5','genome6','genome1','genome4']
    orthogenome=[]
    corecount=0
    dispensiblecount=0
    privatecount=0
    for genome in all_genomes:
        if genome!=g:
            orthogenome.append(genome)
    with open(outfile, 'w') as outf:
        #outf.write(g+'\t'+'\t'.join(orthogenome)+'\tType\tNumberofExonicVariants\tIRGSP_Cloneid\n')
        outf.write(g+'\t'+'\t'.join(orthogenome)+'\tType\tnumberofins\tnumberofdels\tnumberofindels\tIRGSP_Cloneid\n')
        for k,v in indict.items():
            #print(k,'\t', v)
            outstr=k+'\t'
            cloneid=''
            numofexonicvars=0
            #print(k)
            if k in querynumvariants.keys():
                #numofexonicvars=querynumvariants[k]
                numofvars=querynumvariants[k]
                numofvars=numofvars.replace(";", "\t")
            #    print(k, numofexonicvars)
            else:
                print(k, 'not in gff2vcf')
                #numofexonicvars=0
                numofvars=0
            if len(v) > 0:
                valuedict=v

                orthocategory='' # if gene is orthologous to all 16 it is core, 2-15 it is dispensible, 0 it is private
                if len(valuedict.keys())==15:
                    orthocategory="core"
                    corecount+=1
                elif len(valuedict.keys()) >= 1 and len(valuedict.keys()) <=14:
                    orthocategory="dispensible"
                    dispensiblecount+=1

                for orthog in orthogenome:
                    if orthog in valuedict.keys():
                        outstr+=valuedict[orthog]+'\t'
                    else:    
                        outstr+='NA\t'
                if "genome1" in valuedict.keys() and valuedict["genome1"] in clonedict.keys():
                    #print(valuedict['genome1'])
                    cloneid=clonedict[valuedict["genome1"]]
                else:
                    cloneid='NA'
                #print(numofexonicvars)
                #outstr+=orthocategory+'\t'+str(numofexonicvars)+'\t'+cloneid          
                outstr+=orthocategory+'\t'+str(numofvars)+'\t'+cloneid          
            else:
                #print(numofexonicvars)
                #orthocategory="private" 
                privatecount+=1
                #outstr+='\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprivate\t'+str(numofexonicvars)+'\tNA'
                outstr+='\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprivate\t'+str(numofvars)+'\tNA'
            outf.write(outstr+'\n')
            #print(outstr)   
    privatecount-=1        
    return([corecount, dispensiblecount,privatecount])

def GetOrthologsofQueryGenome(querygenome, orthologfile):
    all_genomes=['genome27','genome28',	'genome7','genome8','genome9','genome10','genome11','genome12',	'genome13', 'genome14',	'genome15','genome16','genome5','genome6','genome1','genome4']
    orthogenome=[]
    for g in all_genomes:
        if g!= querygenome:
            orthogenome.append(g)
    #print(querygenome, all_genomes, orthogenome)        
    ortho_dict={}
    with open(orthologfile, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            valuedict={} #create a dict for each orthologous genome as key and orthogene as value
            key='' # query genome's gene
            if row[querygenome]!=".":
                key = row[querygenome]
                for og in orthogenome:
                    if row[og]!=".":
                        valuedict[og]=row[og]
            ortho_dict[key]=valuedict
    return ortho_dict


## call function
vcfdict = GetVcf2Dict(sys.argv[1])
gffdict = GetGff2Dict(sys.argv[2], "gene")
gffwithvcf = GetVariantsForFeature(gffdict, vcfdict, "gene")
varpath = "/home/thimmamp/MAGIC16/tmp/vcf/"
#all_genomes=['genome27','genome28', 'genome7','genome8','genome9','genome10','genome11','genome12', 'genome13', 'genome14', 'genome15','genome16','genome5','genome6','genome1','genome4']
all_genomes=['genome1']
rapgene2clone = GetIRGSPGenewithCloneid(sys.argv[3])
for g in all_genomes:
    #coresharedict = GetOrthologsofQueryGenome(g, sys.argv[1])
    coresharedict = GetOrthologsofQueryGenome(g, sys.argv[4])
#    filepath=varpath+g+'_processed/'+g+'_numvariants_pertxpt.txt'
#    querynumvars = GetNumberofVariantsperGene(filepath)
    #fname='/home/thimmamp/MAGIC16/indels_results/'+g+'VsMAGIC16_Orthologs_withtotalexoniccounts.txt'
    fname='/home/thimmamp/MAGIC16/indels_results/'+g+'VsMAGIC16_Orthologs_withtotalgenicindelcounts_testing.txt'
    #counts = CreateReportofCoreSharedGenes(coresharedict, g, rapgene2clone, querynumvars, fname)
    counts = CreateReportofCoreSharedGenes(coresharedict, g, rapgene2clone, gffwithvcf, fname)
    print(g, counts)

## test calling function
#vcfdict = GetVcf2Dict(sys.argv[1], "Chr01")
#gffdict = GetGff2Dict(sys.argv[2], "Chr01","gene")
#vcfdict = GetVcf2Dict(sys.argv[1])
#gffdict = GetGff2Dict(sys.argv[2], "gene")
#gffwithvcf = GetVariantsForFeature(gffdict, vcfdict, "gene")
#for k,v in gffwithvcf.items():
#    for k1,v1 in gffwithvcf[k1].items():
#        print(k,v,k1,v1)
#PrintSimpleDict(gffwithvcf)
