import os,sys,csv,random,re,subprocess
from pathlib import Path
from subprocess import Popen, PIPE


## Usage:
## thimmamp@kw61043:~/PythonPrograms/utilities$ python GetCoreSharedGenesFromPangenomes.py ~/MAGIC16/annotation/MAGIC16_genomes_orthologs.txt
## Usage: python GetCoreSharedGenesFromPangenomes.py ~/MAGIC16/annotation/MAGIC16_genomes_orthologs.txt ~/MAGIC16/annotation/cloned.gene.list
## above includes clone id from IRGSP and numberof variants for each gene from MAGIC16 genomes

## python GetCoreSharedGenesFromPangenomes.py /home/thimmamp/MAGIC16/tmp/vcf/genome1_processed/genome1.nostar.bimultiallele.50bpindel_annotated_col1_9.vcf ~/MAGIC16/tmp/annotation/genome1_RAP_chr_modified.gff 

## testing 
## Usage: python GetVennFeatureDetails.py ~/MAGIC16/indels_results/ClonedGenes_startlost_stopgainlost_venn_output.txt ~/MAGIC16/indels_results/genome1_stop_gained_details_MAPmanandRAPannotation_withrap_cloneid.txt  ~/MAGIC16/indels_results/genome1_stop_gained_withclone_datasetsvenn_new.txt

## venn output file format
#Names	total	elements
#start_lost stop_gain stop_lost	4	OsFbox665 OsSRP-PSG rNBS41 LRR-RLK
#start_lost stop_gain	30	OsSRP-FLC OsRLCK227 YK19 OsWD40-186 OsDjC77 Xa39 R2R3-MYB OsMIR OsDof13 OsRNS8 OsSWEET15 Os_F0764 OsFbox650 OsPYL7 OsPUB65 OsWAK117 OsWAK101 Pik OsFbox572 OsMED15_1 OsUCL33 OsFbox009 OsAT18 Os_F0718 OsFbox413 OsAT17 OsSLPP1 ETR4 OsRLCK44 OsRLCK338
#start_lost stop_lost	2	HSfA6 DERF3
#stop_gain stop_lost	29	OsFbox604 OsMADS73 OsRLCK131 OsFbox211 OsMB39 OSK21 YR2 OsHSP70 OsCYP734A2 OsSub46 OsTPS13 OsATL16 Nbs1-Pi9 OsWAK10 OsPILS7b OsSTRL18 OsMADS93 Osh359-2 OsRLCK231 OsFbox081 PPR458 OsEnS-106 OsFbox667 OsSTA16 Os_F0739 LP2 OsFbox174 OsFbox028 Bph6
#start_lost	128	OsGA2ox5 Os_F0716 ONAC120 OsREM6.1 OsCPA1 OsERF#080 WIN1 OMTN3 OsKMD4 OsBBX4 OsRH38 OsIPK2 OsWAK111 OsARF4 OsWAK36 OsKMD1 OsHXK8 OsCBSSIS1 eIF-3l PsaE OsRNS6 EIL2 OsCCR21 OsERF#042 FKF1 OsWD40-199 OSRR9 OsCDC48E OsINO2 OsTET14 OsTBL27 Os4bglu10
def GetVennOut2Dict(vennout):
    outDict = {}
    total=0
    with open(vennout, 'r') as f:
        for line in f:
            if not line.startswith("Names"):
                tokens=line.rstrip('\n').split('\t')
                datasets = tokens[0].replace(' ',',')
                total+=int(tokens[1])
                features = tokens[-1].split(' ')
                for feat in features:
                    outDict[feat] = datasets
    #print("total number of features in the Venn :", str(total))
    return outDict



## combine venn out for snpeff event file
## snpeff start_lost file format
#thimmamp@kw61043:~/MAGIC16/indels_results$ head genome1_start_lost_details_MAPmanandRAPannotation_withrap_cloneid.txt
#featureid	gene	transcript	feature_coord	variant	Chromosome	variant_position	Ref	Alt	vcf_details	INDEL_classification	Function	Cellular_Component	RAPGeneSymbol	GeneDescription	RAPid	Cloneid
#Osativa.01G000170_01.exon.1	Osativa.01G000170	Osativa.01G000170_01	Chr01:88986-89204	Chr01_89203_ATCT_ATCTTCT	Chr01	89203	ATCT	ATCTTCTAC=945,423;AF=0.156,0.07;AN=6044;BaseQRankSum=0.238;DP=55783;ExcessHet=-0;FS=0;InbreedingCoeff=0.9252;MLEAC=928,406;MLEAF=0.154,0.067;MQ=60;MQRankSum=0;QD=25.59;ReadPosRankSum=0;SOR=0.628	deletion,insertion	 no annotation	NA	NA	Similar to OSIGBa0075F02.3 protein	Os01g0101900	NA
#Osativa.01G000250_01.exon.1	Osativa.01G000250	Osativa.01G000250_01	Chr01:169384-169656	Chr01_169423_AT_A	Chr01	169423	AT	A	AC=8;AF=0.001324;AN=6042;BaseQRankSum=0.537;DP=42793;ExcessHet=0;FS=3.241;InbreedingCoeff=0.1576;MLEAC=7;MLEAF=0.001159;MQ=60;MQRankSum=0;QD=30.31;ReadPosRankSum=0.397;SOR=0.415	deletion	 Light-regulated protein	 chloroplastic	Lir1|OsLIR1	Light-regulated protein%2C Regulation of light-dependent attachment of LEAF-TYPE FERREDOXIN-NADP+ OXIDOREDUCTASE (LFNR) to the thylakoid membran	Os01g0102900	LIR1
def GetCloneFunctionwithVenn(clonefunc, venndict):
    outDict = {}
    with open(clonefunc, 'r') as f:
        for line in f:
            if not line.startswith("featureid"):
                tokens = line.rstrip('\n').split('\t')
                if tokens[-1]!="NA":
                    key = tokens[-1]
                    #print("cloneid read : ", key)
                    variant = tokens[4]
                    classification = tokens[10]
                    func = tokens[11]
                    gdesc = tokens[14]
                    rapid = tokens[-2]
                    if key in outDict.keys():
                       if 'variants' in outDict[key].keys() and variant not in outDict[key]['variants']:
                           outDict[key]['variants'].append(variant)
                       if 'function' in outDict[key].keys() and func not in outDict[key]['function']:
                           outDict[key]['function'].append(func)
                       if 'genedesc' in outDict[key].keys() and gdesc not in outDict[key]['genedesc']:
                           outDict[key]['genedesc'].append(gdesc)
                       if 'rapid' in outDict[key].keys() and rapid not in outDict[key]['rapid']:
                           outDict[key]['rapid'].append(rapid)
                       outDict[key]['classification'].append(classification)
                    else:
                       outDict[key]={}
                       outDict[key]['variants']=[variant]
                       outDict[key]['function']=[func]
                       outDict[key]['genedesc']=[gdesc]
                       outDict[key]['rapid']=[rapid]
                       outDict[key]['classification']=[classification]
                    if key in venndict.keys():
                       outDict[key]['datasets']=venndict[key]

    #PrintNestedDict(outDict)
    return outDict


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

def PrintNestedDict(indict):
    for k, v in indict.items():
        #print(k,v)
        outstr=k+'\t'
        temp=''
        for k1, v1 in indict[k].items():
            #print(k, k1, v1,'\n')
            print(v1)
            if len(v1) > 1:
               temp+=';'.join(v1)+'\t'
            else:
               temp+=v1[0]+'\t'
        outstr+=temp    
        print(outstr)    

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
#vcfdict = GetVcf2Dict(sys.argv[1])
#gffdict = GetGff2Dict(sys.argv[2], "gene")
#gffwithvcf = GetVariantsForFeature(gffdict, vcfdict, "gene")
#varpath = "/home/thimmamp/MAGIC16/tmp/vcf/"
##all_genomes=['genome27','genome28', 'genome7','genome8','genome9','genome10','genome11','genome12', 'genome13', 'genome14', 'genome15','genome16','genome5','genome6','genome1','genome4']
#all_genomes=['genome1']
#rapgene2clone = GetIRGSPGenewithCloneid(sys.argv[3])
#for g in all_genomes:
#    #coresharedict = GetOrthologsofQueryGenome(g, sys.argv[1])
#    coresharedict = GetOrthologsofQueryGenome(g, sys.argv[4])
##    filepath=varpath+g+'_processed/'+g+'_numvariants_pertxpt.txt'
##    querynumvars = GetNumberofVariantsperGene(filepath)
#    #fname='/home/thimmamp/MAGIC16/indels_results/'+g+'VsMAGIC16_Orthologs_withtotalexoniccounts.txt'
#    fname='/home/thimmamp/MAGIC16/indels_results/'+g+'VsMAGIC16_Orthologs_withtotalgenicindelcounts_testing.txt'
#    #counts = CreateReportofCoreSharedGenes(coresharedict, g, rapgene2clone, querynumvars, fname)
#    counts = CreateReportofCoreSharedGenes(coresharedict, g, rapgene2clone, gffwithvcf, fname)
#    print(g, counts)

mydict = GetVennOut2Dict(sys.argv[1])
#PrintSimpleDict(mydict)
integdict = GetCloneFunctionwithVenn(sys.argv[2], mydict)
header = "CloneId\tVariants\tFunction\tGeneDescription\tRAPid\tINDELS\tINDELType\tEvents\n"
with open(sys.argv[3], 'w') as f:
    f.write(header)
    for k,v in integdict.items():
        outstr=k+'\t'
        temp=';'.join(integdict[k]['variants'])+'\t'
        func= integdict[k]['function'][0].strip()
        gdesc = integdict[k]['genedesc'][0].strip()
        temp+=func+'\t' 
        temp+=gdesc+'\t'
        temp+=integdict[k]['rapid'][0]+'\t' 
        classlist = integdict[k]['classification']
        indeltype=''
        if "insertion" in classlist and "deletion" not in classlist and "insertion,deletion" not in classlist and "deletion,insertion" not in classlist:
            indeltype="INS"
        elif "insertion" not in classlist and "deletion" in classlist and "insertion,deletion" not in classlist and "deletion,insertion" not in classlist:
            indeltype="DEL"
        #elif "insertion" in classlist and "deletion" in classlist:
        #    indeltype="BOTH"
        #elif "insertion,deletion" in classlist or "deletion,insertion" in classlist:
        #    indeltype="BOTH"
        else:
            indeltype="BOTH"
        temp+=';'.join(classlist)+'\t'   
        temp+=indeltype+'\t'
        if 'datasets' in integdict[k].keys():
            temp+=integdict[k]['datasets']
        else:    
            temp+='NA'
        outstr+=temp
        #print(outstr)
        f.write(outstr+'\n')

#PrintNestedDict(integdict)

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
