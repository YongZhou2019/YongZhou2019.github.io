import os,sys,csv,random,re
from pathlib import Path


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
query = Queryids2Dict(sys.argv[1])
#PrintMappingQueryinGff(query, sys.argv[2])
PrintMappingQueryRegioninGff(query, sys.argv[2])
#genetxptcomponents = GetGffAnnotationToDict(sys.argv[1])
#PrintDict(genetxptcomponents)



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
