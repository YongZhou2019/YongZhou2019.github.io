import os,sys,csv,random,re,subprocess
from pathlib import Path
from subprocess import Popen, PIPE
from itertools import islice

# algorithm
# to get absent genes for IRGSP
## combine all dispensable and private genes from all genomes
## exclude IRGSP genes, remaining will be "absent" genes for IRGSP
# 1. Read dispensable table, specific tables
# exclude IRGSP genes, result will be absent in in IRGSP


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

specificgenome2file = {
    "genome10" :    "genome10.127518.KYG.gset.specific",
    "genome11" :    "genome11.127564.Lima.gset.specific",
    "genome12" :    "genome12.127652.NaBo.gset.specific",
    #"genome13" :    "genome13.127742.PR106.gset.specific",
    "genome13" :    "genome13.127742.Pr106.gset.specific",
    "genome14" :    "genome14.128077.KeNa.gset.specific",
    "genome15" :    "genome15.132278.CMeo.gset.specific",
    "genome16" :    "genome16.132424.GoSa.gset.specific",
    "genome1" :     "genome1.Nipponbare.gset.specific",
    #"genome27" :    "genome27.MH63RS3.gset.specific",
    "genome27" :    "genome27.mh63.gset.specific",
    #"genome28" :    "genome28.ZS97RS3.gset.specific",
    "genome28" :    "genome28.zs97.gset.specific",
    "genome33" :    "genome33.IR64RS2.gset.specific",
    "genome4" :     "genome4.117534.N22.gset.specific",
    #"genome5" :     "genome5.Azucena.gset.specific",
    "genome5" :     "genome5.Azu.gset.specific",
    "genome7" :     "genome7.117425.ARC.gset.specific",
    "genome8" :     "genome8.125619.LaMu.gset.specific",
    "genome9" :     "genome9.125827.LiXu.gset.specific",
        }

def GetCoreGenesForM16Orthologs(infile, mode): # mode="unique" or "all", unique will take first gene in duplicated list, all will take all ',' separated genes
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


#def GetAbsentGenesForGenomeFromDispensable(targetgenome, mode):
def GetAbsentGenesForGenomeFromDispensable(targetgenome):
    h = targetgenome ## get the genome name via command line
    absentcount=0
    #print("Number of dispensable genes from each genome with duplication\n")
    genelist=set()
    dispfile ="/home/thimmamp/MAGIC16/gramene_genome33_Only_YongUpdatedFromZhiChao/2.IR64RS2.groupset.sorted.dispensable.16.homologs" 
    with open(dispfile, 'r') as f: ## Read dispensable orthologous file
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            match = re.search(r"^\.$", row[h])
            if match:
            #if "Os" not in row[h]: # if the target genome has not ortholog is the dispensable table
            #if r"\." in row[h]: # if the target genome has not ortholog is the dispensable table
                absentcount+=1
                #print(row)
                for k in row.keys():
                    #if k!=h and row[k] is not None: #ignore target genome and get its orthologue gene from other genome which is "asbent gene"
                    if k!=h and  "Os" in row[k]: #ignore target genome and get its orthologue gene from other genome which is "asbent gene"
                        #print(row[k])
                        #if "Os" in row[k]:
                        if ',' in row[k]:
                            g = row[k].split(',')[0]
                            genelist.add(g)
                            break;
                        else:    
                            genelist.add(row[k])
                            break;
        #print(genelist) 
        #temp='\n'.join(genelist)
        #print(temp)
        genomenum=header2genome[h]
        g2f=specificgenome2file[genomenum]
#        ## now get the absent gene's genome and the number of indels from gff and vcf
        genome2absentgenes={}
        for gene in genelist:
            keylist=re.split('\_|\.', gene)
            temp=keylist[0]
            genome=temp.split('_')[0].replace("Os","")
            #print(gene, genome)
            for k,v in specificgenome2file.items():
                if genome in v:
                    #print(gene, genome, k, v,'\n')
                    if k in genome2absentgenes.keys():
                        #if gene not in genome2absentgenes[k]:
                        genome2absentgenes[k].append(gene)
                    else:
                        genome2absentgenes[k]=[gene]
                #else:
                #    print(genome, 'not in specificgenome2file','\n')
        for k1,v1 in genome2absentgenes.items():
            ## send these genes for annotation
            ## step3 get gene loci and indels within them for each of the dispensable genome's genes
            ## now get the count of indels from genes of respective genome gff and vcf
            print(k1,len(v1),'\n')
            GetindelDetailsandcounts(k1, v1, "no duplication")
        print(h+"\t"+header2genome[h]+'\t'+g2f+'\t'+str(len(genelist)))
        #print(h+"\t"+header2genome[h]+'\tnumberofabsentgenes:\t'+str(absentcount))

def Convert(lst):
    count=1
    res_dct={}
    for item in lst:
        res_dct[item]=count
    return res_dct    

def GetindelDetailsandcounts(genome, genelist, mode):
    g = genome
    genes_dict=Convert(genelist)
    annotfile = "/home/thimmamp/MAGIC16/annotation/"+g+"/genes.gff"
    print("reading annotation file ", annotfile)
    featuredict = GetAnnotationtoDictfromGffForQuery(annotfile, 'gene', genes_dict)
    print("Number of genes from gff reading: ", len(featuredict.keys()))
    vcfile = "/home/thimmamp/MAGIC16/tmp/vcf/"+g+"_processed/"+g+".nostar.bimultiallele.50bpindel_annotated_col1_9.vcf"
    print("reading vcf file ", vcfile)
    vcfdict = GetVariantClassifiction(vcfile) # input annotated vcf file
    PrintDict(vcfdict)
    ## get indel counts for each gene in the querygenesdict
    gene2svcount=GetIndelCountsForGenes(featuredict, vcfdict) # genes in a dict
    #outfile = "/home/thimmamp/MAGIC16/indels_results/"+g+"_Dispensable_with_"+mode+"__orthogenes_genic_InDelsStats_Details_testing.txt"
    
    print("Number of genes from genes2svcount reading: ", len(gene2svcount.keys()))
    #PrintDict(gene2svcount)
    #WriteDict(gene2svcount, outfile)
    #print("Done writing to file ", outfile, "...")
    print("Done indel details and counts for ", genome, " with ", len(genelist), " genes\n")


mode="no duplication"
## call functions
#GetAbsentGenesForGenomeFromDispensable(sys.argv[1], sys.argv[2])
#for g in ["MH63RS3"]:
for g in ["Nipponbare"]:
#for g in ["MH63RS3", "ZS97RS3", "117425", "125619",  "125827", "127518",	"127564", "127652", "127742", "128077", "132278","132424",  "Azucena",	"IR64RS2",  "Nipponbare","117534"]:
#GetAbsentGenesForGenomeFromDispensable(sys.argv[1])
    GetAbsentGenesForGenomeFromDispensable(g)

