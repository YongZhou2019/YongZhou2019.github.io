import os
from pathlib import Path


folder = Path('/home/thimmamp/MAGIC16/indels_results/')

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



def GetGenicregionBedForQueryId(indict, annotfile, genome, outfile):
    # Chr12 cshl_gene   gene    24526624    24527466    .   -   .   ID=gene:OsIR64_12g0020410;biotype=protein_coding;logic_name=cshl_gene
    with open(annotfile, 'r') as f:
        with open(outfile, 'w') as of:
            for line in f:
                if not line.startswith('#'):
                    tokens=line.rstrip('\n').split('\t')
                    #if "gene" in tokens[2]:# for gene
                    if "mRNA" in tokens[2]:# for gene
                        txpt=''
                        gene=''
                        for temp in tokens[-1].split(';'):
                            if "ID=" in temp:
                                txpt = temp.split(':')[1]
                            if "Parent=" in temp:
                                gene = temp.split(':')[1]
                        for k in indict.keys():
                            if k==gene:
                            #if k in tokens[-1]:
                                #print(line)
                                #annotid=''
                                annotid=txpt
                                #if 'gene' in tokens[-1].split(';')[0]:
                                #    annotid = tokens[-1].split(';')[0].split(':')[1]
                                #else:
                                #    annotid = tokens[-1].split(';')[0].split('=')[1]
                                    
                                #if k==annotid:
                                annotid = genome+"_"+annotid
                                of.write(tokens[0]+'\t'+tokens[3]+'\t'+tokens[4]+'\t'+annotid+'\t0\t'+tokens[6]+'\n')



outdict={}
for f in folder.glob("Nipponbare_AbsentGenes_from_genome*_with_mode_duplication*.txt"):
    filename = os.path.basename(f)
    genome = filename.split('_')[3]
    #print(f, '\n', filename,'\n', genome)
    with open(f , "r") as f1:
        for line in f1:
            if not line.startswith('Geneid'):
                gene = line.rstrip('\n').split('\t')[0]
                ## to map gff gene name with ortholog short name
                for k ,v in gff2shortname.items():
                    if v in gene:
                        gene=gene.replace(v , k)
                        
                if genome not in outdict.keys():
                    outdict[genome] = [gene]
                else:
                    outdict[genome].append(gene)

                    

for k,v in outdict.items():
    print("Processing Absent in Nipponebare genome,but present in genome ", k) 
    print("Genome ", k, " has ", len(v), "genes")
    gene_dict = {value: len(value) for value in v}
    outputfile = "/home/thimmamp/MAGIC16/indels_results/Nipponbare_absent_genesin_transcripts_"+k+".bed"
    GetGenicregionBedForQueryId(gene_dict, "/home/thimmamp/MAGIC16/annotation/"+k+"/genes.gff", k, outputfile)
    print("Processing Absent in Nipponebare genome,but present in genome ", k, " done\n") 


