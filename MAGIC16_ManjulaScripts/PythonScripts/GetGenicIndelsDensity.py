import os,sys,csv,random,re,subprocess
from pathlib import Path
from subprocess import Popen, PIPE
from itertools import islice

## Usage
## thimmamp@KW61043:~/PythonPrograms/utilities$ python GetGenicIndelsDensity.py /home/thimmamp/MAGIC16/indels_results/genome1_specefic_genic_InDelsStats_Details_testing.txt > MergedAgriTraitGenes_ingenome1Specificgenes_zerofiltered.log
## python GetGenicIndelsDensity.py ~/MAGIC16/annotation/Oryza_sativa.IRGSP-1.0.48.gff3 ~/MAGIC16/indels_results/genome1_Dispensable_with_all__orthogenes_genic_InDelsStats_Details_testing.txt testingout.txt
## thimmamp@KW61043:~/PythonPrograms/utilities$ python GetGenicIndelsDensity.py ~/MAGIC16/annotation/Oryza_sativa.IRGSP-1.0.48.gff3 ~/MAGIC16/indels_results/genome1_specefic_genic_InDelsStats_Details_testing.txt ~/MAGIC16/indels_results/genome1_Specific_withRAP_MapmanAnnot_vareffects.txt

## thimmamp@KW61043:~/PythonPrograms/utilities$ python GetGenicIndelsDensity.py ~/MAGIC16/annotation/Oryza_sativa.IRGSP-1.0.48.gff3 ~/MAGIC16/annotation/geneKeyword.table.txt ~/MAGIC16/indels_results/genome1_specefic_genic_InDelsStats_Details_testing.txt ~/MAGIC16/indels_results/genome1_Specific_withRAP_MapmanAnnot_vareffects.txt


## thimmamp@KW61043:~/PythonPrograms/utilities$ python GetGenicIndelsDensity.py ~/MiAGIC16/annotation/geneKeyword.table.txt ~/MAGIC16/indels_results/genome1_Specific_withRAP_MapmanAnnot.txt
## thimmamp@KW61043:~/PythonPrograms/utilities$ python GetGenicIndelsDensity.py ~/MAGIC16/annotation/geneKeyword.table.txt ~/MAGIC16/indels_results/genome1_Specific_withRAP_MapmanAnnot.txt > temp1.txt

## thimmamp@KW61043:~/PythonPrograms/utilities$ python GetGenicIndelsDensity.py ~/MAGIC16/annotation/geneKeyword.table.txt ~/MAGIC16/indels_results/genome1_specific_temp  > temp2.txt

## thimmamp@KW61043:~/PythonPrograms/utilities$ python GetGenicIndelsDensity.py ~/MAGIC16/annotation/Oryza_sativa.IRGSP-1.0.48.gff3 ~/MAGIC16/annotation/geneKeyword.table.txt ~/MAGIC16/indels_results/genome1_specefic_genic_InDelsStats_Details_testing.txt ~/MAGIC16/indels_results/g1_Specific_withRAP_MapmanAnnot_vareffects_Keywords_Evidences_testing.txt

mergefile="/home/thimmamp/MAGIC16/MergedListMappedtoIRGSP.txt"
rapnipfile="/home/thimmamp/MAGIC16/annotation/Nipponbare_IRGSPRAP.one2one"
rapgramenefile="/home/thimmamp/MAGIC16/annotation/IRGSPgramene_IRGSPRAP.one2one"
mapmanfile="/home/thimmamp/MAGIC16/indels_results/mapman_annotation_forparsing.csv"
varianteffectfile = "/home/thimmamp/MAGIC16/tmp/SnpEffResults/genome1_snpEff/genome1.nostar.bimultiallele.50bpindel_annotated_col1_9.csv2.genes.txt"

## contents of IRGSPgramene_RAP
#Osativa.01G000010       Os01g0100100    RBH     ChrN
#Osativa.01G000020       Os01g0100200    RBH     ChrN
#Osativa.01G000030       Os01g0100300    RBH     ChrN

## contentst of /home/thimmamp/MAGIC16/annotation/Nipponbare_IRGSPRAP.one2one
#OsNip_01G000010	Os01g0100100	RBH	ChrN
#OsNip_01G000020	Os01g0100200	RBH	ChrN
#OsNip_01G000030	Os01g0100300	RBH	ChrN
#OsNip_01G000040	Os01g0100400	RBH	ChrN

# input file format
#GeneName        GeneId  TranscriptId    BioType variants_impact_HIGH    variants_impact_LOW     variants_impact_MODERATE        variants_impact_MODIFIER        variants_effect_3_prime_UTR_truncation  variants_effect_3_prime_UTR_variant     variants_effect_5_prime_UTR_truncation  variants_effect_5_prime_UTR_variant     variants_effect_conservative_inframe_deletion   variants_effect_conservative_inframe_insertion  variants_effect_disruptive_inframe_deletion     variants_effect_disruptive_inframe_insertion    variants_effect_downstream_gene_variant variants_effect_exon_loss_variant       variants_effect_frameshift_variant      variants_effect_intron_variant  variants_effect_non_coding_transcript_variant   variants_effect_splice_acceptor_variant variants_effect_splice_donor_variant    variants_effect_splice_region_variant   variants_effect_start_lost      variants_effect_start_retained_variant  variants_effect_stop_gained     variants_effect_stop_lost       variants_effect_stop_retained_variant   variants_effect_upstream_gene_variant
#Osir64_01g0000010.01.exon.1     GENE_Osir64_01g0000010.01.exon.1        transcript:Osir64_01g0000010.01 protein_coding  0       0       0       96      0       1       0       0       0       0       0       0       33      0       0       21      0       0       0       0       0       0       0       0       0       41
#Osir64_01g0000020.01.exon.1     GENE_Osir64_01g0000020.01.exon.1        transcript:Osir64_01g0000020.01 protein_coding  7       0       5       38      0       0       0       5       0       4       0       1       17      0       7       0       0       0       0       0       0       0       1       0       0       16
#Osir64_01g0000030.01.exon.1     GENE_Osir64_01g0000030.01.exon.1        transcript:Osir64_01g0000030.01 protein_coding  10      0       7       33      0       0       0       0       1       5       1       0       16      0       10      0       0       0       0       0       0       0       0       0       0       17
## function parses this files, looks for variants_impact_HIGH > 1 and non-zero value for each variants_effects and write one file per variants_effect
## Usage : thimmamp@kw61043:~/PythonPrograms/utilities$ python AnnotateMagic16Genomes_indels_dev.py ~/MAGIC16/tmp/vcf/genome6_processed/genome6.nostar.bimultiallele.50bpindel_annotated_col1_9.csv2.genes.txt

#def GetVariantEffectwithHighImpact(inputgenesfile, variant_effect_type, genome):
def GetVariantEffectwithHighImpact(inputgenesfile, variant_effect_type):
    outDict = {}
    with open(inputgenesfile, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
                #print(row['variants_impact_HIGH'], row[effect])
            if int(row['variants_impact_HIGH']) > 1 and int(row[variant_effect_type]) > 0:
                temp = row['GeneName'].split("_")[0]
                #outDict[temp]=row['variants_impact_HIGH']
                outDict[temp]=row[variant_effect_type]
    return outDict   

def ReadMergedList(infile):
    outDict = {}
    with open(infile, 'r') as f:
        for line in f:
            #print(line)
            token=line.rstrip('\n').split('\t')[-1]
            id='';name='';desc='';
            for item in token.split(';'):
                if "ID=gene" in item:
                    id = item.split(':')[1]
                if "Name=" in item:
                    name = item.split('=')[1]
                if "description=" in item:
                    desc = item.split('=')[1]
            if name:        
               outDict[id]=name+';'+desc
            else:
                outDict[id]='NA;'+desc
    return outDict

def ReadGffFile(gffile):
#    awk -F'\t' '{if($1==10&&$3=="gene")print}' ../annotation/Oryza_sativa.IRGSP-1.0.48.gff3 | head
#10	RAP2018-11-26	gene	44902	45692	.	+	.	ID=gene:Os10g0100200;biotype=protein_coding;description=Hypothetical protein;gene_id=Os10g0100200;logic_name=rapdb_genes
#10	RAP2018-11-26	gene	46031	58278	.	-	.	ID=gene:Os10g0100300;Name=TIC62;biotype=protein_coding;description=NAD(P)-binding domain containing protein;gene_id=Os10g0100300;logic_name=rapdb_genes
#10	RAP2018-11-26	gene	67119	72969	.	+	.	ID=gene:Os10g0100500;Name=OsPUB61;biotype=protein_coding;description=Serine/threonine protein kinase-related domain containing protein;gene_id=Os10g0100500;logic_name=rapdb_genes
#10	RAP2018-11-26	gene	81102	82243	.	-	.	ID=gene:Os10g0100700;Name=OsPDX1.3b;biotype=protein_coding;description=Vitamin B6 biosynthesis protein family protein;gene_id=Os10g0100700;logic_name=rapdb_genes
    outDict = {}
    with open(gffile, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                tokens = line.rstrip('\n').split('\t')
                if tokens[2]=="gene":
                    chromosome=''
                    annotsource = tokens[1]
                    start = tokens[3]
                    end = tokens[4]
                    strand = tokens[6]
                    details=tokens[-1].split(';')
                    geneid = ''
                    name=''
                    description=''
                    if tokens[0].isnumeric():
                    #if int(tokens[0])>=1 and int(tokens[0])<=12:
                        if int(tokens[0])>=1 and int(tokens[0])<=9:
                            chromosome="Chr0"+tokens[0]
                        else:
                            chromosome="Chr"+tokens[0]
                        for temp in details:
                            if "ID=" in temp:
                                geneid = temp.split(':')[1]
                            if "Name=" in temp:
                                name = temp.split("=")[1]
                            if "description=" in temp:
                                description = temp.split("=")[1]
                        if name :        
                            outDict[geneid]=chromosome+';'+annotsource+';'+start+';'+end+';'+strand+';'+name+';'+description
                        else:    
                            outDict[geneid]=chromosome+';'+annotsource+';'+start+';'+end+';'+strand+';NA;'+description
    return outDict

def ReadFile2Dict(infile, refindex, queryindex):
    outDict={}
    with open(infile, 'r') as f:
        for line in f:
            #print(line)
            tokens=line.rstrip('\n').split('\t')
            #outDict[tokens[1]]=tokens[0]
            outDict[tokens[int(refindex)]]=tokens[int(queryindex)]
    return outDict

def PrintDict(indict):
    for k, v in indict.items():
        print(k, v)

def ReadInDelDetails(infile, rap2nipdict, rap2gramenedict):
    outDict = {}
    with open(infile, 'r') as f:
        for line in f:
            if not line.startswith("Geneid"):
                tokens = line.rstrip('\n').split('\t')
                temp='NA'

                for k,v in rap2nipdict.items():
                    if v==tokens[0]:
                        temp=k
                        break;
                #for k,v in rap2gramenedict.items(): #if rap2gramene
                    #if v==tokens[0]:
                        #temp=k
                for k,v in rap2gramenedict.items(): #if gramene2rap
                    if k==tokens[0]:
                        temp=v
                        break;
                
                outDict[tokens[0]] = temp+'::'+'::'.join(tokens[1:])
    return outDict        


def GetMapmanAnnot2Dict(annotfile):
    outDict = {}
    with open(annotfile, 'r') as f:
        for line in f:
            if not line.startswith("BINCODE"):
                tokens = line.rstrip('\n').split('\t')
                if tokens[2]!="":
                    function=''
                    #cellcomponent=''
                    tokens[2]= tokens[2].replace('o','O')
                    tokens[2]= tokens[2].replace('g','G')
                    gene = tokens[2]
                    temp = tokens[3].split(':')
                    funcs= temp[3].split(' & ')
                    if ',' in funcs[0]:
                        function = funcs[0].split(',')[0]
                        #cellcomponent = funcs[0].split(',')[1]
                    else:
                        function = funcs[0]
                        #cellcomponent = 'NA'
                    #outDict[gene]=function+';'+cellcomponent
                    outDict[gene]=function
    #PrintDict(outDict)
    return outDict


def GetKeywordEvidenceForRAPId(infile):
    outdict={}
    #with open(infile, 'r') as f:
    with open(infile, 'r', encoding='unicode_escape') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            #print("rapid:" ,row['RAPdb'], "keyword:", row['Keyword'], "evidence:", row['Evidence'])
            key=row['RAPdb']
            kw = row['Keyword']
            evi = row['Evidence']
            if key not in outdict.keys():
                outdict[key]={}
                outdict[key]['keywords']=[kw]
                outdict[key]['evidences']=[evi]
            else:
                if kw not in outdict[key]['keywords']:
                    outdict[key]['keywords'].append(kw)
                if evi not in outdict[key]['evidences']:
                    outdict[key]['evidences'].append(evi)
    return outdict

def GetValueforqueryinDict(query, indict):
    val = '0'
    if query in indict.keys():
        val = str(indict[query])
    return val

## call functions
#mergedict = ReadMergedList(mergefile)
#gramene2rapid = ReadFile2Dict(rapgramenefile, 0, 1)

#rap2nipid = ReadFile2Dict(rapnipfile, 1, 0) #RAPid to OsNip
#rap2gramene = ReadFile2Dict(rapgramenefile, 1, 0) # RAPid to gramene Osativa.xxxx
#mapmandict = GetMapmanAnnot2Dict(mapmanfile)

#gffdict = ReadGffFile(sys.argv[1]) ## gff3 file
rap2kwevi = GetKeywordEvidenceForRAPId(sys.argv[1])
#PrintDict(rap2kwevi)
#indelstatsdict = ReadInDelDetails(sys.argv[3], rap2nipid, rap2gramene) # core, dispensable, specific file , rapid as key, Nip/gramene in the val, with indels

#PrintDict(indelstatsdict)

###################
## variant effect
## working with snpeff output genes file for variant effect
#variant_effects = ['variants_effect_3_prime_UTR_truncation', 'variants_effect_3_prime_UTR_variant', 'variants_effect_5_prime_UTR_truncation', 'variants_effect_5_prime_UTR_variant', 'variants_effect_conservative_inframe_deletion', 'variants_effect_conservative_inframe_insertion', 'variants_effect_disruptive_inframe_deletion', 'variants_effect_disruptive_inframe_insertion', 'variants_effect_downstream_gene_variant', 'variants_effect_exon_loss_variant', 'variants_effect_frameshift_variant', 'variants_effect_intron_variant', 'variants_effect_non_coding_transcript_variant', 'variants_effect_splice_acceptor_variant', 'variants_effect_splice_donor_variant', 'variants_effect_splice_region_variant', 'variants_effect_start_lost', 'variants_effect_start_retained_variant', 'variants_effect_stop_gained', 'variants_effect_stop_lost', 'variants_effect_stop_retained_variant', 'variants_effect_upstream_gene_variant']
#variant_effects = ['variants_effect_exon_loss_variant',  'variants_effect_start_lost', 'variants_effect_start_retained_variant', 'variants_effect_stop_gained', 'variants_effect_stop_lost', 'variants_effect_stop_retained_variant']

exon_loss_var_dict = GetVariantEffectwithHighImpact(varianteffectfile, 'variants_effect_exon_loss_variant')
start_lost_var_dict = GetVariantEffectwithHighImpact(varianteffectfile, 'variants_effect_start_lost')
start_retained_var_dict = GetVariantEffectwithHighImpact(varianteffectfile, 'variants_effect_start_retained_variant')
stop_gained_var_dict = GetVariantEffectwithHighImpact(varianteffectfile, 'variants_effect_stop_gained')
stop_lost_var_dict = GetVariantEffectwithHighImpact(varianteffectfile, 'variants_effect_stop_lost')
threeputr_trunc_dict = GetVariantEffectwithHighImpact(varianteffectfile, 'variants_effect_3_prime_UTR_truncation')
threeputr_variant_dict = GetVariantEffectwithHighImpact(varianteffectfile, 'variants_effect_3_prime_UTR_variant')
fiveputr_trunc_dict = GetVariantEffectwithHighImpact(varianteffectfile, 'variants_effect_5_prime_UTR_truncation')
fiveputr_variant_dict = GetVariantEffectwithHighImpact(varianteffectfile, 'variants_effect_5_prime_UTR_variant')
cons_inframe_del_dict = GetVariantEffectwithHighImpact(varianteffectfile, 'variants_effect_conservative_inframe_deletion')
cons_inframe_ins_dict = GetVariantEffectwithHighImpact(varianteffectfile,'variants_effect_conservative_inframe_insertion')
disp_inframe_del_dict = GetVariantEffectwithHighImpact(varianteffectfile, 'variants_effect_disruptive_inframe_deletion')
disp_inframe_ins_dict = GetVariantEffectwithHighImpact(varianteffectfile, 'variants_effect_disruptive_inframe_insertion')
downstream_genevar_dict = GetVariantEffectwithHighImpact(varianteffectfile, 'variants_effect_downstream_gene_variant')
frameshift_var_dict = GetVariantEffectwithHighImpact(varianteffectfile, 'variants_effect_frameshift_variant')
intron_var_dict = GetVariantEffectwithHighImpact(varianteffectfile, 'variants_effect_intron_variant')
nctxpt_var_dict = GetVariantEffectwithHighImpact(varianteffectfile, 'variants_effect_non_coding_transcript_variant')
splice_acc_var_dict = GetVariantEffectwithHighImpact(varianteffectfile, 'variants_effect_splice_acceptor_variant')
splice_don_var_dict = GetVariantEffectwithHighImpact(varianteffectfile, 'variants_effect_splice_donor_variant')
splice_region_var_dict = GetVariantEffectwithHighImpact(varianteffectfile, 'variants_effect_splice_region_variant')
stop_retained_var_dict = GetVariantEffectwithHighImpact(varianteffectfile, 'variants_effect_stop_retained_variant')
upstream_genevar_dict = GetVariantEffectwithHighImpact(varianteffectfile, 'variants_effect_upstream_gene_variant')


######################
Type="Specific"

### Add additional variant effects, Keywords and Evidences to Specific group of genes
## start with file  ~/MAGIC16/indels_results/genome1_Specific_withRAP_MapmanAnnot.txt
## format of file genome1_Specific_withRAP_MapmanAnnot.txt
#NiporGrameneId	RAPId	Chromosome	start	end	Name	Description	FunctionalAnnotation	INS	DELS
#Osativa.01G003200	NA	No Chr	No start	No end	No name	No description	NA	14	6
#Osativa.01G003280	Os01g0146801	Chr01	2541950	2542484	NA	Conserved hypothetical protein	NA	2	1
#Osativa.01G005370	NA	No Chr	No start	No end	No name	No description	NA	3	6
#Osativa.01G010840	Os01g0258100	Chr01	8630485	8633592	NA	Hypothetical conserved gene	NA	3	7
#Osativa.01G012640	Os01g0830500	Chr01	35525560	35527592	Os2OG-Fe	2OG-Fe(II) oxygenase domain containing protein	NA	95
#Osativa.01G012900	NA	No Chr	No start	No end	No name	No description	NA	6	13

print("NiporGrameneId\tRAPId\tINDELDetails\tChromosome\tstart\tend\tName\tDescription\tFunctionalAnnotation\tINS\tDELS\texon_loss\tstart_lost\tstop_gained\tstop_lost\tstart_retained\tthreeputr_trunc\tthreeputr_variant\tfiveputr_trunc\tfiveputr_variant\tcons_inframe_del\tcons_inframe_ins\tdisp_inframe_del\tdisp_inframe_ins\tdownstream_genevar\tframeshift\tintron\tnctxpt\tsplice_acc\tsplice_don\tsplice_region\tstop_retained\tupstream_genevar\tType\tKeywords\tEvidences\n")
with open(sys.argv[2], 'r') as f:
    for line in f:
        if not line.startswith('NiporGrameneId'):
            tokens = line.rstrip('\n').split('\t')
            r2g = tokens[0]
            kw = "No Keywords"
            evi = "No Evidences"
            if 'Os' in tokens[1]:
                if tokens[1] in rap2kwevi.keys():
                    kw = ','.join(rap2kwevi[tokens[1]]['keywords'])
                    evi = ','.join(rap2kwevi[tokens[1]]['evidences'])
                
            exon_loss = GetValueforqueryinDict(r2g, exon_loss_var_dict)    
            start_lost = GetValueforqueryinDict(r2g, start_lost_var_dict)
            stop_gain = GetValueforqueryinDict(r2g, stop_gained_var_dict)
            stop_lost = GetValueforqueryinDict(r2g, stop_lost_var_dict)
            start_retained = GetValueforqueryinDict(r2g, start_retained_var_dict) 
            threeputr_trunc = GetValueforqueryinDict(r2g, threeputr_trunc_dict) 
            threeputr_variant = GetValueforqueryinDict(r2g, threeputr_variant_dict) 
            fiveputr_trunc = GetValueforqueryinDict(r2g, fiveputr_trunc_dict) 
            fiveputr_variant = GetValueforqueryinDict(r2g, fiveputr_variant_dict) 
            cons_inframe_del = GetValueforqueryinDict(r2g, cons_inframe_del_dict) 
            cons_inframe_ins = GetValueforqueryinDict(r2g, cons_inframe_ins_dict) 
            disp_inframe_del = GetValueforqueryinDict(r2g, disp_inframe_del_dict) 
            disp_inframe_ins = GetValueforqueryinDict(r2g, disp_inframe_ins_dict)
            downstream_genevar = GetValueforqueryinDict(r2g, downstream_genevar_dict) 
            frameshift = GetValueforqueryinDict(r2g, frameshift_var_dict) 
            intron = GetValueforqueryinDict(r2g, intron_var_dict) 
            nctxpt = GetValueforqueryinDict(r2g, nctxpt_var_dict) 
            splice_acc = GetValueforqueryinDict(r2g, splice_acc_var_dict) 
            splice_don = GetValueforqueryinDict(r2g, splice_don_var_dict) 
            splice_region = GetValueforqueryinDict(r2g, splice_region_var_dict) 
            stop_retained = GetValueforqueryinDict(r2g, stop_retained_var_dict)
            upstream_genevar = GetValueforqueryinDict(r2g, upstream_genevar_dict)
            print('\t'.join(tokens[0:11]) +'\t'+exon_loss+'\t'+start_lost+'\t'+stop_gain+'\t'+stop_lost+'\t'+start_retained+'\t'+threeputr_trunc+'\t'+threeputr_variant+'\t'+fiveputr_trunc+'\t'+fiveputr_variant+'\t'+cons_inframe_del+'\t'+cons_inframe_ins+'\t'+disp_inframe_del+'\t'+disp_inframe_ins+'\t'+downstream_genevar+'\t'+frameshift+'\t'+intron+'\t'+nctxpt+'\t'+splice_acc+'\t'+splice_don+'\t'+splice_region+'\t'+stop_retained+'\t'+upstream_genevar+'\t'+Type+'\t'+kw+'\t'+evi)

###################################
'''

## Get indels with annot details for Core, Dispensable, and Specific group of genes
#Geneid  ListofVariants  NumberofINS     NumberofDELS
#OsNip_01G000030 Chr01_11454_TCGTTTAGCATG_T;TCGTTTAGCATG;T;bi_allele;del;316;;Chr01_11501_G_GTCGA;G;GTCGA;bi_allele;ins;65;;Chr01_11524_TCACCTCCTC_T;TCACCTCCTC;T;bi_allele;del;6;;Chr01_11663_C_CCA;C;CCA;bi_allele;ins;839;;Chr01_11693_GCCGCCGCCTTCACA_G;GCCGCCGCCTTCACA;G;bi_allele;del;69;;Chr01_11780_GGCGCGGTCC_G;GGCGCGGTCC;G;bi_allele;del;49;;Chr01_12249_T_TCG;T;TCG,TCA;;ins,ins;17,10;;Chr01_12250_T_TGTGGCG;T;TGTGGCG,TGTGGCGGAATACGGCGGGCTTGCCCA,TGGCGGAA,TGGCG;;ins,ins,ins,ins;154,52,4,1;;Chr01_12251_G_GTGGCGGAATACGGCGGGCTTGCCCAGAGCTTGCCGA;G;GTGGCGGAATACGGCGGGCTTGCCCAGAGCTTGCCGA,GGCGGA,GTGGCGGAATACGGCGGGCTTGCCC,GT,GTGGCGGAAT;;ins,ins,ins,ins,ins;66,13,35,5,7;;Chr01_12253_G_GA;G;GA,GGCGGAATA,GGCGGAATACGGCGGGCTTGCC,GAGCTTGCCGAAGGCGCCGA,GAGCTTGCCGAAGGCGC,GGCGGAATACGGCGGGCTTGC;;ins,ins,ins,ins,ins,ins;4,12,50,7,28,8    2       4

with open(sys.argv[4], 'w') as of:
    #of.write("RAPId\tNiporGrameneId\tName\tDescription\tINS\tDELS\n")
    #of.write("NiporGrameneId\tRAPId\tINDELDetails\tChromosome\tstart\tend\tName\tDescription\tFunctionalAnnotation\tINS\tDELS\texon_loss\tstart_lost\tstop_gained\tstop_lost\tstart_retained\tthreeputr_trunc\tthreeputr_variant\tfiveputr_trunc\tfiveputr_variant\tcons_inframe_del\tcons_inframe_ins\tdisp_inframe_del\tdisp_inframe_ins\tdownstream_genevar\tframeshift\tintron\tnctxpt\tsplice_acc\tsplice_don\tsplice_region\tstop_retained\tupstream_genevar\tType\tKeywords\tEvidences\n")
    print("NiporGrameneId\tRAPId\tINDELDetails\tChromosome\tstart\tend\tName\tDescription\tFunctionalAnnotation\tINS\tDELS\texon_loss\tstart_lost\tstop_gained\tstop_lost\tstart_retained\tthreeputr_trunc\tthreeputr_variant\tfiveputr_trunc\tfiveputr_variant\tcons_inframe_del\tcons_inframe_ins\tdisp_inframe_del\tdisp_inframe_ins\tdownstream_genevar\tframeshift\tintron\tnctxpt\tsplice_acc\tsplice_don\tsplice_region\tstop_retained\tupstream_genevar\tType\tKeywords\tEvidences\n")
    #of.write("NiporGrameneId\tRAPId\tINDELDetails\tChromosome\tstart\tend\tName\tDescription\tFunctionalAnnotation\tINS\tDELS\texon_loss\tstart_lost\tstop_gained\tstop_lost\n")
    #for k, v in mergedict.items(): # k is RAPid, v is Name, description
    for k, v in indelstatsdict.items(): # k is Nip or Grameneid, v is rapid'::'detail'::'ins'::'dels
        #print(k, v)
        rapid, detail,ins,dels = v.split('::')
        function='NA'
        r2g = ''
        exon_loss='0'
        start_lost='0'
        stop_gain='0'
        stop_lost='0'
        if rapid in rap2gramene.keys():
            r2g=rap2gramene[rapid]
            if r2g in mapmandict.keys():
                function = mapmandict[r2g]
            exon_loss = GetValueforqueryinDict(r2g, exon_loss_var_dict)    
            start_lost = GetValueforqueryinDict(r2g, start_lost_var_dict)
            stop_gain = GetValueforqueryinDict(r2g, stop_gained_var_dict)
            stop_lost = GetValueforqueryinDict(r2g, stop_lost_var_dict)
            start_retained = GetValueforqueryinDict(r2g, start_retained_var_dict) 
            threeputr_trunc = GetValueforqueryinDict(r2g, threeputr_trunc_dict) 
            threeputr_variant = GetValueforqueryinDict(r2g, threeputr_variant_dict) 
            fiveputr_trunc = GetValueforqueryinDict(r2g, fiveputr_trunc_dict) 
            fiveputr_variant = GetValueforqueryinDict(r2g, fiveputr_variant_dict) 
            cons_inframe_del = GetValueforqueryinDict(r2g, cons_inframe_del_dict) 
            cons_inframe_ins = GetValueforqueryinDict(r2g, cons_inframe_ins_dict) 
            disp_inframe_del = GetValueforqueryinDict(r2g, disp_inframe_del_dict) 
            disp_inframe_ins = GetValueforqueryinDict(r2g, disp_inframe_ins_dict)
            downstream_genevar = GetValueforqueryinDict(r2g, downstream_genevar_dict) 
            frameshift = GetValueforqueryinDict(r2g, frameshift_var_dict) 
            intron = GetValueforqueryinDict(r2g, intron_var_dict) 
            nctxpt = GetValueforqueryinDict(r2g, nctxpt_var_dict) 
            splice_acc = GetValueforqueryinDict(r2g, splice_acc_var_dict) 
            splice_don = GetValueforqueryinDict(r2g, splice_don_var_dict) 
            splice_region = GetValueforqueryinDict(r2g, splice_region_var_dict) 
            stop_retained = GetValueforqueryinDict(r2g, stop_retained_var_dict)
            upstream_genevar = GetValueforqueryinDict(r2g, upstream_genevar_dict)

        if int(ins)!=0 and int(dels)!=0:
            print("RAPId: ", rapid)
            if rapid in gffdict.keys():
                chromosome, source, start, end, strand, name, description = gffdict[rapid].split(';')
                #of.write(k+'\t'+rapid+'\t'+detail+'\t'+chromosome+'\t'+start+'\t'+end+'\t'+name+'\t'+description+'\t'+function+'\t'+ins+'\t'+dels+'\t'+exon_loss+'\t'+start_lost+'\t'+stop_gain+'\t'+stop_lost+'\t'+start_retained+'\t'+threeputr_trunc+'\t'+threeputr_variant+'\t'+fiveputr_trunc+'\t'+fiveputr_variant+'\t'+cons_inframe_del+'\t'+cons_inframe_ins+'\t'+disp_inframe_del+'\t'+disp_inframe_ins+'\t'+downstream_genevar+'\t'+frameshift+'\t'+intron+'\t'+nctxpt+'\t'+splice_acc+'\t'+splice_don+'\t'+splice_region+'\t'+stop_retained+'\t'+upstream_genevar+'\t'+Type+'\t'+kw+'\t'+evi+'\n')
                print(k+'\t'+rapid+'\t'+detail+'\t'+chromosome+'\t'+start+'\t'+end+'\t'+name+'\t'+description+'\t'+function+'\t'+ins+'\t'+dels+'\t'+exon_loss+'\t'+start_lost+'\t'+stop_gain+'\t'+stop_lost+'\t'+start_retained+'\t'+threeputr_trunc+'\t'+threeputr_variant+'\t'+fiveputr_trunc+'\t'+fiveputr_variant+'\t'+cons_inframe_del+'\t'+cons_inframe_ins+'\t'+disp_inframe_del+'\t'+disp_inframe_ins+'\t'+downstream_genevar+'\t'+frameshift+'\t'+intron+'\t'+nctxpt+'\t'+splice_acc+'\t'+splice_don+'\t'+splice_region+'\t'+stop_retained+'\t'+upstream_genevar+'\t'+Type+'\t'+kw+'\t'+evi+'\n')
                #of.write(k+'\t'+rapid+'\t'+detail+'\t'+chromosome+'\t'+start+'\t'+end+'\t'+name+'\t'+description+'\t'+function+'\t'+ins+'\t'+dels+'\t'+exon_loss+'\t'+start_lost+'\t'+stop_gain+'\t'+stop_lost+'\n')
            else:    
                #of.write(k+'\t'+rapid+'\t'+detail+'\t'+'No Chr'+'\t'+'No start'+'\t'+'No end'+'\t'+'No name'+'\t'+'No description'+'\t'+function+'\t'+ins+'\t'+dels+'\t'+exon_loss+'\t'+start_lost+'\t'+stop_gain+'\t'+stop_lost+'\t'+start_retained+'\t'+threeputr_trunc+'\t'+threeputr_variant+'\t'+fiveputr_trunc+'\t'+fiveputr_variant+'\t'+cons_inframe_del+'\t'+cons_inframe_ins+'\t'+disp_inframe_del+'\t'+disp_inframe_ins+'\t'+downstream_genevar+'\t'+frameshift+'\t'+intron+'\t'+nctxpt+'\t'+splice_acc+'\t'+splice_don+'\t'+splice_region+'\t'+stop_retained+'\t'+upstream_genevar+'\t'+Type+'\t'+kw+'\t'+evi+'\n')
                print(k+'\t'+rapid+'\t'+detail+'\t'+'No Chr'+'\t'+'No start'+'\t'+'No end'+'\t'+'No name'+'\t'+'No description'+'\t'+function+'\t'+ins+'\t'+dels+'\t'+exon_loss+'\t'+start_lost+'\t'+stop_gain+'\t'+stop_lost+'\t'+start_retained+'\t'+threeputr_trunc+'\t'+threeputr_variant+'\t'+fiveputr_trunc+'\t'+fiveputr_variant+'\t'+cons_inframe_del+'\t'+cons_inframe_ins+'\t'+disp_inframe_del+'\t'+disp_inframe_ins+'\t'+downstream_genevar+'\t'+frameshift+'\t'+intron+'\t'+nctxpt+'\t'+splice_acc+'\t'+splice_don+'\t'+splice_region+'\t'+stop_retained+'\t'+upstream_genevar+'\t'+Type+'\t'+kw+'\t'+evi+'\n')
                #of.write(k+'\t'+rapid+'\t'+detail+'\t'+'No Chr\tNo start\tNo end\tNo name\tNo description\t'+function+'\t'+ins+'\t'+dels+'\t'+exon_loss+'\t'+start_lost+'\t'+stop_gain+'\t'+stop_lost+'\n')
'''
