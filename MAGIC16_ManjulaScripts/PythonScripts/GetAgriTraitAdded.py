import os,sys,csv,random,re,subprocess
#sys.setdefaultencoding("ISO-8859-1")

## Usage:
## thimmamp@KW61043:~/PythonPrograms/utilities$ python GetAgriTraitAdded.py ~/MAGIC16/annotation/geneKeyword.table.txt ~/MAGIC16/indels_results/genome1_Core_Dispensible_Specific_temp > ~/MAGIC16/indels_results/genome1_Core_Dispensible_Specific_withKeywordEvidences.txt
## thimmamp@KW61043:~/PythonPrograms/utilities$ python GetAgriTraitAdded.py ~/MAGIC16/annotation/Grain_clonedgenes.txt ~/MAGIC16/indels_results/genome1_Core_Dispensible_Specific_temp > ~/MAGIC16/indels_results/genome1_Core_Dispensible_Specific_withGrainTrait.txt

## Usage for filtering impacted variant events
##python GetAgriTraitAdded.py ~/MAGIC16/indels_results/genome1_Core_Dispensible_Specific_withAnnotation_KW_Evidences.txt ~/MAGIC16/indels_results/genome1_Core_Dispensible_Specific_withAnnotation_KW_Evidences_filteredforImpactVariantevents.txt

## Usage for combining indels, blast hits with nr, keywords, and evidences
## thimmamp@KW61043:~/PythonPrograms/utilities$ python GetAgriTraitAdded.py ~/MAGIC16/indels_results/genome33_specefic_genic_InDelsStats_Details_testing.txt ~/MAGIC16/indels_results/IR64SpecificVsnr_blastx_percidentity_80_evalue_1e-3_target10.results.txt ~/MAGIC16/annotation/geneKeyword.table.txt ~/MAGIC16/indels_results/temp_report.txt 

##########
'''
exon_loss
12 start_lost
13 stop_gained
14 stop_lost
15 start_retained
16 threeputr_trunc
17 threeputr_variant
18 fiveputr_trunc
19 fiveputr_variant
20 cons_inframe_del
21 cons_inframe_ins
22 disp_inframe_del
23 disp_inframe_ins
24 downstream_genevar
25 frameshift
26 intron
27 nctxpt
28 splice_acc
29 splice_don
30 splice_region
31 stop_retained
32 upstream_genevar

{'NiporGrameneId': 'OsNip_01G000100', 'RAPId': 'Os01g0137950', 'INDELDetails': 'Chr01_57787_AC_A;AC;A;bi_allele;del;82;;Chr01_57989_GTGC_G;GTGC;G;bi_allele;del;6;;Chr01_58140_CGCTA_C;CGCTA;C;bi_allele;del;79;;Chr01_58172_C_CTAACAACAT;C;CTAACAACAT;bi_allele;ins;29;;Chr01_58200_A_ATG;A;ATG;bi_allele;ins;52;;Chr01_58218_T_TG;T;TG;bi_allele;ins;359;;Chr01_58247_A_AT;A;AT;bi_allele;ins;876', 'Chromosome': 'Chr01', 'start': '2037202', 'end': '2038198', 'Name': 'OsRLCK23', 'Description': 'Similar to predicted protein', 'FunctionalAnnotation': 'NA', 'INS': '4', 'DELS': '3', 'exon_loss': '0', 'start_lost': '0', 'stop_gained': '1', 'stop_lost': '0', 'start_retained': '0', 'threeputr_trunc': '0', 'threeputr_variant': '1', 'fiveputr_trunc': '0', 'fiveputr_variant': '0', 'cons_inframe_del': '0', 'cons_inframe_ins': '0', 'disp_inframe_del': '0', 'disp_inframe_ins': '1', 'downstream_genevar': '114', 'frameshift': '2', 'intron': '2', 'nctxpt': '0', 'splice_acc': '0', 'splice_don': '0', 'splice_region': '0', 'stop_retained': '0', 'upstream_genevar': '61', 'Type': 'Core', 'Keywords': 'NA', 'Evidences': 'NA'}

'''
##########

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

def ReadOmicsAnnotOutput(infile):
    outdict={}
    with open(infile, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        ## keys in the file ['Status', 'Tags', 'SeqName', 'Description', 'Length', '#Hits', 'e-Value', 'sim mean', '#GO', 'GO IDs', 'GO Names', 'Enzyme Codes', 'Enzyme Names', 'InterPro IDs', 'InterPro GO IDs', 'InterPro GO Names']
        for row in reader:
            #print(row.keys())
            gonames='NA'
            mf_str='NA'
            if row['GO Names']:
                gonames = row['GO Names']
                if gonames:
                    mf_list = set()
                    for item in gonames.split('; '):
                        if 'F:' in item: ## get only molecular function out of BP, MF, CC
                            tmp = item
                            tmp = tmp.replace('F:','')
                            mf_list.add(tmp)
                            #mf_str+=item.replace('F:','')+';'
                    mf_str = ';'.join(mf_list)        
            #outdict[row['SeqName']] = gonames
            outdict[row['SeqName']] = mf_str
    return outdict        


def ReadAgriTrait(infile):
    outdict={}
    with open(infile, 'r') as f:
        for line in f:
            if not line.startswith("RAP"):
                key, msu, trait = line.rstrip('\n').split('\t')
                outdict[key]=msu+";;"+trait
                #print(key, msu, trait)
    return outdict            


def GetGenicIndelDetailstoDict(infile):
    outdict = {}
    with open(infile, 'r') as f:
        for line in f:
            if not line.startswith('Geneid'):
                tokens = line.rstrip('\n').split('\t')
                outdict[tokens[0]] = '\t'.join(tokens[1:])
    return outdict            

def GetBlastOutputtoDict(blastoutfile):
    outdict = {}
    pat = r'(?<=\[).+?(?=\])' 
    with open(blastoutfile, 'r') as f:
        for line in f:
            tokens = line.rstrip('\n').split('\t')
            #key = tokens[0].split('::')[0]
            key = tokens[0]
            details = tokens[5]
            hit_protein=tokens[3]
            rapids = []
            hit_species = set(re.findall(pat, details))
            hit_species = ','.join(hit_species)
            for item in details.split(' '):
                if item.startswith('Os'):
                   item = item.replace(',','')
                   if item not in rapids:
                       rapids.append(item)
            rapid='NA'
            if rapids:
                rapid = ','.join(rapids)
            if key in outdict.keys():
                outdict[key]['details'].append(details)
                outdict[key]['rapids'].append(rapid)
                if hit_protein not in outdict[key]['hit_protein']:
                   outdict[key]['hit_protein'].append(hit_protein)
                if hit_species not in outdict[key]['hit_species']:   
                   outdict[key]['hit_species'].append(hit_species)
            else:
                outdict[key] = {}
                outdict[key]['details'] = [details]
                outdict[key]['rapids'] = [rapid]
                outdict[key]['hit_protein'] = [hit_protein]
                outdict[key]['hit_species'] = [hit_species]
            #outdict[key] = details+';;'+rapid+';;'+hit_protein+';;'+hit_species
    return outdict            


def FilterGeneswithINDELSVarimpactKWEvidences(infile, outfile):
    with open(outfile, 'w') as f1:
        f1.write("NiporGrameneId\tRAPId\tName\tDescription\tFunctionalAnnotation\tINS\tDELS\tType\tKeywords\tEvidences\tvarimpact\n")
        with open(infile, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                #print(row)
                varimpact=''
                for num, name in enumerate(row.keys()):
                    if num >= 11 and num <= 32:
                        if int(row[name]) > 0:
                            varimpact += name+'='+row[name]+';'
                if varimpact=='':            
                    varimpact = 'NA'
                f1.write(row['NiporGrameneId']+'\t'+row['RAPId']+'\t'+row['Name']+'\t'+row['Description']+'\t'+row['FunctionalAnnotation']+'\t'+row['INS']+'\t'+row['DELS']+'\t'+row['Type']+'\t'+row['Keywords']+'\t'+row['Evidences']+'\t'+varimpact+'\n')

##~/MAGIC16/indels_results$ less IR64SpecificVsnr_blastx_percidentity_80_evalue_1e-3_new.results.txt
#(base) thimmamp@KW61043:~/MAGIC16/indels_results$ head genome33_specefic_genic_InDelsStats_Details_testing.txt


## call functions
#FilterGeneswithINDELSVarimpactKWEvidences(sys.argv[1], sys.argv[2])

## go through loop of indel detail files, blastout
#:~/MAGIC16/indels_results$ ls -altr Nipponbare_AbsentGenes_from_genome*_mode_duplication*.txt
directory = "/home/thimmamp/MAGIC16/indels_results/"
blast2details = GetBlastOutputtoDict(directory+"IRGSPabsentGenesVsnr_blastx_percidentity_95_evalue_1e-3_target10.results.txt")
irgsp2kwevi = GetKeywordEvidenceForRAPId("/home/thimmamp/MAGIC16/annotation/geneKeyword.table.txt")

outfile = directory+"GenesAbsentinNipponbare_presentin_othergenomes_indeldet_blasthits_irgsp_kwevi_withUniqueGoMolecularFunction.txt"
omicsout = ReadOmicsAnnotOutput(directory+"Nipponbare_Absent_Genes_GOFuncAnnotation.txt")

with open(outfile, 'w') as of:
    of.write("Geneid\tINDEL_Details\tINS\tDELS\tBlastquery\tBlastdetails\tRAPid\tHitProtein\tHitSpecies\tkeywords\tevidences\tGONames\n")
    for name in os.listdir(directory):
        if 'Nipponbare_AbsentGenes_from_genome' in name and '_mode_duplication' in name:
            print("processing file: ", name) 
            genome=name.split('_')[3]
            genomeindelfile = directory+"Nipponbare_AbsentGenes_from_"+genome+"_with_mode_duplication_genic_InDelsStats_Details_testing.txt"
            genedetails = GetGenicIndelDetailstoDict(genomeindelfile)
            for k1, v1 in genedetails.items():
                geneid=''
                if '.' in k1:
                    geneid = k1.split('.')[-1]
                elif '_' in k1:
                    geneid = k1.split('_')[-1]
                for k,v in blast2details.items():
                    outStr='' 
                    blastgeneid=''
                    blastgenome=''
                    if '.' in k:
                        blastgeneid = k.split('.')[-1]
                        blastgenome = k.split('_')[0]
                    else:
                        blastgenome = k.split('_')[0]
                        blastgeneid = k.split('_')[-1]

                    if genome==blastgenome and geneid==blastgeneid:
                        blastdet = ';;'.join(blast2details[k]['details'])
                        hit_protein = ';;'.join(blast2details[k]['hit_protein'])
                        hit_species = ';;'.join(blast2details[k]['hit_species'])
                        #rapid = str( set( blast2details[k]['rapids'] ) )        
                        rapid = ';;'.join(set( blast2details[k]['rapids'] ) )        
                        kw='NA';evi='NA'
                        if rapid in irgsp2kwevi.keys():
                            kw = ','.join(irgsp2kwevi[rapid]['keywords'])
                            evi = ','.join(irgsp2kwevi[rapid]['evidences'])
                        outStr=k1+'\t'+v1+'\t'+k+'\t'+ blastdet+'\t'+ rapid +'\t'+ hit_protein +'\t'+ hit_species+'\t'+kw+'\t'+evi
                        #print(outStr)
                        if k in omicsout.keys():
                            outStr+='\t'+omicsout[k]+'\n'
                            #print(outStr, k, omicsout[k])
                            of.write(outStr)




#######################        
'''
genedetails = GetGenicIndelDetailstoDict(sys.argv[1]) ## read genome33 Specific genes and their indels
blast2details = GetBlastOutputtoDict(sys.argv[2]) ## read genome33 Specific genes blast with nr and their details
irgsp2kwevi = GetKeywordEvidenceForRAPId(sys.argv[3])

## Now combine the above dicts for analysis
with open(sys.argv[4], 'w') as of:
    of.write("Geneid\tIndelDetails\tINS\tDELS\tBlasthit\tHitProtein\tHitSpecies\tkeywords\tevidences\n")
    for k,v in genedetails.items():
        if k in blast2details.keys():
            blastdet = ','.join(blast2details[k]['details'])
            rapid = str( set( blast2details[k]['rapids'] ) )        
            hit_protein = ','.join(blast2details[k]['hit_protein'])
            hit_species = ','.join(blast2details[k]['hit_species'])
            #blastdet, rapid, hit_protein, hit_species = blast2details[k].split(';;')
            kw='NA';evi='NA'
            if rapid in irgsp2kwevi.keys():
                kw = ','.join(irgsp2kwevi[rapid]['keywords'])
                evi = ','.join(irgsp2kwevi[rapid]['evidences'])
            of.write(k+'\t'+v+'\t'+blastdet+'\t'+hit_protein+'\t'+hit_species+'\t'+kw+'\t'+evi+'\n')    
            #print(k+'\t'+v+'\t'+blastdet+'\t'+hit_protein+'\t'+hit_species+'\t'+kw+'\t'+evi+'\n')    
        else:
            of.write(k+'\t'+v+'\tNo blast hit\tNA\tNA\tNA\tNA\n')
            #print(k+'\t'+v+'\tNo blast hit\tNA\tNA\tNA\tNA\n')
'''

'''
#rapid2trait = ReadAgriTrait(sys.argv[1])
rap2kwevi = GetKeywordEvidenceForRAPId(sys.argv[1])
## Map rap id in the input file and get the trait
with open(sys.argv[2], 'r') as varf:
    for line in varf:
        tokens = line.rstrip('\n').split('\t')
        if "Os" in tokens[1]:
            rapid = tokens[1]
            #if rapid in rapid2trait.keys():
            if rapid in rap2kwevi.keys():
                traitkw = ','.join(rap2kwevi[rapid]['keywords'])
                traitevi = ','.join(rap2kwevi[rapid]['evidences'])
                out = '\t'.join(tokens)+'\t'+traitkw+'\t'+traitevi
                print(out)
            else:
                out = '\t'.join(tokens)+'\tNA\tNA'
                print(out)
        else:
            out = '\t'.join(tokens)+'\tNA\tNA'
            print(out)


#for k,v in rap2kwevi.items():
#    outStr=k
#    print(k,'\t',','.join(rap2kwevi[k]['keywords']),'\t', ','.join(rap2kwevi[k]['evidences']))
    #for k1 in rap2kwevi[k].keys():
    #    print(k, k1, rap2kwevi[k][k1])
'''    
