import os,sys,csv,random,re,subprocess
#sys.setdefaultencoding("ISO-8859-1")

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

## call functions
rap2kwevi = GetKeywordEvidenceForRAPId(sys.argv[1])

for k,v in rap2kwevi.items():
    outStr=k
    print(k,'\t',','.join(rap2kwevi[k]['keywords']),'\t', ','.join(rap2kwevi[k]['evidences']))
    #for k1 in rap2kwevi[k].keys():
    #    print(k, k1, rap2kwevi[k][k1])

