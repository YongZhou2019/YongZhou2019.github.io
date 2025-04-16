import os,sys,csv

## how to use
## python GetIncompleteSplits.py split_mod.txt missing_splits > incomplete_split.txt
## python GetIncompleteSplits.py split.txt tobe_tbi > split.pending.txt
def GetIncompleteSplits(mastersplit, completedsplit):
    #outdict={}
    with open(mastersplit, 'r') as f:
         with open(completedsplit, 'r') as inputf:
             for line in f:
                 temp = line.rstrip('\n').split('\t')
                 print(temp)
                 #print(chro, chunk)
                 for line1 in inputf:
                     tokens = line1.rstrip('\n').split('\t')
                     print(tokens)
                     #if tokens[0]==temp[0] and tokens[1]==temp[1]:
                     #   print(line1)
def GetTbi2Dict(infile):
    outdict={}
    with open(infile, 'r') as f:
        for line in f:
            tokens = line.rstrip('\n').split("\t") #. for tbi files
            k = tokens[0]+'_'+tokens[1] # 1 and 2 for tbi files
            v = 1
            outdict[k] = v
    return outdict            


def GetSplit2Dict(infile):
    outdict={}
    with open(infile, 'r') as f:
        for line in f:
            tokens = line.rstrip('\n').split('\t')
            if len(tokens)==4:
                k = tokens[0]+'_'+tokens[1]
                v = '\t'.join(tokens[2:])
                outdict[k] = v
            else:
                k = tokens[0]+'_'+tokens[1]
                v = 1
                outdict[k] = v
    return outdict            

master = GetSplit2Dict(sys.argv[1])
incomplete = GetSplit2Dict(sys.argv[2])
#incomplete = GetTbi2Dict(sys.argv[2])

#print("Number of entries in master : ", len(master.keys()),'\n') 
#print("Number of entries in processed : ", len(incomplete.keys()),'\n') 

for k,v in master.items():
    if not k in incomplete.keys():
        k = k.replace("_",'\t')
        print(k+'\t'+v)

#GetIncompleteSplits(sys.argv[1], sys.argv[2])


