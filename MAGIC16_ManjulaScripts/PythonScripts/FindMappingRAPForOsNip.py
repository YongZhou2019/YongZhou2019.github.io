import os,sys,csv,random,re,subprocess
from pathlib import Path
from subprocess import Popen, PIPE

## Usage:
## thimmamp@kw61043:~/PythonPrograms/utilities$ python FindMappingRAPForOsNip.py OsNip_geneids_fromdispensabletable.txt ~/MAGIC16/annotation/IRGSPRAP_Nipponbare.one2one

def ReadFile2Dict(infile):
    outDict={}
    with open(infile, 'r') as f:
        for line in f:
            line=line.rstrip("\n")
            if ',' in line: # there could be multiple Nip ids in a line
                tokens = line.split(",")
                for item in tokens:
                    outDict[item]=1
            else:
                outDict[line]=1
    return outDict


def File2Dict(infile):
    outDict={}
    with open(infile, 'r') as f:
        for line in f:
            tokens=line.rstrip('\n').split('\t')
            outDict[tokens[0]]=tokens[1]
    return outDict  

def FindQueryinTarget(querydict, targetdict):
    for k,v in querydict.items():
        if k in targetdict.keys():
            print(k, targetdict[k])
        else:
            print(k, "No RAPid")

## use functions
nipid_dict = ReadFile2Dict(sys.argv[1])
targetdict = File2Dict(sys.argv[2])
FindQueryinTarget(nipid_dict, targetdict)

#dispNip2RAP=FindQueryinTarget(nipid_dict, sys.argv[2])
#print("############################")
#print("Unmapped Nip ids")
#for k,v in dispNip2RAP.items():
#    print(k,'\t',v)
#print("############################")
