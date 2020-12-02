# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 15:56:07 2019

@author: gottince
"""

#========================================================
# This program extract LRR motifs from blast results 
# file
#========================================================

# Packages and Classes
#=====================================

import argparse
#from Bio import SearchIO
from Bio import SeqIO
import Class_HMMhit as HH
#import hmmsearch_results as HR

# Arguments
#=====================================

parser = argparse.ArgumentParser(description='This program extract LRR motifs from interLRR sequences.',prog='Extract_interLRR_motifs.py')

parser.add_argument('-s', '--sequences', nargs=1, type=str, required=True, help='interLRR sequences in fasta format')

parser.add_argument('-b', '--blastfile', nargs=1, type=str, required=True, help='blast results file')

parser.add_argument('-o', '--output', nargs=1, type=str, required=True, help='name of results file')

args = parser.parse_args()


# Script
#=====================================

## 1 . Importing data
##--------------------

## a) interLRR sequences fasta file       
fastaFile = args.sequences[0]
sequences = SeqIO.to_dict(SeqIO.parse(fastaFile,"fasta"))

## b) HMMsearch results file
blastFile = args.blastfile[0]

listElm={}

with open(blastFile) as f:
    for line in f:
        pname=line.split(';')[3]
        #txt=line.split(';')[3:6]
        #pname=txt[0]+";"+txt[1]+";"+txt[2]
        if(pname not in listElm.keys()):
            listElm[pname]=HH.hmmHit()
        
        listElm[pname].complete_with_blastline(line,sequences[pname].seq)

f.close()

## 2 . Processing data
##---------------------

## a) File to print annotation results
resFile = args.output[0]

myFile=open(resFile,"w")
myFile.write("Protein;Domain;Index;Start;End;Length;Sequence\n")
myFile.close()

for elm in listElm.keys() :
    listElm[elm].processed_blastRes()
    listElm[elm].save_hit_blast(resFile,elm)


#========================================================
#                    END of SCRIPT
#======================================================== 
