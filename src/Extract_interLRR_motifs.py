# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 16:56:07 2019

@author: gottince
"""

#========================================================
# This program extract LRR motifs from HMMsearch results 
# file
#========================================================

# Packages and Classes
#=====================================

import argparse
from Bio import SearchIO
from Bio import SeqIO
import Class_HMMhit as HH
import hmmsearch_results as HR

# Arguments
#=====================================

parser = argparse.ArgumentParser(description='This program extract LRR motifs from interLRR sequences.',prog='Extract_interLRR_motifs.py')

parser.add_argument('-s', '--sequences', nargs=1, type=str, required=True, help='interLRR sequences in fasta format')

parser.add_argument('-t', '--tblfile', nargs=1, type=str, required=True, help='hmmsearch domtblout file')

parser.add_argument('-o', '--output', nargs=1, type=str, required=True, help='name of results file')

args = parser.parse_args()


# Script
#=====================================

## 1 . Importing data
##--------------------

## a) interLRR sequences fasta file       
fastaFile = args.sequences[0]
interLRR_sequences = SeqIO.to_dict(SeqIO.parse(fastaFile,"fasta"))

## b) HMMsearch results file
tblFile = args.tblfile[0]
myRes=SearchIO.read(tblFile,"hmmsearch3-domtab")

## 2 . Processing data
##---------------------

## a) File to print annotation results
resFile = args.output[0]

myFile=open(resFile,"w")
myFile.write("Protein;Domain;Index;Start;End;Length;Sequence\n")
myFile.close()

## b) Processed Hit
for hit in myRes :
    #filtering data
    hit=HR.hsp_start_filter(hit,8)
    if(hit):
        hit=HR.hsp_end_filter(hit,10)
    
    ##Il reste des Hit?
    if(hit):
        #save start postion of interLRR in protein and protein id
        pId=hit.id.split(';')[0]
        sPos=int(hit.id.split(';')[1])
        #Object hmmHit
        myHit=HH.hmmHit()
        myHit.complete_with_hmmsearch(hit, myRes.id, myRes.seq_len, str(interLRR_sequences[hit.id].seq))
        myHit.processed_interLRR()
        myHit.save_hit_interLRR(resFile,pId,sPos)

#========================================================
#                    END of SCRIPT
#======================================================== 
