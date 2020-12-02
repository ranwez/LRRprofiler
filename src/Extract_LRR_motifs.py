# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 16:45:13 2019

@author: gottince
"""

#========================================================
# This program extract LRR motifs from HMMsearch results 
# file for LRR_CC and NLR
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

parser = argparse.ArgumentParser(description='This program extract LRR motifs from HMMsearch results.',prog='Extract_LRR_motifs.py')

parser.add_argument("-s", "--sequences", nargs=1, type=str, required=True, help='protein sequences in fasta format')

parser.add_argument("-t", "--tblfile", nargs=1, type=str, required=True, help='hmmsearch domtblout file')

parser.add_argument("-o", "--output", nargs=1, type=str, required=True, help='name of results file')

args = parser.parse_args()

# Script
#=====================================

## 1 . Importing data
##--------------------

## a) Protein sequences fasta file
fastaFile = args.sequences[0]
protein_sequences = SeqIO.to_dict(SeqIO.parse(fastaFile,"fasta"))

## b) HMMsearch results file
tblFile = args.tblfile[0]
myRes=SearchIO.read(tblFile,"hmmsearch3-domtab")


## 2 . Processing data
##---------------------

## a) File to print annotation results
resFile = args.output[0]

myFile=open(resFile,"w")
myFile.write("Protein;Domain;Index;Start;End;Length;eval;Sequence\n")
myFile.close()

## b) Processed Hit
for hit in myRes :
    #if hit.evalue < 10 :
    ## filtering hsp
    hit=HR.hsp_start_filter(hit, 8)
    if(hit):
        hit=HR.hsp_end_filter(hit, 10)
    ## processed
    if(hit):
        myHit=HH.hmmHit()
        myHit.complete_with_hmmsearch(hit, myRes.id, myRes.seq_len, str(protein_sequences[hit.id].seq))
        myHit.processed_data()
        myHit.save_hit_to_file(resFile)
        
        
#========================================================
#                    END of SCRIPT
#========================================================        
