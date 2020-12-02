# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 17:11:37 2019

@author: gottince
"""

#==================================================
# This program concatenate all LRR motifs 
#==================================================

# Packages and Classes
#=====================================

import argparse
import os
import re
from Bio import SeqIO
import ProteinDataClass as PDC
import csv

# Arguments
#=====================================

parser = argparse.ArgumentParser(description='This program concatenate all LRR motifs from several sources without duplication.',prog='Concat_all_motifs.py')

parser.add_argument('-s', '--sequences', nargs=1, type=str, required=True, help='protein sequences in fasta format')

parser.add_argument('-o', '--output', nargs=1, type=str, required=True, help='name of results file')

parser.add_argument('-d', '--directory', nargs=1, type=str, required=True, help='directory of extract LRR files')

args = parser.parse_args()

# Script
#=====================================

## 1 . Importing data
##--------------------   

## a) Protein sequences fasta file
fastaFile = args.sequences[0]
protein_sequences = SeqIO.to_dict(SeqIO.parse(fastaFile,"fasta"))

## b) create instance proteome
DATA=PDC.Proteome("Proteome")

## c) csv files
myDir=args.directory[0]
fileList=os.listdir(myDir)
           

## Fichier LRR
for files in fileList :
    if(re.search(".csv",files)!=None):
        #csvFile="%s/%s"%(myDir,files)
        #csv_file=open(csvFile, mode='r')
        csv_file=open(files, mode='r')
        csv_reader = csv.reader(csv_file, delimiter=';')
        line_count = 0
        for row in csv_reader:
            #skip header
            if line_count == 0 :
                line_count += 1
            else :
                if (row[0] not in DATA.proteins) :
                    ##Creation de la proteine dans le proteome
                    DATA.add_protein(row[0],0,DATA.id)
                # Ajout des motifs
                if (row[1]!="interLRR") :
                    DATA.proteins[row[0]].add_motif(row[1],int(row[3]),int(row[4]),float(row[6]))
        csv_file.close()


## 2 . Processing data
##---------------------
# Pour chaque proteine, on ordonne les motifs en fonction de la position de depart
# puis on supprime les elements ayant la meme position
for prot in DATA.proteins :  
    DATA.proteins[prot].rm_duplicate()
    DATA.proteins[prot].order_motifs()
    DATA.proteins[prot].correct_pos()
    DATA.proteins[prot].exclude_blast_outlier()
    DATA.proteins[prot].order_motifs()
    DATA.proteins[prot].extract_inter_regions()

## 3 . Save data
##---------------

resFile = args.output[0]
DATA.save_to_file(resFile,protein_sequences)

#========================================================
#                    END of SCRIPT
#========================================================
