# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 15:05:35 2019

@author: gottin celine
"""

#=================================================
#  Project : LRR proteins annotation
#
# Class for sequence and domain data
#=================================================

import re

#__________________________________________________
#                   Sequence Class
#__________________________________________________

class Sequence :

    def __init__(self,Identifier,Size=0) :
        self.id=Identifier
        self.size=Size

    def get_id(self) :
        return self.id

    def set_id(self, NewId) :
        self.id=NewId

    def get_size(self) :
        return self.size

    def set_size(self,Size) :
        self.size=Size

#__________________________________________________
#                   Protein Class
#__________________________________________________

class Protein(Sequence) :

    def __init__(self,Identifier,Size=0,Organism="Unknown",Chromosome="Unknown") :
        super().__init__(Identifier,Size)
        self.organism=Organism
        self.chromosome=Chromosome
        self.motifs=[]
        self.intermotifs=[]
        #self.domains=[]

    def __str__(self) :
        line=" Protein id : %s \n Size : %i aa \n Organism : %s \n Chromosome : %s" % (self.id,self.size,self.organism,self.chromosome)
        return line

## Accessors and modificators

    def get_organism(self) :
        return self.organism

    def set_organism(self,Organism) :
        self.organism=Organism

    def get_chromosome(self) :
        return self.chromosome

    def set_chromosome(self,Chrm) :
        self.chromosome=Chrm

## Methods

    def add_motif(self,motif,Start=0,End=0,Eval=0) :
        index=len(self.motifs)+1
        self.motifs.append(Motif(motif,index,Start,End,Eval))


    def add_interMotif(self, motif,Start=0,End=0):
        index=len(self.intermotifs)+1
        self.intermotifs.append(Motif(motif,index,Start,End))


    def insert_motif(self,motif,index,Start=0,End=0,Eval=0) :
        self.motifs.insert(index,Motif(motif,index,Start,End,Eval))


    def order_motifs(self) :
        ##Sort motifs according to their start positions
        ##then update indexes and handle duplicate/overlap

        pos=0
        if (len(self.motifs)>1) :
            for idx in range(1,len(self.motifs)) :
                if (self.motifs[idx].start<self.motifs[idx-1].start) :
                    elm=self.motifs[idx]
                    self.motifs.pop(idx)
                    pos=0
                    while (self.motifs[pos].start<=elm.start) :
                        pos=pos+1
                    self.motifs.insert(pos,elm)

        for i in range(len(self.motifs)) :
            self.motifs[i].set_index(i+1)


    def rm_duplicate(self) :
        #For every LRR motif of the protein, if two motifs have the same
        #start position --> suppr one of them
        #Suppr LRR motifs that are inside another

        allStart=[]
        allEnd=[]
        for i in range(len(self.motifs)):
            allStart.append(self.motifs[i].start)
            allEnd.append(self.motifs[i].end)

        # Check redondance
        for i in range(len(self.motifs)-1,-1,-1):
            if(allStart.count(allStart[i])>1) : ## if several elm with same start position
                #index of first duplicated elm
                id1=allStart.index(allStart[i])
                #index of second duplicated elm
                id2=allStart.index(allStart[i],id1+1)

                ##if f-box vs other -->priority to other
                if(re.search("LRR_Fbox",self.motifs[id2].type)!=None):
                    ##suppr second elm
                    self.rm_motif(id2)
                    allStart.pop(id2)
                    allEnd.pop(id2)
                elif(re.search("LRR_Fbox",self.motifs[id1].type)!=None):
                    # Suppr first elm
                    self.rm_motif(id1)
                    allStart.pop(id1)
                    allEnd.pop(id1)
                elif(self.motifs[id1].eval<=self.motifs[id2].eval):
                    #Suppr second elm
                    self.rm_motif(id2)
                    allStart.pop(id2)
                    allEnd.pop(id2)
                else :
                    # Suppr first elm
                    self.rm_motif(id1)
                    allStart.pop(id1)
                    allEnd.pop(id1)

        # Check overlapping
        exclude=[]
        for i in range(len(allStart)-1,0,-1):
            for j in range(0,i):
                #if j inside i or i inside j
                if((allStart[i]<=allStart[j] and allEnd[i]>=allEnd[j]) or (allStart[i]>=allStart[j] and allEnd[i]<=allEnd[j])):
                    # Check eval
                    if(self.motifs[i].eval <= self.motifs[j].eval):
                        if not j in exclude :
                            exclude.append(j)
                    else:
                        if not i in exclude :
                            exclude.append(i)
                # if i inside j
                #if(allStart[i]>allStart[j] and allEnd[i]<allEnd[j]):
                    #if not i in exclude:
                        #exclude.append(i)
                # if overlap
                elif((allStart[i]>allStart[j] and allStart[i]<allEnd[j]) or (allStart[j]>allStart[i] and allStart[j]<allEnd[i])):
                    if(self.motifs[i].eval <= self.motifs[j].eval):
                        if not j in exclude:
                            exclude.append(j)
                    else:
                        if not i in exclude:
                            exclude.append(i)
        # To avoid "out of range", supr index in revrse order
        exclude.sort(reverse=True)
        for ind in exclude :
            self.rm_motif(ind)


    def correct_pos(self) :
        ## Corretc motif position if overlapping
        for i in range(len(self.motifs)-1) :
            if(self.motifs[i].end>=self.motifs[i+1].start):
               self.motifs[i].set_end(self.motifs[i+1].start-1)


    def rm_motif(self,index) :
        self.motifs.pop(index)


    def exclude_blast_outlier(self) :
        ##First motif : if BLAST and over 10 AA from next motif -> remove
        if(self.motifs[0].type=="BLAST" and self.motifs[0].end+10<self.motifs[1].start) :
            self.rm_motif(0)

         ##Last motif
        if(self.motifs[-1].type=="BLAST" and self.motifs[-1].start-10>self.motifs[-2].end) :
            self.rm_motif(-1)

    def find_motif_in_range(self,start,end) :
        "TODO"


    def start_to_list(self) :
        L=[]
        for elm in self.motifs :
            L.append(elm.start)
        return L


    def end_to_list(self) :
        L=[]
        for elm in self.motifs :
            L.append(elm.end)
        return L


    def extract_inter_regions(self) :
        for i in range(len(self.motifs)-1) :
            # InterLRR of 3 aa at least
            if(self.motifs[i+1].start > self.motifs[i].end+3):
                self.add_interMotif("interLRR",self.motifs[i].end+1,self.motifs[i+1].start-1)


#__________________________________________________
#                   Domain Class
#__________________________________________________

class Domain(Sequence) :

    def  __init__(self,Name,Start,End) :
        super().__init__(Name,(End-Start+1))
        self.start=Start
        self.end=End
#__________________________________________________
#                   Motif Class
#__________________________________________________

class Motif() :

    def __init__(self,Motif,index,Start,End,Score=0) :
        self.type=Motif
        self.index=index
        self.size=End-Start+1
        self.start=Start
        self.end=End
        self.eval=Score
        #!! check if size and Start/End are coherent
        #if self.start!=0 and self.end!=0 :
         #   if self.size==0 :
         #       self.size=(self.end-self.start+1)
         #   else:
          #      if self.size!=(self.end-self.start+1) :
          #          print("error in domain size")

    def __str__(self) :
        return " Motif Type : %s \n Index : %i \n Size : %i aa \n Start : %i \n End : %i" % (self.type,self.index,self.size,self.start,self.end)

    def get_type(self) :
        return self.type

    def set_type(self,NewType) :
        self.type=NewType

    def get_start(self) :
        return self.start

    def set_start(self,NewStart) :
        self.start=NewStart

    def get_end(self) :
        return self.end

    def set_end(self,NewEnd) :
        self.end=NewEnd

    def get_index(self) :
        return self.index

    def set_index(self,NewIndex) :
        self.index=NewIndex

#__________________________________________________
#               Class Proteome
#__________________________________________________

class Proteome:
    def __init__(self,ID,Organism="unknown") :
        self.id=ID
        self.organism=Organism
        self.size=0
        self.proteins={}

    def add_protein(self,Identifier,Size=0,Chromosome="Unknown"):
        self.proteins[Identifier]=Protein(Identifier,Size,self.organism,Chromosome)
        self.size=len(self.proteins)

    def save_to_file(self,filename,sequences) :
        file=open(filename,"w")

        file.write("Protein;Domain;Index;Start;End;Length;Eval;Sequence\n")
        for prot in self.proteins :
            protSeq=sequences[prot].seq

            for mot in self.proteins[prot].motifs :
                file.write("%s;%s;%i;%i;%i;%i;%f;%s\n" % (prot,mot.type,mot.index,mot.start,mot.end,mot.end-mot.start+1,mot.eval,protSeq[mot.start-1:mot.end]))

            for imot in self.proteins[prot].intermotifs :
                file.write("%s;%s;%i;%i;%i;%i;%f;%s\n" % (prot,imot.type,imot.index,imot.start,imot.end,imot.end-imot.start+1,imot.eval,protSeq[imot.start-1:imot.end]))

        file.close()
