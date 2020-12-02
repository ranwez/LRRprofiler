# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 15:26:58 2019

@author: gottince
"""
#=================================================
#  Project : LRR proteins annotation
#
# Class for HMMsearch result hit processed
#=================================================

#import re

class hmmHit() :
    def __init__(self) :
        ## Récupération des info sur un Hit (protéine) HMMsearch
        self.id="" ## id protéine
        self.seq="" ## sequence protéine
        self.profile="" ## profile name
        self.profile_len="" ## taille motif
        self.score=[]
        self.eval=[]
        self.startHsps=[]
        self.endHsps=[]
        self.startQuery=[]
        self.endQuery=[]
        self.interMotifLen=[]
        self.interMotifPos=[]
        self.typeMotif=[]
        self.limitSup=[] ##pour blast : borne sup et inf du interLRR à ne pas dépasser
        self.limitInf=[]
     
        
    def complete_with_hmmsearch(self, hmmsearchHit, profilename, profLen, sequence) :
        for elm in hmmsearchHit.hsps :
            #if(elm.evalue<1.0): ## filtre sur evalue du motif
            self.id=hmmsearchHit.id
            self.seq=sequence
            self.profile=profilename
            self.profile_len=profLen
            self.score.append(elm.bitscore)
            self.eval.append(elm.evalue)
            self.startHsps.append(elm.hit_start+1)
            self.endHsps.append(elm.hit_end)
            self.startQuery.append(elm.query_start+1)
            self.endQuery.append(elm.query_end)
            self.typeMotif.append(self.profile)
  
  
    def complete_with_blastline(self, resline, sequence) :
        liste=resline.split(';')
        self.id=liste[3]
        self.limitSup.append(int(liste[5]))
        self.limitInf.append(int(liste[4]))
        self.seq=sequence
        self.profile="BLAST"
        self.profile_len=int(liste[16])
        self.score.append(float(liste[10]))
        self.eval.append(float(liste[11]))
        self.startHsps.append(int(liste[4])+int(liste[14])-1)
        self.endHsps.append(int(liste[4])+int(liste[15])-1)
        self.startQuery.append(int(liste[12]))
        self.endQuery.append(int(liste[13]))
        self.typeMotif.append(self.profile)
        
        
    def extend_nter(self) :
        # Etendre les LRR incomplets au début du motif
        self.startHsps[0]=self.startHsps[0]-(self.startQuery[0]-1)  # Pour le premier motif   
        ##si out of range (debut motif avant début sequence) debut motif=1
        if (self.startHsps[0]<1) :
            self.startHsps[0]=1
        # Pour les autres motifs
        for lrr in range(1,len(self.startHsps)) :
            # Si debut du motif i après/chevauche moins de deux résidus du motif i-1
            # Debut motif i = début profil
            if (self.startHsps[lrr]-(self.startQuery[lrr]-1)) > self.endHsps[lrr-1]-2 :
                self.startHsps[lrr]=self.startHsps[lrr]-(self.startQuery[lrr]-1)
            else :
                self.startHsps[lrr]=self.endHsps[lrr-1]+1   
                
                
                
    def extend_cter(self) :
        # Etendre la fin du LRR jusqu'a fin du motifs (+2 aa si motif suivant juste derrière) 
        #(ou jusqu'à motif suivant si chevauchement)
        for lrr in range(len(self.endHsps)-1) :
            ## pos fin du motif si profil complet
            endMot=self.endHsps[lrr]+(self.profile_len-self.endQuery[lrr])
            if endMot < (self.startHsps[lrr+1]-3) :
                ## Si prochain motif au dela de la fin du profil +2
                ## etendre motif jusqu'à fin du profil
                self.endHsps[lrr]=self.endHsps[lrr]+(self.profile_len-self.endQuery[lrr])
            else:
                ## Si prochain LRR avant fin du profil -> pos=prochain motif -1
                self.endHsps[lrr]=self.startHsps[lrr+1]-1 
    
        self.endHsps[-1]=self.endHsps[-1]+(self.profile_len-self.endQuery[-1])
        ## Si out of range (fin motif après fin séquence)
        if (self.endHsps[-1]>len(self.seq)) :
            self.endHsps[-1]=len(self.seq)
            
            
            
    def inter_motif(self) :
        # Récupération des régions entre deux motifs LRR
        for i in range(1,len(self.startHsps)) :
            self.interMotifLen.append(self.startHsps[i]-self.endHsps[i-1]-1)
            self.interMotifPos.append(self.endHsps[i-1]+1)
            
            
    def inter_motif2(self) :
        ##stocker interMotif après recherche de LRR dans interMotif
        ##si pas de lrr en pos zero --> intermotif de 0 jusqu'au premier LRR
        i=0
        ## Ajout interLRR au début
        self.interMotifLen.append(self.startHsps[0]-1)
        self.interMotifPos.append(1)
        if (len(self.startHsps)>1) :
            for i in range(1,len(self.startHsps)) :
                self.interMotifLen.append(self.startHsps[i]-self.endHsps[i-1]-1)
                self.interMotifPos.append(self.endHsps[i-1]+1)       
        ##Ajout interLRR à la fin
        self.interMotifLen.append(len(self.seq)-self.endHsps[i])
        self.interMotifPos.append(self.endHsps[i]+1)
        
        
        
    def exclude_outliers(self) :
        ## exclude first and/or last lrr motif if it is far from next/previous
        ## motif and has low score
        ## outScore = distance(i,j)/mean(score(i),score(j))
        ## liste moyenneScore*Dist entre deux LRR
        M=[]
        if len(self.score)>2 :
            for i in range(len(self.score)-1):
                meanscore=(self.score[i]+self.score[i+1])/2
                M.append(self.interMotifLen[i]*meanscore)
            #meanDist=sum((self.interMotif[i] for i in self.interMotif))/len(self.interMotif)
            #while(exclude) :
            meanM=sum(M)/len(M)
            #Pour premier elm
            #si metric[i] > x fois la moyenne --> exclusion
            if M[0]>=3*meanM:
                ##exclude first
                self.score.pop(0)
                self.eval.pop(0)
                self.startHsps.pop(0)
                self.endHsps.pop(0)
                self.startQuery.pop(0)
                self.endQuery.pop(0)
                self.interMotifLen.pop(0)
                self.interMotifPos.pop(0)
                self.typeMotif.pop(0)
                
            if M[-1]>3*meanM:
                ##exclude last
                self.score.pop()
                self.eval.pop()
                self.startHsps.pop()
                self.endHsps.pop()
                self.startQuery.pop()
                self.endQuery.pop()
                self.interMotifLen.pop()
                self.interMotifPos.pop()
                self.typeMotif.pop()
                
                
                
#    def find_LRR_in_interLRR(self):
#        "Trouver des motifs non détectés par HMM dans les régions interLRR par regexp"
#        LRRexp="[LAIVFM]..[LAIVF]..[LAIVF].[LAIVF]..[NCST]"
#        for i in range(len(self.interMotifLen)):
#            if self.interMotifLen[i]>=15 :
#                start=self.interMotifPos[i]
#                #print(start)
#                sequence=self.seq[start:(self.interMotifPos[i]+self.interMotifLen[i]+1)]
#                #print(sequence)
#                ma=re.finditer(LRRexp,sequence)
#                for elm in ma :
#                    pos=start+elm.start()
#                    size=elm.end()-elm.start()
#                    #print(pos)
#                    #print(size)
#                    i=0
#                    while self.startHsps[i] < pos : ##trouver l'indice pour insertion
#                        i=i+1
#                    ##insertion dans les listes
#                    #print(i)
#                    self.score.insert(i,0)
#                    self.startHsps.insert(i,pos)
#                    self.endHsps.insert(i,pos+size)
#                    self.startQuery.insert(i,0)
#                    self.endQuery.insert(i,12)
#                    self.typeMotif.insert(i,"regex")
                
                
       
    def save_hit_to_file(self,filename):
        file=open(filename,"a")
        c=0
        ##Ecrire les LRR et InterLRR dans l'ordre des positions dans la prot avec 
        ## indice d'apparition
        for i in range(len(self.startHsps)-1) : ## -1 pour range interLRR (1 de moins que LRR)
            c=c+1
            file.write("%s;%s;%i;%i;%i;%i;%f;%s\n" % (self.id,self.typeMotif[i],c,self.startHsps[i],self.endHsps[i],self.endHsps[i]-self.startHsps[i]+1,self.eval[i],self.seq[self.startHsps[i]-1:self.endHsps[i]]))
            if (self.interMotifLen[i]) > 0 :
                c=c+1
                file.write("%s;%s;%i;%i;%i;%i;%f;%s\n" % (self.id,"interLRR",c,self.interMotifPos[i],self.interMotifPos[i]+self.interMotifLen[i]-1,self.interMotifLen[i],0,self.seq[self.interMotifPos[i]-1:self.interMotifPos[i]+self.interMotifLen[i]-1]))                
        # extraction dernier LRR
        file.write("%s;%s;%i;%i;%i;%i;%f;%s\n" % (self.id,self.typeMotif[-1],c+1,self.startHsps[-1],self.endHsps[-1],self.endHsps[-1]-self.startHsps[-1]+1,self.eval[-1],self.seq[self.startHsps[-1]-1:self.endHsps[-1]]))
        file.close()       
        
        
    def save_hit_interLRR(self,filename,protId,startPos):
        file=open(filename,"a")
        c=0
        for i in range(len(self.startHsps)) :
            c=c+1
            file.write("%s;%s;%i;%i;%i;%i;%f;%s\n" % (protId,self.typeMotif[i],c,startPos+self.startHsps[i]-1,startPos+self.endHsps[i]-1,self.endHsps[i]-self.startHsps[i]+1,self.eval[i],self.seq[self.startHsps[i]-1:self.endHsps[i]]))
            
        file.close()
        
        
    def save_hit_blast(self,filename,protId):
        file=open(filename,"a")
        c=0
        for i in range(len(self.startHsps)) :
            c=c+1
            file.write("%s;%s;%i;%i;%i;%i;%f;%s\n" % (protId,self.typeMotif[i],c,self.startHsps[i],self.endHsps[i],self.endHsps[i]-self.startHsps[i]+1,self.eval[i],self.seq[self.startHsps[i]-1:self.endHsps[i]]))
            
        file.close()
        
        
    def processed_data(self) :
        ##Execution de toutes les fonctions pour etendre les LRR et extraire les ilots
        self.inter_motif()
        self.exclude_outliers()
        #self.find_LRR_in_interLRR()
        self.extend_nter()
        self.extend_cter()
        self.interMotifLen=[]
        self.interMotifPos=[]
        self.inter_motif()
        
        
    def processed_interLRR(self) :
        ##execution des methodes pour etendre les LRR des interLRR
        self.extend_nter()
        self.extend_cter()
        #self.inter_motif2()

    def processed_blastRes(self) :
        ##Peut-importe les chevauchements, ils seront traites par le script Concat.
        for i in range(len(self.startHsps)) :
            ## extend n-ter sur la base du match 
            if(self.startHsps[i]-(self.startQuery[i]-1)<self.limitInf[i]) :
                self.startHsps[i]=self.limitInf[i]
            else :
                self.startHsps[i]=self.startHsps[i]-(self.startQuery[i]-1)

            ## extend c-ter sur la base du match 
            if(self.endHsps[i]+(self.profile_len-self.endQuery[i])>self.limitSup[i]):
                self.endHsps[i]=self.limitSup[i]
            else:
                self.endHsps[i]=self.endHsps[i]+(self.profile_len-self.endQuery[i])
    
    
#_________________
#
# end class hmmHit
#_________________   
