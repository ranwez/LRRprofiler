# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 13:45:01 2019

@author: gottince
"""

#=================================================
#  Project : LRR proteins annotation
#
#hmmsearch results exctraction, compilation and treatment
#=================================================

    
    # filtering results on protein e-value
def hit_eval_filter(res, hit_ths) :
    filt = lambda hit: hit.evalue<=hit_ths
    res=res.filter(filt)   
    
def hsp_eval_filter(hit_res, eval_ths) :
    filt = lambda hsp: hsp.evalue<eval_ths
    hit_res=hit_res.filter(filt)
    return hit_res
    
    # filtering hsp on motif start position
def hsp_start_filter(hit_res, hsp_ths) :
    filt = lambda hsp: hsp.query_start<hsp_ths
    hit_res=hit_res.filter(filt)
    return hit_res
    
    # filtering hsp on motif end position
def hsp_end_filter(hit_res, hsp_ths) :
    filt = lambda hsp: hsp.query_end>hsp_ths
    hit_res=hit_res.filter(filt)
    return hit_res
    
##
    # extract hit proteins
def extract_hit(proteome, results) :
    for hit in results :
        proteome.add_protein(hit.id,hit.seq_len)
    
    
    
    