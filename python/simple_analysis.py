#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 10:24:10 2016

@author: kenelm
"""

# Simple analysis of datasets

from math import floor
import igraph
import pandas as pd

from reader_ev import set_graph


def read_file(filename):
    df = pd.read_csv(filename, sep='\t')
    return df
    
def progress_bar(current_num, total_num):
    if ( current_num%floor(total_num/100) == 0 ):
        progress = current_num//round(total_num/100)
        pb = '['+'*'*round(progress)+'-'*(100-progress)+'] ' + str(progress) + '%'
        print(pb, end = '\r')

def get_genes_list(df, set_of_ensg):

    # Make lists of genes
    all_vertices = []
    for ensg in set_of_ensg:
        all_vertices += [row[ensg] for index, row in df.iterrows()]
        
    all_vertices = set(all_vertices) # Set of unique genes

    return all_vertices

if __name__ == '__main__':
    df = read_file("../data/alz_intact_int_PPI.txt")
    genes_in_alz_PPI = get_genes_list(df,['ensg1', 'ensg2'])
    print( 'Genes in Alzheimer PPIs:', len(genes_in_alz_PPI) )
    
    df = read_file("../data/intact_int.txt")
    genes_in_all_PPI = get_genes_list(df,['ensg1', 'ensg2'])
    print( 'Genes in all PPIs: ', len(genes_in_all_PPI) )
    print( 'Common genes in Alzheimer PPI and all PPI: ', len( genes_in_alz_PPI.intersection(genes_in_all_PPI) ) )
    
    df = read_file("../data/alzcoexp_int_00001_10102016.txt")
    genes_in_all_CoExp = get_genes_list(df, ['ensg1', 'ensg2'])
    print( 'Genes in co-expression network: ', len(genes_in_all_CoExp) )
    print( 'Common genes in Alzheimer PPI and co-expression network: ', len( genes_in_alz_PPI.intersection(genes_in_all_CoExp) ) )    

    df = read_file("../data/alzheimer_pathway_genes.txt")
    genes_in_alz_ptw = get_genes_list(df, ['ensg'])
    print( 'Genes in Alzheimer pathways: ', len(genes_in_alz_ptw) )  
    print( 'Common genes in Alzheimer PPI and Alzheimer pathways: ', len( genes_in_alz_PPI.intersection(genes_in_alz_ptw) ) )
    
    df = read_file("../data/GWAS_genes.txt")
    genes_in_alz_gws = get_genes_list(df, ['ensg'])
    print( 'Genes in GWAS: ', len(genes_in_alz_gws) )  
    print( 'Common genes in Alzheimer PPI and GWAS: ', len( genes_in_alz_PPI.intersection(genes_in_alz_gws) ) )

    df = read_file("../data/synapse_intact_int.txt")
    genes_in_syn_PPI = get_genes_list(df, ['ensg1', 'ensg2'])
    print( 'Genes in synaps PPI: ', len(genes_in_syn_PPI) )  
    print( 'Common genes in Alzheimer PPI and synaps PPI: ', len( genes_in_alz_PPI.intersection(genes_in_syn_PPI) ) )
    
    df = read_file("../data/parkinson_intact_int_PPI.txt")
    genes_in_par_PPI = get_genes_list(df, ['ensg1', 'ensg2'])
    print( 'Genes in Parkinson PPI: ', len(genes_in_par_PPI) )  
    print( 'Common genes in Alzheimer PPI and Parkinson PPI: ', len( genes_in_alz_PPI.intersection(genes_in_par_PPI) ) )
    