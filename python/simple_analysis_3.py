#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 10:24:10 2016

@author: kenelm
"""

# Simple analysis of datasets

from math import floor
import numpy as np
import igraph
import pandas as pd
import matplotlib.pyplot as plt
import random

from reader_ev2 import set_graph
from reader_ev2 import read_file


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


def plot_heatmap(all_vectors, ax):
    plt.imshow(all_vectors, aspect='auto', interpolation="nearest")
    plt.colorbar(orientation='vertical')
    ax.set_xticks( np.arange(0, 15, 1) )                                                       
    ax.set_xticklabels(['Int Alz PPI(d)',
                        'Int Alz PPI(j)',
                        'Int all PPI(d)',
                        'Int all PPI(j)',
                        'CoExp(d)',
                        'CoExp(j)',
                        'Int syn(d)',
                        'Int syn(j)',
                        'Int Par PPI(d)',
                        'Int Par PPI(j)',
                        'Dif Co(d)',
                        'Dif Co(j)',
                        'Alz ptw 1',
                        'Alz ptw 2',
                        'Alz DISEASES 1',
                        'Alz DISEASES 2'])
    plt.xticks(rotation='vertical')
    pos1 = ax.get_position() # get the original position 
    pos2 = [pos1.x0, pos1.y0 + 0.06,  pos1.width, pos1.height] 
    ax.set_position(pos2) # set a new position    

    
def pair_analysis(gene1, gene2, set_of_graphs, set_of_lists):
    analysis_vector = []
    for current_graph in set_of_graphs:
        if ( gene1 in current_graph.vs["name"] ) and ( gene2 in current_graph.vs["name"] ):
            sp = current_graph.shortest_paths(gene1, gene2, weights=None, mode='ALL')[0][0]
            inv_dist = 1/max(sp, 1)
            jaccard = current_graph.similarity_jaccard(vertices=[gene1,gene2], pairs=None, mode='ALL', loops=True)[0][1]
        else:
            inv_dist = 0            
            jaccard = 0
            
        analysis_vector.append(inv_dist)
        analysis_vector.append(jaccard)

    for current_list in set_of_lists:
        presence_code = [(gene1 in current_list) , (gene2 in current_list)]
        analysis_vector += presence_code
    
    return analysis_vector       


def save_features(all_vectors, graph, edges, file='features.csv'):
    with open(file, 'w') as f:
        output = ['%s\n' % ','.join([
                        'Gene1',
                        'Gene2',
                        'Int Alz PPI(d)',
                        'Int Alz PPI(j)',
                        'Int all PPI(d)',
                        'Int all PPI(j)',
                        'CoExp(d)',
                        'CoExp(j)',
                        'Int syn(d)',
                        'Int syn(j)',
                        'Int Par PPI(d)',
                        'Int Par PPI(j)',
                        'Dif Co(d)',
                        'Dif Co(j)',
                        'Alz ptw 1',
                        'Alz ptw 2',
                        'Alz DISEASES 1',
                        'Alz DISEASES 2'])]
        for i, edge in enumerate(edges):
            output.append(
                '{},{},{}\n'.format(graph.vs[edge[0]]["name"],
                                    graph.vs[edge[1]]["name"],
                                    ','.join([str(feature) for feature in all_vectors[i]])))
        f.writelines(output)


if __name__ == '__main__':
        
    # TODO: + Biogrid
    # TODO: + Diff co-expression
    # TODO: + STRING   
    
    df = read_file('../data/alz_intact_int_PPI.txt')
    gr_int_alz_PPI = set_graph(df)
    print('0. alz_intact_int_PPI.txt ' + 'is loaded...')
    
    df = read_file("../data/intact_int.txt")
    gr_int_all_PPI = set_graph(df)
    print('1. intact_int.txt ' + 'is loaded...')
    
    df = read_file("../data/alzcoexp_int_00001_10102016.txt")
    gr_int_coexpr = set_graph(df)    
    print('2. alzcoexp_int_00001_10102016.txt ' + 'is loaded...')
    
    df = read_file("../data/synapse_intact_int.txt")
    gr_int_syn_PPI = set_graph(df)
    print('3. synapse_intact_int.txt ' + 'is loaded...')
    
    df = read_file("../data/parkinson_intact_int_PPI.txt")
    gr_int_par_PPI = set_graph(df)
    print('4. parkinson_intact_int_PPI.txt ' + 'is loaded...')

    df = read_file('../data/dc_pairs_all_v2.txt')
    gr_dc_pairs_all = set_graph(df)
    print('5. dc_pairs_all_v2.txt ' + 'is loaded...')
    
    df = read_file("../data/alzheimer_pathway_genes.txt")
    ls_alz_ptw = get_genes_list(df, ['ensg'])
    print('6. alzheimer_pathway_genes.txt ' + 'is loaded...')

    df = read_file("../data/alz_DISEASES_ENSG")
    ls_alz_disease = get_genes_list(df, ['ensg'])
    print('7. alz_DISEASES_ENSG ' + 'is loaded...', len(ls_alz_disease))
    
    set_of_graphs = [gr_int_alz_PPI, gr_int_all_PPI, gr_int_coexpr, gr_int_syn_PPI, gr_int_par_PPI, gr_dc_pairs_all]
    set_of_lists = [ls_alz_ptw, ls_alz_disease]
        
    fig = plt.figure(1)
    
    # Analysis of known PPIs in a graph of Alzheimer PPIs
    all_alz_edges = gr_int_alz_PPI.get_edgelist()
    all_vectors = [pair_analysis(gr_int_alz_PPI.vs[edge[0]]["name"], gr_int_alz_PPI.vs[edge[1]]["name"], set_of_graphs, set_of_lists) for edge in all_alz_edges]

    ax = fig.add_subplot(1,3,1)
    plot_heatmap(all_vectors, ax)
    save_features(all_vectors, gr_int_alz_PPI, all_alz_edges, file='../data/alz1_features.csv')

    # Analysis of unknown PPIs in a graph of Alzheimer PPIs
    all_vertices = gr_int_alz_PPI.vs["name"]
    all_quazi_alz_edges = [random.sample(all_vertices, 2) for k in range(len(all_alz_edges))]
    all_vectors = [pair_analysis(couple[0], couple[1], set_of_graphs, set_of_lists) for couple in all_quazi_alz_edges]

    ax = fig.add_subplot(1,3,2)
    plot_heatmap(all_vectors, ax)
    save_features(all_vectors, gr_int_alz_PPI,
                  [(gr_int_alz_PPI.vs.find(couple[0]).index, gr_int_alz_PPI.vs.find(couple[1]).index) for couple in all_quazi_alz_edges],
                  file='../data/alz2_quazi_features.csv')
    
    # Analysis of unknown PPIs in a graph of all PPIs
    all_vertices = gr_int_all_PPI.vs["name"]
    all_quazi_alz_edges = [random.sample(all_vertices, 2) for k in range(len(all_alz_edges))]
    all_vectors = [pair_analysis(couple[0], couple[1], set_of_graphs, set_of_lists) for couple in all_quazi_alz_edges]
    
    ax = fig.add_subplot(1, 3, 3)
    plot_heatmap(all_vectors, ax)
    save_features(all_vectors, gr_int_all_PPI,
                  [(gr_int_all_PPI.vs.find(couple[0]).index, gr_int_all_PPI.vs.find(couple[1]).index) for couple in all_quazi_alz_edges],
                  file='../data/alz3_unknown_features.csv')


    plt.show()
