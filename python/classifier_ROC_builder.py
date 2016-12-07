#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 10:24:10 2016

@author: kenelm
"""

# Simple analysis of datasets
#from IPython import get_ipython
#get_ipython().magic('reset -sf')

import sys
from math import floor
import numpy as np
import igraph
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import random
from sklearn import svm
from sklearn.model_selection import LeaveOneOut, StratifiedKFold
from sklearn.metrics import roc_curve, auc
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
        
    all_vertices = list( set(all_vertices) ) # Set of unique genes

    return all_vertices

def plot_heatmap(all_vectors, ax):
    plt.imshow(all_vectors,aspect = 'auto', interpolation="nearest")
    plt.colorbar(orientation='vertical')
    ax.set_xticks( np.arange(0, 15, 1) )                                                       
    ax.set_xticklabels(['Int Alz PPI(d)','Int Alz PPI(j)','Int all PPI(d)','Int all PPI(j)','CoExp(d)','CoExp(j)','Int syn(d)','Int syn(j)','Int Par PPI(d)','Int Par PPI(j)','Dif Co(d)','Dif Co(j)','Alz ptw 1','Alz ptw 2'])
    plt.xticks(rotation='vertical')
    pos1 = ax.get_position() # get the original position 
    pos2 = [pos1.x0, pos1.y0 + 0.06,  pos1.width, pos1.height] 
    ax.set_position(pos2) # set a new position    

    
def pair_analysis(gene1, gene2, set_of_graphs, set_of_lists, opt_del_edge):
    analysis_vector = []
    for index, current_graph in enumerate( set_of_graphs ): 
        if ( gene1 in current_graph.vs["name"] ) and ( gene2 in current_graph.vs["name"] ):
            gene1_ID = current_graph.vs.find(name=gene1)
            gene2_ID = current_graph.vs.find(name=gene2)                
            edge_id = current_graph.get_eid(gene1_ID.index, gene2_ID.index, directed=False, error=False)

            if opt_del_edge[index] and (edge_id >= 0):
                new_graph = current_graph.copy()
                new_graph.delete_edges(edge_id)
                sp = new_graph.shortest_paths(gene1, gene2, weights=None,mode='ALL')[0][0]
                inv_dist = 1/max(sp, 1)
                jaccard = new_graph.similarity_jaccard(vertices=[gene1,gene2], pairs=None, mode='ALL', loops=True)[0][1]
            else:
                sp = current_graph.shortest_paths(gene1, gene2, weights=None,mode='ALL')[0][0]
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
        
    
if __name__ == '__main__':
    
              
    create_new_picle = False # Get from txt files and make new pickle with graphs - True; just get from pickle - False
           
    # List of used sources: File name / Variable name / Graph - 1, List - 0 /
    sources = [
               ['alz_intact_int_PPI.txt'         , 'gr_int_alz_PPI' , 1],
               ['intact_int.txt'                 , 'gr_int_all_PPI' , 1],
               ['alzcoexp_int_00001_10102016.txt', 'gr_int_coexpr'  , 1],
               ['synapse_intact_int.txt'         , 'gr_int_syn_PPI' , 1],
               ['parkinson_intact_int_PPI.txt'   , 'gr_int_par_PPI' , 1],
               ['dc_pairs_all_v2.txt'            , 'gr_dc_pairs_all', 1],
               ['alzheimer_pathway_genes.txt'    , 'ls_alz_ptw'     , 0]]    
        
    use_source_in_class = [0, 1, 1, 1, 1, 1, 1] # Flags - whether to use data from this source as features in classifier
       
    FNR_oneClass = 0.1 # False negative rate for one-class SVM     
    
    remove_repeats = False # Flag -- remove the duplicate entries in positives or negatives lists is required
    
    # ----------- Long and boring re-loading of required data --------------------
    # Here I take a raw data from txt (csv) files, convert them into iGraph and
    # save these graphs into a "all_graphs.pickle" file. After that program breaks.
    
    if create_new_picle:
        set_of_all_graphs = [] # Set of all the loaded graphss
        set_of_all_lists  = [] # Set of all the loaded lists

        for index, source in enumerate(sources):        
             
            df = read_file('../data/' + source[0])
            if source[2]:
                command_load   = source[1] + "= set_graph(df)"
                command_append = "set_of_all_graphs.append (" + source[1] + ")"                    
            else:
                command_load = source[1] + "= get_genes_list(df, ['ensg'])"
                command_append = "set_of_all_lists.append (" + source[1] + ")"
                                                       
            exec(command_load)
            exec(command_append)
            print (str(index)+'. ' + source[0] + ' is loaded and converted.')
                                
        with open('all_graphs.pickle', 'wb') as f:
            pickle.dump([set_of_all_graphs, set_of_all_lists], f)   
            
        sys.exit("New pickle file with graphs and lists is created.")
    # ----------- End of long and boring re-loading of required data --------------------
 
    # --------- Get the graphs and lists from picle file --------------
    with open('all_graphs.pickle', 'rb') as f:
        set_of_all_graphs, set_of_all_lists = pickle.load(f)
    
    set_of_graphs_class = [] # Set of graphs which will be used in a classifier
    set_of_lists_class  = [] # Set of lists which will be used in a classifier

    all_graphs = '' # Temporary string for names of all variables storing graphs
    all_lists  = '' # Temporary string for names of all variables storing lists
    
    graphs_for_class = '' # Temporary string for names of variables storing graphs, used in classifier
    lists_for_class = ''  # Temporary string for names of variables storing lists, used in classifier
    
    for index, source in enumerate(sources):                    
                
        if source[2]:
            all_graphs += ( source[1] + ',' )
            if use_source_in_class[index]:
                graphs_for_class += ( source[1] + ',' )
                print (str(index)+'. ' + source[0] + ' will be used in classifier!')
            else:
                print (str(index)+'. ' + source[0] + ' WILL NOT be used in classifier!')
        else:
            all_lists += ( source[1] + ',' )
            if use_source_in_class[index]:
                lists_for_class += ( source[1] + ',' )
                print (str(index)+'. ' + source[0] + ' will be used in classifier!')
            else:
                print (str(index)+'. ' + source[0] + ' WILL NOT be used in classifier!')
               
    graphs_for_class = ''.join(list(graphs_for_class)[:-1]) # Remove extra comma
    lists_for_class  = ''.join(list(lists_for_class)[:-1])  # Remove extra comma
    all_graphs = ''.join(list(all_graphs)[:-1]) # Remove extra comma
    all_lists  = ''.join(list(all_lists)[:-1])  # Remove extra comma
                
    exec('['+all_graphs+'] = set_of_all_graphs') 
    exec('['+all_lists +'] = set_of_all_lists')                     

    exec('set_of_graphs_class = [' + graphs_for_class + ']' )        
    exec('set_of_lists_class  = [' + lists_for_class + ']' )
        
    # --------- End of getting the graphs and lists from picle file -----------------------
                     
        
    # Analysis of known PPIs in a graph of Alzheimer PPIs
    all_alz_edges = gr_int_alz_PPI.get_edgelist()
    opt_del_edge = [0, 0, 0, 0, 0]
    all_vectors1 = [pair_analysis(gr_int_alz_PPI.vs[edge[0]]["name"], gr_int_alz_PPI.vs[edge[1]]["name"], set_of_graphs_class.copy(), set_of_lists_class.copy(), opt_del_edge) for edge in all_alz_edges]
    if remove_repeats:
        all_vectors1 = [list(x) for x in set(tuple(x) for x in all_vectors1)] # Leave only unique entires
    positives    = [analysis_vector for analysis_vector in all_vectors1 if sum(analysis_vector)>0]
    clf = svm.OneClassSVM(nu=FNR_oneClass, kernel="linear")
    clf.fit(positives)    
    
    # Analysis of unknown PPIs in a graph of all PPIs
    all_vertices = gr_int_all_PPI.vs["name"]
    set_of_random_edges = [random.sample(all_vertices, 2) for _ in range(10**3)]
    set_of_random_edges = [pair for pair in set_of_random_edges if ( pair[0] not in gr_int_alz_PPI.vs["name"] ) and ( pair[1] not in gr_int_alz_PPI.vs["name"] )]

    opt_del_edge = [0, 0, 0, 0, 0]
    all_vectors2 = [pair_analysis(edge[0], edge[1], set_of_graphs_class.copy(), set_of_lists_class.copy(), opt_del_edge) for edge in set_of_random_edges]
    all_vectors2 = [analysis_vector for analysis_vector in all_vectors2 if sum(analysis_vector)>0]
    if remove_repeats:
        all_vectors2 = [list(x) for x in set(tuple(x) for x in all_vectors2)] # Leave only unique entires
    predictions  = clf.predict(all_vectors2)
    negatives    = [vector for index, vector in enumerate(all_vectors2) if predictions[index] == -1]

    X = np.array( positives + negatives )
    y = np.array([1]*len(positives) + [0]*len(negatives))
    
    # ROC 
    loo = LeaveOneOut()
    loo.get_n_splits(X)
    
    alpha = []
    beta  = []

    wt_list = [1/100, 11.5, 100] # Skew between two classes
    for wt in wt_list:
        
        classifier = svm.SVC(kernel='rbf', class_weight={1: wt}, probability=False)
        FNR = 0
        TNR = 0    
        for train, test in loo.split(X):
            if y[test]: 
                FNR += (classifier.fit(X[train], y[train]).predict(X[test]) == 0)
            else:
                TNR += (classifier.fit(X[train], y[train]).predict(X[test]) == 0)
        
        alpha.append(FNR[0]/len(positives))
        beta.append(TNR[0]/len(negatives))
        
        print((FNR/len(positives),TNR/len(negatives)))
    
    plt.plot(alpha,beta)
    plt.ylabel('True negative rate')
    plt.xlabel('False negative rate')
    plt.plot([0,1],[0,1])
    plt.xlim([0,1])
    plt.ylim([0,1])
    plt.grid()
    