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
np.set_printoptions(precision=3)

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
               ['intact_int_filtered.txt'        , 'gr_int_all_PPI' , 1],
               ['alzcoexp_int_00001_10102016.txt', 'gr_int_coexpr'  , 1],
               ['synapse_intact_int.txt'         , 'gr_int_syn_PPI' , 1],
               ['parkinson_intact_int_PPI.txt'   , 'gr_int_par_PPI' , 1],
               ['dc_pairs_all_v2.txt'            , 'gr_dc_pairs_all', 1],
               ['alzheimer_pathway_genes.txt'    , 'ls_alz_ptw'     , 0]]    
        
    use_source_in_class = [0, 0, 1, 1, 1, 1, 1] # Flags - whether to use data from this source as features in classifier
       
    FNR_oneClass = 0.1 # False negative rate for one-class SVM     
    
    remove_repeats = False # Flag -- remove the duplicate entries in positives or negatives lists is required
    
    build_ROC = False # If True, builds ROC curve via LOO cross-validation. If False, just builds classifier
    
    two_class_SVM_wt = np.array([3]) # Skew between two classes
    # two_class_SVM_wt = 2**np.arange( -7, 8, 1, dtype = np.float) 
    
    max_distance_predict = 5 # Maximum distance in neighbours_of_gwas.pickle
    
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
    set_of_random_edges = [random.sample(all_vertices, 2) for _ in range(10**4)]
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
    
    print ('Number of positives: ' + str(len(positives)))
    print ('Number of negatives: ' + str(len(negatives)))

    # ROC 
    loo = LeaveOneOut()
    loo.get_n_splits(X)
    
    FPR = []
    TPR = []
    ERR = []

    if build_ROC:        
        for wt in two_class_SVM_wt:
        
            classifier = svm.SVC(kernel='rbf', class_weight={1: wt}, probability=False)
            FP_num = 0
            TP_num = 0
            ER_num = 0
            for train, test in loo.split(X):
                prediction = classifier.fit(X[train], y[train]).predict(X[test]) 
                FP_num += (prediction == 1)*(1 - y[test])
                TP_num += (prediction == 1)*y[test]            
                ER_num += not(prediction == y[test])
        
            TPR.append(TP_num[0]/len(positives))
            FPR.append(FP_num[0]/len(negatives))
            ERR.append(ER_num/len(X))
        
            print((wt,FP_num[0],TP_num[0],ER_num))
    
        plt.plot(FPR, TPR)
        plt.xlabel('False postitive rate')
        plt.ylabel('True positive rate')    
        plt.plot([0,1],[0,1])
        plt.xlim([0,1])
        plt.ylim([0,1])
        plt.grid()
    else:
        classifier = svm.SVC(kernel='rbf', class_weight={1: two_class_SVM_wt}, probability=False)
        classifier.fit(X, y)
        with open('neighbours_of_gwas.pickle', 'rb') as f:
            neighbours_of_gwas = pickle.load(f)
            
        counter_all = 0
        counter_pr = 0
        x = 0
        for item in neighbours_of_gwas:
            x += 1
            print( (x,len(neighbours_of_gwas)) )
            
            neighbours_list = neighbours_of_gwas[item]
            for nb_tuple in neighbours_list:
                if ( nb_tuple[1] <= max_distance_predict ):
                    counter_all += 1
                    features = np.array( [pair_analysis(item, nb_tuple[0], set_of_graphs_class.copy(), set_of_lists_class.copy(), opt_del_edge)] )
                    if classifier.predict(features):
                        counter_pr += 1
#                        print(item+' '+nb_tuple[0])
            
        print((counter_all, counter_pr))           
    
#(12.0, 0, 217)
#(12.1, 0, 217)
#(12.199999999999999, 0, 217)
#(12.299999999999999, 0, 217)
#(12.399999999999999, 0, 217)
#(12.499999999999998, 0, 217)
#(12.599999999999998, 0, 217)
#(12.699999999999998, 0, 217)
#(12.799999999999997, 0, 219)
#(12.899999999999997, 0, 221)
#(12.999999999999996, 0, 221)
#(13.099999999999996, 0, 221)
#(13.199999999999996, 0, 221)
#(13.299999999999995, 0, 221)
#(13.399999999999995, 0, 221)
#(13.499999999999995, 0, 221)
#(13.599999999999994, 0, 221)
#(13.699999999999994, 0, 221)
#(13.799999999999994, 0, 221)
#(13.899999999999993, 0, 221)
#(13.999999999999993, 0, 221)
#(14.099999999999993, 0, 221)
#(14.199999999999992, 0, 221)
#(14.299999999999992, 0, 221)
#(14.399999999999991, 1, 221)
#(14.499999999999991, 1, 221)
#(14.599999999999991, 1, 222)
#(14.69999999999999, 2, 222)
#(14.79999999999999, 2, 229)
#(14.89999999999999, 2, 229)
#(14.999999999999989, 2, 229)
#(15.099999999999989, 298, 229)
#(15.199999999999989, 298, 229)
#(15.299999999999988, 299, 229)
#(15.399999999999988, 300, 229)
#(15.499999999999988, 300, 251)
#(15.599999999999987, 300, 251)
#(15.699999999999987, 565, 251)
#(15.799999999999986, 565, 251)
#(15.899999999999986, 565, 251)
#(15.999999999999986, 681, 251)
#(16.099999999999987, 681, 262)
#(16.199999999999985, 709, 262)
#(16.299999999999983, 720, 262)
#(16.399999999999984, 722, 262)
#(16.499999999999986, 722, 262)
#(16.599999999999984, 722, 262)