from math import floor
import igraph
import pandas as pd


def read_file(filename):
    df = pd.read_csv(filename, sep='\t')
    return df

def progress_bar(current_num, total_num):
    if ( current_num%floor(total_num/100) == 0 ):
        progress = current_num//round(total_num/100)
        pb = '['+'*'*round(progress)+'-'*(100-progress)+'] ' + str(progress) + '%'
        print(pb, end = '\r')

def set_graph(df):

    # Make lists of all parameters
    all_vertices = [];
    all_edges = []
    all_interactions = []
    all_scores = []
    for index, row in df.iterrows():
        all_vertices.append( row['ensg1'] )
        all_vertices.append( row['ensg2'] )
        all_edges.append( ( row['ensg1'] , row['ensg2'] ) )
#        all_interactions.append(row['interaction'])
#        all_scores.append(row['score'])                  
        
    all_vertices = list(set(all_vertices)) # List of unique vertices
        
    graph = igraph.Graph() # Create graph
    graph.vs["name"] = []
    graph.add_vertices(all_vertices) # Add all vertices
    graph.add_edges(all_edges) # Add all edges
#    graph.es["interaction"] = all_interactions # Add attribute "Interaction"
#    graph.es["score"] = all_scores # Add attribute "Score"
        
    graph.simplify(multiple = True, loops = False) # Remove multiple edges, remain loops
    
    return graph

if __name__ == '__main__':
    df = read_file("../data/alz_intact_int_PPI.txt")
    #df = read_file("../data/intact_int.txt")
    g = set_graph(df)
  #  print( g.shortest_paths('ENSG00000142192', 'ENS00000166313', weights=None,mode='ALL') )
    igraph.plot(g)    