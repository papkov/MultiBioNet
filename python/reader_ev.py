from math import floor
import igraph
import pandas as pd


def read_file(filename):
    df = pd.read_csv(filename, sep='\t',
                     usecols=["ensg1", "ensg2", "score", "interaction"])
    return df


def progress_bar(current_num, total_num):
    if current_num % floor(total_num / 100) == 0:
        progress = current_num // round(total_num / 100)
        pb = '[' + '*' * round(progress) + '-' * (100 - progress) + '] ' + str(progress) + '%'
        print(pb, end='\r')


def set_graph(df):
    # Make lists of all parameters
    all_vertices = []
    all_edges = []
    all_interactions = []
    all_scores = []
    for index, row in df.iterrows():
        all_vertices.append(row['ensg1'])
        all_vertices.append(row['ensg2'])
        all_edges.append((row['ensg1'], row['ensg2']))
        all_interactions.append(row['interaction'])
        all_scores.append(row['score'])

    all_vertices = list(set(all_vertices))  # List of unique vertices

    graph = igraph.Graph()  # Create graph
    graph.vs["name"] = []
    graph.add_vertices(all_vertices)  # Add all vertices
    graph.add_edges(all_edges)  # Add all edges
    graph.es["interaction"] = all_interactions  # Add attribute "Interaction"
    graph.es["score"] = all_scores  # Add attribute "Score"

    igraph.summary(graph)
    return graph


if __name__ == '__main__':
    #df = read_file("../data/alz_intact_int_PPI.txt")
    #df = read_file("../data/intact_int.txt")
    #df = read_file("../data/alzcoexp_int_00001_10102016.txt")
    #df = read_file("../data/parkinson_intact_int_PPI.txt")
    df = read_file("../data/synapse_intact_int.txt")

    g = set_graph(df)
    igraph.plot(g)
