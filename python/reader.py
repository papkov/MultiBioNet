import igraph
import pandas as pd


def read_file(filename):
    df = pd.read_csv(filename, sep='\t',
                     usecols=["ensg1", "ensg2", "score", "interaction"])
    return df


def set_graph(df, verbal=False):
    graph = igraph.Graph()
    graph.vs["name"] = []
    graph.es["interaction"] = []
    graph.es["score"] = []

    for index, row in df.iterrows():
        ensg1 = row['ensg1']
        ensg2 = row['ensg2']
        score = row['score']
        interaction = row['interaction']

        if ensg1 not in graph.vs["name"]:
            graph.add_vertex(name=ensg1)
            ensg1_index = len(graph.vs["name"]) - 1
        else:
            ensg1_index = graph.vs["name"].index(ensg1)

        if ensg2 not in graph.vs["name"]:
            graph.add_vertex(name=ensg2)
            ensg2_index = len(graph.vs["name"]) - 1
        else:
            ensg2_index = graph.vs["name"].index(ensg2)

        if verbal:
            print("Add edge between %s, %s" % (ensg1_index, ensg2_index))

        graph.add_edges([(ensg1_index, ensg2_index)])
        graph.es["interaction"] = interaction
        graph.es["score"] = score

    igraph.summary(graph)
    return graph


df = read_file("../data/alz_intact_int_PPI.txt")
g = set_graph(df)
print(g.vs["name"])
igraph.plot(g)