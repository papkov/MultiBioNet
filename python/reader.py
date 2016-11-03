import igraph
import pandas as pd


def read_file(filename):
    df = pd.read_csv(filename, sep='\t',
                     usecols=["ensg1", "ensg2", "score", "interaction"])
    return df


def set_graph(df, verbal=False):
    graph = igraph.Graph()
    graph.vs["name"] = []

    for index, row in df.iterrows():

        new_edge = []

        for ensg in [row['ensg1'], row['ensg2']]:
            if ensg not in graph.vs["name"]:
                graph.add_vertex(name=ensg)
                new_edge.append(len(graph.vs["name"]) - 1)
            else:
                new_edge.append(graph.vs["name"].index(ensg))

        graph.add_edge(new_edge[0], new_edge[1],
                       interaction=row['interaction'], score=row['score'])
        if verbal:
            print("Add edge between %s, %s" % (new_edge[0], new_edge[1]))

    igraph.summary(graph)
    return graph


if __name__ == '__main__':
    df = read_file("../data/alz_intact_int_PPI.txt")
    g = set_graph(df)
    #print(g.vs["name"])
    igraph.plot(g)
