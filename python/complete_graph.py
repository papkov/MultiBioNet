import igraph
import pandas as pd


def read_interactions(filename):
    df = pd.read_csv(filename, sep='\t',
                     usecols=["ensg1", "ensg2", "score", "interaction"])
    return df


def set_of_ensg(file_path):
    s = set()
    with open(file_path, "r") as f:
        for gene in f.readlines()[1:]:
            s.add(gene.strip())
    return s


def add_edges(graph, df, label):

    all_vertices = [row['ensg1'] for index, row in df.iterrows()]
    all_vertices += [row['ensg2'] for index, row in df.iterrows()]
    new_vertices = list(set(all_vertices) - set(graph.vs["name"]))
    graph.add_vertices(new_vertices)

    all_edges = [(row['ensg1'], row['ensg2']) for index, row in df.iterrows()]

    #all_interactions = [row['interaction'] for index, row in df.iterrows()]
    #all_scores = [row['score'] for index, row in df.iterrows()]
    #all_labels = [label for index, row in df.iterrows()]

    graph.add_edges(all_edges)
    #TODO
    #graph.es["interaction"] = all_interactions
    #graph.es["score"] = all_scores
    #graph.es["label_t"] = all_labels

    print("+ ", label)
    igraph.summary(graph)
    print()

    return graph


def add_vertices(graph, set_of_vertices, source):

    new_vertices = list(set_of_vertices - set(graph.vs["name"]))
    graph.add_vertices(new_vertices)

    #TODO
    #for ensg in set_of_vertices:
    #    g.vs.find(ensg)[source] = True

    print('+', source)
    igraph.summary(graph)
    print()

    return graph


def add_alz(graph):
    graph = add_edges(graph, read_interactions("../data/alz_intact_int_PPI.txt"),
                      "alz_intact")
    graph = add_vertices(graph, set_of_ensg("../data/GWAS_genes_2.txt"),
                         "GWAS")
    graph = add_vertices(graph, set_of_ensg("../data/alzheimer_pathway_genes_2.txt"),
                         "pathway")
    graph.vs["alz"] = [True for i in range(graph.vcount())]
    graph.vs["color"] = ["blue" for i in range(graph.vcount())]

    return graph


def add_non_alz(graph):
    graph = add_edges(graph, read_interactions("../data/alzcoexp_int_00001_10102016.txt"),
                      "alz_coexp")
    graph = add_edges(graph, read_interactions("../data/parkinson_intact_int_PPI.txt"),
                      "parkinson_int")
    graph = add_edges(graph, read_interactions("../data/synapse_intact_int.txt"),
                      "synapse_int")
    graph = add_edges(graph, read_interactions("../data/intact_int.txt"),
                      "all_intact")

    return graph


def delete_far_vertices(graph):
    dists = [float("Inf") for i in range(graph.vcount())]
    for v in graph.vs.select(alz_eq=True):
        cur_dists = graph.get_shortest_paths(v)
        cur_dists = [len(path) - 1 for path in cur_dists]
        for i in range(graph.vcount()):
            if cur_dists[i] > -1:
                dists[i] = min(dists[i], cur_dists[i])

    graph.vs["dist"] = dists
    graph.delete_vertices(graph.vs.select(dist_gt=2))

    print("- far vertices")
    igraph.summary(graph)
    print()

    return graph


if __name__ == '__main__':

    g = igraph.Graph()
    g.vs["name"] = []

    g = add_alz(g)
    g = add_non_alz(g)

    g = delete_far_vertices(g)

    g.simplify()
    print("- multiedges and loops")
    igraph.summary(g)

    #igraph.plot(g)
