import numpy as np

def adjacency_matrix(graph, directed=False, weighted=False):
    ord_nodes = list(graph.nodes.keys())
    ord_dict = {v:k for k,v in enumerate(ord_nodes)}
    amat = np.zeros([len(ord_nodes)]*2, dtype=int)
    for i,e in enumerate(graph.edges):
        edge = graph.edges[e]
        src, inc = (ord_dict[k.identifier] for k in (edge.source, edge.incident))
        amat[src, inc] += 1 if not weighted else edge.weight
    if not directed:
        nzi = amat.nonzero()
        amat[nzi[::-1]] = amat[nzi]
    return amat

