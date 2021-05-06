import numpy as np
import pandas as pd
from typing import Union
from imgrasp.graph.core import WeightedEdge, Node, GenericGraph

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

def graph_from_dataframe(frame: Union[pd.DataFrame, np.ndarray], **kwargs):
    if isinstance(frame, np.ndarray):
        node_names, edge_names = [kwargs[k] if k in kwargs.keys() else list(map(str,range(frame.shape[0])))
                                  for k in ['nodes', 'edges']]
    elif isinstance(frame, pd.DataFrame):
        node_names, edge_names = [k.tolist() for k in [frame.index, frame.columns]]
        frame = frame.values

    else:
        raise Exception('Invalid frame argument')

    edge_dict = [[] for k in range(len(edge_names))]

    node_nodes = [Node(identifier=k, annotation={}) for k in node_names]
    edge_nodes = [Node(identifier=k, annotation={}) for k in edge_names]

    edges_to_add = []
    for sense, mat in zip([True, False],[frame < 0, frame > 0]):
        vals = frame[mat]
        rowind, colind = np.where(mat)
        for r,c,v in zip(rowind, colind, vals):
            e_name = '_'.join([edge_names[c],str(len(edge_dict[c]))])
            n_to_add = node_nodes[r], edge_nodes[c]
            if sense:
                n_to_add = n_to_add[::-1]
            e_to_add = WeightedEdge(source=n_to_add[0], incident=n_to_add[1],
                                    weight=v, identifier=e_name, annotation={})
            edge_dict[c].append(e_to_add.identifier)
            edges_to_add.append(e_to_add)
    print('x')

    return GenericGraph(node_nodes+edge_nodes, edges_to_add), dict(zip(edge_names,edge_dict))

