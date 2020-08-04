import unittest

from imgrasp.graph.core import Node, WeightedEdge, GenericGraph
from imgrasp.graph.matrix import adjacency_matrix
from imgrasp.graph.manipulation import graph_intersection, graph_union


class MyTestCase(unittest.TestCase):
    def test_something(self):
        node_dict = {s: Node(s, {}) for s in ['I1', 'I2', 'N1', 'N2', 'M1', 'M2']}
        edge_dict = {x + '->' + y: WeightedEdge(node_dict[x], node_dict[y], 0, x + '->' + y, {}) for x, y in
                     [['I1', 'N1'], ['I2', 'N2'], ['N1', 'M1'], ['N1', 'M2'], ['N2', 'M1'], ['N2', 'M2']]}
        weights = [2, 3, 5, 1, 8, 10]
        for k,v in zip(edge_dict.items(), weights): k[1].weight = v

        graph = GenericGraph(node_dict.values(), edge_dict.values())

        print(adjacency_matrix(graph))
        print(adjacency_matrix(graph, directed=True))
        print(adjacency_matrix(graph, weighted=True))
        print(adjacency_matrix(graph, directed=True, weighted=True))

        other_node_ids = ['I1','N2','N1']
        other_nodes = [node_dict[k] for k in other_node_ids]
        other_edges = [e for k,e in edge_dict.items() if
                       len(set(map(lambda x: x.identifier, e.involved_nodes)) & set(other_node_ids)) > 1]
        another_graph = GenericGraph(other_nodes, other_edges)
        intersect = graph_intersection([graph, another_graph])
        intersect.remove_unconnected_edges()

if __name__ == '__main__':
    unittest.main()