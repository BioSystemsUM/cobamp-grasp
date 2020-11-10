from functools import reduce
import pandas as pd
from imgrasp.graph.core import GenericGraph, WeightedEdge, Node, Edge


class SIFParser(object):
    @staticmethod
    def parse(edge_table:pd.DataFrame, interaction_columns=('source','target'), weight_column='weight'):
        srccol, tarcol = interaction_columns
        available_nodes = reduce(lambda x, y: x | y, (set(col.unique())
                                                      for coln, col in
                                                      edge_table[[srccol, tarcol]].iteritems()))

        graph_nodes = {g: Node(g, {}) for g in available_nodes}
        if weight_column in edge_table.columns:
            graph_edges = {str(i): WeightedEdge(source=graph_nodes[row[srccol]], incident=graph_nodes[row[tarcol]],
                                                weight=row[weight_column], identifier=str(i), annotation={})
                           for i, row in edge_table.iterrows()}
        else:
            graph_edges = {str(i): Edge(source=graph_nodes[row[srccol]], incident=graph_nodes[row[tarcol]],
                                        identifier=str(i), annotation={}) for i, row in edge_table.iterrows()}

        return GenericGraph(graph_nodes.values(), graph_edges.values())
