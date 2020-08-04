from imgrasp.graph.core import GenericGraph


class SIFParser(object):
    @staticmethod
    def parse(node_table, edge_table, id_columns=('id', 'id'), interaction_columns=('source','target')):
        pass

    def __generate_nodes(self, node_table, edge_table, id_columns=('id', 'id'), interaction_columns=('source','target'), edge_weight_col='weight'):
        srccol, tarcol = interaction_columns
        if node_table is None:
            node_names = set(edge_table[srccol].unique()) | set(edge_table[tarcol].unique())

