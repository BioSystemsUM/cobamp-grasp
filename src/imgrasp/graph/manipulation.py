from imgrasp.graph.core import GenericGraph
from typing import Iterable
from copy import deepcopy
from imgrasp.utilities.dict import dict_intersect, dict_union
from functools import reduce

def __manipulate(graphs: Iterable[GenericGraph], dict_func):
    args = [reduce(dict_func,(getattr(g, attr) for g in graphs)).values() for attr in ['nodes', 'edges']]
    return GenericGraph(*args)


def graph_intersection(graphs: Iterable[GenericGraph]): return __manipulate(graphs, dict_intersect)
def graph_union(graphs: Iterable[GenericGraph]): return __manipulate(graphs, dict_union)
