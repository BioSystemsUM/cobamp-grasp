from indexed import IndexedOrderedDict

class IdentifiableEntity(object):
    def __init__(self, identifier):
        self.identifier = str(identifier)

    @property
    def identifier(self):
        return self._id

    @identifier.setter
    def identifier(self, value):
        self._id = value

    def __repr__(self):
        return self.__class__.__name__+' '+self.identifier+' at '+str(hex(id(self)))

def do_something_if_valid(obj: IdentifiableEntity, func_to_do, instance, where, should_be_there=False):
    if isinstance(obj, instance):
        proceed = (should_be_there == (obj.identifier in where.keys()))
        if proceed:
            func_to_do(obj)
        else:
            raise Exception(obj.identifier + (' does not exist.' if should_be_there else ' already exists.'))
    else:
        raise TypeError(str(type(obj)), ' is not a valid instance of '+ str(instance) +'.')

def add_identifiable_to_dict(x,d): d[x.identifier] = x
def rm_identifiable_from_dict(x,d): del d[x.identifier]

class AnnotateableEntity(object):
    def __init__(self, annotation):
        self.annotation = annotation
    @property
    def annotation(self):
        return self._annotation

    @annotation.setter
    def annotation(self, value):
        self._annotation = value

class GenericGraph(object):
    def __init__(self, nodes, edges):
        self._nodes = IndexedOrderedDict()
        self._edges = IndexedOrderedDict()
        self.add(nodes, edges)

    @property
    def nodes(self):
        return self._nodes

    @property
    def edges(self):
        return self._edges

    def add(self, nodes: IndexedOrderedDict = None, edges: IndexedOrderedDict = None):
        if nodes is None: nodes = ()
        if edges is None: edges = ()

        for l, d, typ in zip([nodes, edges],[self._nodes, self._edges],[Node, Edge]):
            for item in l: do_something_if_valid(item, lambda x: add_identifiable_to_dict(x, d), typ, d, False)


    def remove(self, nodes = None, edges = None, remove_unconnected_edges=True):
        if nodes is None: nodes = ()
        if edges is None: edges = ()

        for l, d, typ in zip([nodes, edges],[self._nodes, self._edges],[Node, Edge]):
            for item in l: do_something_if_valid(item, lambda x: add_identifiable_to_dict(x, d), typ, d, True)

        if remove_unconnected_edges: self.remove_unconnected_edges()

    def remove_unconnected_edges(self):
        rmv = [e for k,e in self.edges.items() if len({e.source.identifier, e.incident.identifier} & set(self.nodes.keys())) < 2]
        if len(rmv) > 0:
            self.remove(edges=rmv, remove_unconnected_edges=True)

class Node(IdentifiableEntity, AnnotateableEntity):
    def __init__(self, identifier, annotation):
        IdentifiableEntity.__init__(self, identifier)
        AnnotateableEntity.__init__(self, annotation)


class Edge(IdentifiableEntity, AnnotateableEntity):
    def __init__(self, source, incident, identifier, annotation):
        IdentifiableEntity.__init__(self, identifier)
        AnnotateableEntity.__init__(self, annotation)
        self.source = source
        self.incident = incident


    @property
    def source(self) -> Node:
        return self._source
    
    @source.setter
    def source(self, value):
        self._source = value
    
    @property
    def incident(self) -> Node:
        return self._incident
    
    @incident.setter
    def incident(self, value):
        self._incident = value

    @property
    def involved_nodes(self):
        return (self.source, self._incident)

    def __repr__(self):
        return super(IdentifiableEntity, self).__repr__() + ' connecting ' + \
               self._source.identifier + ' to ' + self._incident.identifier

class WeightedEdge(Edge):
    def __init__(self, source, incident, weight, identifier, annotation):
        super().__init__(source, incident, identifier=identifier, annotation=annotation)
        self.weight = weight

    @property
    def weight(self):
        return self._weight

    @weight.setter
    def weight(self, value):
        self._weight = value
