#!/usr/bin/env python3
# -*- coding: utf-8 -*-
def getGraphEdgeList(graph):
    return graph.edges

def getGraphEdgeLabelList(graph):
    graphEdgeList=[]
    for edge in graph.edges:
        vertex1, vertex2=edge
        graphEdgeList.append((graph.nodes[vertex1]['label'], graph.nodes[vertex2]['label']))
    return graphEdgeList

def formatEdgeTuple(edgeTuple):
    return '('+str(edgeTuple[0])+','+str(edgeTuple[1])+')'

def isGraphEdgeDifferenceEmpty(graph1,graph2):    
    if graph1.number_of_edges()!=graph2.number_of_edges():
        return False
    
    for e in graph1.edges():
        if not graph2.has_edge(*e):
            return True

    return False

def isGraphEdgeIntersectionEmpty(graph1,graph2):    
    if graph1.number_of_edges()<=graph2.number_of_edges():
        for e in graph1.edges():
            if graph2.has_edge(*e):
                return False
    else:
        for e in graph2.edges():
            if graph1.has_edge(*e):
                return False

    return True

def isGraphVertexIntersectionEmpty(graph1,graph2, omitList=[]):    
    if graph1.number_of_nodes()<=graph2.number_of_nodes():
        for vertex in graph1.nodes():
            if (not vertex in omitList) and graph2.has_node(vertex):
                return False
    else:
        for vertex in graph2.nodes():
            if (not vertex in omitList) and graph1.has_node(vertex):
                return False

    return True

def graphUnion(graph1, graph2):
    union = graph1.fresh_copy()
    # add graph attributes, graph2 attributes take precedent over graph1 attributes
    union.graph.update(graph1.graph)
    union.graph.update(graph2.graph)

    union.add_nodes_from(graph1)
    union.add_edges_from(graph1.edges(data=True))
    
    union.add_nodes_from(graph2)
    union.add_edges_from(graph2.edges(data=True))
    
    for n in graph1:
        union.nodes[n].update(graph1.nodes[n])
    for n in graph2:
        union.nodes[n].update(graph2.nodes[n])

    return union
