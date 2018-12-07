#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from graphUtils import *
from collections import *

def delete_nth(d, n):
    d.rotate(-n)
    d.popleft()
    d.rotate(n)

class TdsSubtemplateConstraint():
    def __init__(self, subgraph, fromMerge):
        
        self.size=subgraph.number_of_nodes()
        self.subgraph=subgraph
        self.fromMerge=fromMerge
                
    def __str__(self):
        return 'length : \t{}\nedge label list : \t{}\nedge vertex list : \t{}\n'.format(\
                str(self.size),\
                ' '.join(map(formatEdgeTuple, getGraphEdgeLabelList(self.subgraph))),\
                ' '.join(map(formatEdgeTuple, getGraphEdgeList(self.subgraph))))
    def __rpr__(self):
        return self.__str__()
    
    def getCsvHeader():        
        return 'length ;edge label list; edge vertex list\n'
        
    def getCsv(self):
        return '{};{};{}\n'.format(\
                str(self.size),\
                ' '.join(map(formatEdgeTuple, getGraphEdgeLabelList(self.subgraph))),\
                ' '.join(map(formatEdgeTuple, getGraphEdgeList(self.subgraph))))

def generateTdsSubtemplateConstraint(tdsCyclicConstraintList, tdsPathConstraintList):
    tdsSubtemplateConstraintDeque=deque(\
                                   [TdsSubtemplateConstraint(tdsCyclicConstraint.subgraph, False)\
                                    for tdsCyclicConstraint in tdsCyclicConstraintList]+\
                                    [TdsSubtemplateConstraint(tdsPathConstraint.subgraph, False)\
                                    for tdsPathConstraint in tdsPathConstraintList])
    
    currentIndex=0
    while currentIndex!= len(tdsSubtemplateConstraintDeque):
        nextIndex=currentIndex+1
        while nextIndex!= len(tdsSubtemplateConstraintDeque):
            if not isGraphVertexIntersectionEmpty(tdsSubtemplateConstraintDeque[currentIndex].subgraph,\
                                         tdsSubtemplateConstraintDeque[nextIndex].subgraph):
                mergedSubgraph=graphUnion(tdsSubtemplateConstraintDeque[currentIndex].subgraph,\
                                            tdsSubtemplateConstraintDeque[nextIndex].subgraph)
                if not isGraphEdgeDifferenceEmpty(mergedSubgraph,tdsSubtemplateConstraintDeque[currentIndex].subgraph) \
                    and not isGraphEdgeDifferenceEmpty(mergedSubgraph,tdsSubtemplateConstraintDeque[nextIndex].subgraph):
                    delete_nth(tdsSubtemplateConstraintDeque, nextIndex)
                    delete_nth(tdsSubtemplateConstraintDeque, currentIndex)
                    tdsSubtemplateConstraintDeque.appendleft(TdsSubtemplateConstraint(mergedSubgraph, True))
                    break
            nextIndex+=1   
        else:
            currentIndex+=1
    
    currentIndex=0
    while currentIndex!= len(tdsSubtemplateConstraintDeque):
        if not tdsSubtemplateConstraintDeque[currentIndex].fromMerge:
            delete_nth(tdsSubtemplateConstraintDeque, currentIndex)            
        else:
            currentIndex+=1
    
    return list(tdsSubtemplateConstraintDeque)
    
def writeTdsSubtemplatecConstraint(outputResultDirectory, tdsSubtemplateConstraintList):
    outputResultTdsSubtemplateConstraint=outputResultDirectory + 'tds_subtemplate_constraint.txt'
    with open(outputResultTdsSubtemplateConstraint, 'w') as file:
        file.write(TdsSubtemplateConstraint.getCsvHeader())
        for tdsSubtemplateConstraint in tdsSubtemplateConstraintList:
            file.write(tdsSubtemplateConstraint.getCsv())
    file.close()   
