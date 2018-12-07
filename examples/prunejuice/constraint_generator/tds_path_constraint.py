#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from collections import *
from graphUtils import *

def delete_nth(d, n):
    d.rotate(-n)
    d.popleft()
    d.rotate(n)

class TdsPathConstraint():
    def __init__(self, subgraph, fromMerge, initialVertexList):
        
        self.size=subgraph.number_of_nodes()
        self.subgraph=subgraph.copy()
        self.fromMerge=fromMerge
        self.initialVertexList=list(initialVertexList)
                
    def __str__(self):
        return 'length : \t{}\nedge label list : \t{}\nedge vertex list : \t{}\ninitial vertex list : \t{}\n'.format(\
                str(self.size),\
                ' '.join(map(formatEdgeTuple, getGraphEdgeLabelList(self.subgraph))),\
                ' '.join(map(formatEdgeTuple, getGraphEdgeList(self.subgraph))),\
                ' '.join(map(str, self.initialVertexList)))
    def __rpr__(self):
        return self.__str__()
    
    def getCsvHeader():        
        return 'length ;edge label list; edge vertex list; initial vertex list\n'
        
    def getCsv(self):
        return '{};{};{};{}\n'.format(\
                str(self.size),\
                ' '.join(map(formatEdgeTuple, getGraphEdgeLabelList(self.subgraph))),\
                ' '.join(map(formatEdgeTuple, getGraphEdgeList(self.subgraph))),\
                ' '.join(map(str, self.initialVertexList)))


def generateTdsPathConstraint(pattern, pathConstraintList):
    tdsPathConstraintDict={}
    
    for pathConstraint in pathConstraintList:
        if not pathConstraint.initialLabel in tdsPathConstraintDict:
            tdsPathConstraintDict[pathConstraint.initialLabel]=\
                deque([TdsPathConstraint(pathConstraint.getSubgraph(pattern),False,\
                                         [pathConstraint.vertexList[0], pathConstraint.vertexList[-1]])])
        else:
            tdsPathConstraintDict[pathConstraint.initialLabel].append(\
                                 TdsPathConstraint(pathConstraint.getSubgraph(pattern),False,\
                                         [pathConstraint.vertexList[0], pathConstraint.vertexList[-1]]))
        
    for label in tdsPathConstraintDict:        
        currentIndex=0
        while currentIndex!= len(tdsPathConstraintDict[label]):
            nextIndex=currentIndex+1
            while nextIndex!= len(tdsPathConstraintDict[label]):
                omitList=list(set(tdsPathConstraintDict[label][currentIndex].initialVertexList\
                    + tdsPathConstraintDict[label][nextIndex].initialVertexList))
                if not isGraphVertexIntersectionEmpty(tdsPathConstraintDict[label][currentIndex].subgraph,\
                                             tdsPathConstraintDict[label][nextIndex].subgraph, omitList):
                    mergedSubgraph=graphUnion(tdsPathConstraintDict[label][currentIndex].subgraph,\
                                                tdsPathConstraintDict[label][nextIndex].subgraph)
                    delete_nth(tdsPathConstraintDict[label], nextIndex)
                    delete_nth(tdsPathConstraintDict[label], currentIndex)
                    tdsPathConstraintDict[label].appendleft(TdsPathConstraint(mergedSubgraph, True, omitList))
                    break
                nextIndex+=1   
            else:
                currentIndex+=1
                
    tdsPathConstraintList=[]    
    for label in tdsPathConstraintDict: 
        for tdsPathConstraint in tdsPathConstraintDict[label]:
            if  tdsPathConstraint.fromMerge:
                tdsPathConstraintList.append(tdsPathConstraint)            
    
    return tdsPathConstraintList

def writeTdsPathConstraint(outputResultDirectory, tdsPathConstraintList):
    outputResultTdsPathConstraint=outputResultDirectory + 'tds_path_constraint.txt'
    with open(outputResultTdsPathConstraint, 'w') as file:
        file.write(TdsPathConstraint.getCsvHeader())
        for tdsPathConstraint in tdsPathConstraintList:
            file.write(tdsPathConstraint.getCsv())
    file.close()   
    