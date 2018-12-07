#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import networkx as nx
from collections import *

def readPatternEdges(inputPatternEdges):
    pattern=nx.DiGraph()
    with open(inputPatternEdges) as file:
        for line in file:
            vertexList=line.replace("\n", "").split(' ')
            if len(vertexList)==2:
                pattern.add_edge(int(vertexList[0]), int(vertexList[1]))
    file.close()
    return pattern

def readPatternVertexData(pattern, inputPatternVertexData):
    with open(inputPatternVertexData) as file:
        for line in file:
            vertexList=line.replace("\n", "").split(' ')
            if len(vertexList)==2:
                pattern.nodes[int(vertexList[0])]['label']=vertexList[1]
    file.close()

def findLeafVertexWithUniqueLabel(pattern):
    labelCount=defaultdict(int) #Default value is 0
    
    for vertex in pattern.node():
        label=pattern.nodes[vertex]['label']
        
        labelCount[label]+=1
        
    leafVertexWithUniqueLabelList=[]
    for vertex in pattern.node():
        if(labelCount[pattern.nodes[vertex]['label']]==1 and len(list(pattern.neighbors(vertex)))==1):
            leafVertexWithUniqueLabelList.append(vertex)
    
    return leafVertexWithUniqueLabelList

