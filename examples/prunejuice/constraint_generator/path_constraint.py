#!/usr/bin/env python3
# -*- coding: utf-8 -*-

class PathConstraint():
    def __init__(self, labelList, vertexList):
        self.size=len(vertexList)
        self.initialLabel=labelList[0]
                
        # Re-order the vertex list such as :
        #  + vertexList[0]<vertexList[-1] 
        # This is used to test equality between constraint
            
        reverse=(vertexList[0]>vertexList[-1])
        if reverse :  
            self.labelList=list(reversed(labelList))
            self.vertexList=list(reversed(vertexList))
        else:
            self.labelList=list(labelList)
            self.vertexList=list(vertexList)
            
    def getSubgraph(self, pattern):
        return  pattern.subgraph(self.vertexList)
                
    def __eq__(self, other):
        return self.size==other.size and self.vertexList==other.vertexList
    def __ne__(self, other):
        return not __eq__(self, other)
    def __hash__(self):
        return hash(self.vertexList)
    
    
    def __str__(self):
        return 'length : \t{}\ninitial label : \t{}\nlabel list : \t{}\nvertex list : \t{}\n'.format(\
                str(self.size), str(self.initialLabel),\
                ' '.join(map(str, self.labelList)), ' '.join(map(str, self.vertexList)))
    def __rpr__(self):
        return self.__str__()
    
    def getCsvHeader():        
        return 'length ; initial label ; label list ; vertex list\n'
        
    def getCsv(self):
        return '{};{};{};{}\n'.format(\
                str(self.size), str(self.initialLabel),\
                ' '.join(map(str, self.labelList)), ' '.join(map(str, self.vertexList)))
        
        
def runPathConstraintSearch(pattern, pathConstraintList, originalVertex, \
                              currentVertex, historyLabelList, historyVertexList, currentDepth):
    if currentDepth>pattern.number_of_nodes():
        return
    
    if currentDepth>=2 :
        labelOriginalVertex=pattern.nodes[originalVertex]['label']
        labelCurrentVertex=pattern.nodes[currentVertex]['label']
        if labelOriginalVertex==labelCurrentVertex:
            pathConstraint = PathConstraint(historyLabelList, historyVertexList)
            
            if not (pathConstraint in pathConstraintList):
                pathConstraintList.append(pathConstraint)
                
            return
    
    for neighborVertex in pattern.neighbors(currentVertex):
        # Avoid back-tracking
        if neighborVertex in historyVertexList:
            continue
        
        neighborLabel=pattern.nodes[neighborVertex]['label']
        
        historyVertexList.append(neighborVertex)
        historyLabelList.append(neighborLabel)
        
        runPathConstraintSearch(pattern, pathConstraintList, originalVertex, \
                                  neighborVertex, historyLabelList, historyVertexList, currentDepth+1)
        
        historyVertexList.pop()
        historyLabelList.pop()
    
    
    
def generatePathConstraint(pattern):
    pathConstraintList=[]
    
    for vertex in pattern.node():
        label=pattern.nodes[vertex]['label']
        
        historyLabelList=[label]
        historyVertexList=[vertex]
        
        runPathConstraintSearch(pattern, pathConstraintList, vertex, \
                                  vertex, historyLabelList, historyVertexList, 0)
        
    return pathConstraintList

def writePathConstraint(outputResultDirectory, pathConstraintList):
    outputResultPathConstraint=outputResultDirectory + 'path_constraint.txt'
    with open(outputResultPathConstraint, 'w') as file:
        file.write(PathConstraint.getCsvHeader())
        for pathConstraint in pathConstraintList:
            file.write(pathConstraint.getCsv())
    file.close()