#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from collections import *


class LocalConstraint():
    def __init__(self, originalLabel, originalVertex, labelConstraintDict):
        self.originalLabel=originalLabel
        self.originalVertex=originalVertex
        self.labelList=list(labelConstraintDict.keys())
        self.labelNumberList=list(labelConstraintDict.values())
    
    def __str__(self):
        return 'original label : \t{}\noriginal vertex : \t{}\nlabel list : \t{}\nlabel number list : \t{}\n'.format(\
                self.originalLabel, self.originalVertex,\
                ' '.join(map(str, self.labelList)), ' '.join(map(str, self.labelNumberList)))
    def __rpr__(self):
        return self.__str__()
    
    def getCsvHeader():        
        return 'original label;original vertex;label list;label number list\n'
        
    def getCsv(self):
        return '{};{};{};{}\n'.format(\
                self.originalLabel, self.originalVertex,\
                ' '.join(map(str, self.labelList)), ' '.join(map(str, self.labelNumberList)))

def generateLocalConstraint(pattern):
    localConstraintList=[]
    for vertex in pattern.node():
        label=pattern.nodes[vertex]['label']
        
        labelConstraintDict=defaultdict(int) #Default value is 0
        for neighborVertex in pattern.neighbors(vertex):
            labelConstraintDict[pattern.nodes[neighborVertex]['label']]+=1
            
        localConstraintList.append(LocalConstraint(label, vertex, labelConstraintDict))
    return localConstraintList

def writeLocalConstraint(outputResultDirectory, localConstraintList):
    outputResultLocalConstraint=outputResultDirectory + 'local_constraint.txt'
    with open(outputResultLocalConstraint, 'w') as file:
        file.write(LocalConstraint.getCsvHeader())
        for localConstraint in localConstraintList:
            file.write(localConstraint.getCsv())
    file.close()
    