#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import argparse 

import utils
import local_constraint
import cycle_constraint
import path_constraint
import tds_cyclic_constraint
import tds_path_constraint
import tds_subtemplate_constraint

def main(inputPatternEdges, inputPatternVertexData, outputResultDirectory):
    pattern=utils.readPatternEdges(inputPatternEdges)
    utils.readPatternVertexData(pattern, inputPatternVertexData)
    
    # Generate local constraint
    localConstraintList=local_constraint.generateLocalConstraint(pattern)
    local_constraint.writeLocalConstraint(outputResultDirectory, localConstraintList)
    
    # Find leaf with unique label and prune the pattern
    leafVertexWithUniqueLabelList=utils.findLeafVertexWithUniqueLabel(pattern)
    pattern.remove_nodes_from(leafVertexWithUniqueLabelList)
    
    # Generate cycle constraint
    cycleConstraintList=cycle_constraint.generateCycleConstraint(pattern)
    cycle_constraint.writeCycleConstraint(outputResultDirectory, cycleConstraintList)
    
    # Generate path constraint
    pathConstraintList=path_constraint.generatePathConstraint(pattern)
    path_constraint.writePathConstraint(outputResultDirectory, pathConstraintList)
    
    # Generate TDS edge monocyclic constraint
    tdsCyclicConstraintList=tds_cyclic_constraint.generateTdsCyclicConstraint(pattern, cycleConstraintList)
    tds_cyclic_constraint.writeTdsCyclicConstraint(outputResultDirectory, tdsCyclicConstraintList)
    
    # Generate TDS path constraint
    tdsPathConstraintList=tds_path_constraint.generateTdsPathConstraint(pattern, pathConstraintList)
    tds_path_constraint.writeTdsPathConstraint(outputResultDirectory, tdsPathConstraintList)
    
    # Generate TDS partial constraint
    tdsSubtemplateConstraintList=tds_subtemplate_constraint.generateTdsSubtemplateConstraint(tdsCyclicConstraintList, tdsPathConstraintList)
    tds_subtemplate_constraint.writeTdsSubtemplatecConstraint(outputResultDirectory, tdsSubtemplateConstraintList)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='main.py', description='Generate constraint for pattern matching.')
    parser.add_argument('-ie', '--input_pattern_edges', help='Input pattern edges', required=True)
    parser.add_argument('-id', '--input_pattern_data', help='Input pattern vertex data', required=True)
    parser.add_argument('-o', '--output_directory', help='Output directory', required=True)
    try :
        args=parser.parse_args()
        
        inputPatternEdges=args.input_pattern_edges
        inputPatternVertexData=args.input_pattern_data
        outputDirectory=args.output_directory
        if(outputDirectory[-1]!='/' or outputDirectory[-1]!='\\' ):
            outputDirectory=outputDirectory+'/'
            
        main(inputPatternEdges, inputPatternVertexData, outputDirectory)
    except SystemExit:
        pass
        
    
