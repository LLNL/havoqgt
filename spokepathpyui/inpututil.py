import os
import sys

import json
import math
import re

from collections import defaultdict
from pprint import pformat

import matplotlib.pyplot as plt

import networkx as nx
from networkx.algorithms import bipartite

import ipywidgets as widgets
from ipywidgets import Output

from IPython.display import display, HTML, IFrame

import ipycytoscape

out = Output()

print("inpututil.py")

# vertex data

include_json_string = True

vertex_type_map = dict()
vertex_data_map = defaultdict(tuple)

def parse_json_data_2(json_data, line, include_json_string = False) :
    id_int = int(json_data["id"])
    label = json_data["labels"][0]
    properties_identifier = ""
    properties_name = ""
    
    if "properties" in json_data :
        if "identifier" in json_data["properties"] :
            properties_identifier = json_data["properties"]["identifier"]
        if "name" in json_data["properties"] :
            properties_name = json_data["properties"]["name"][:15]
    
    # Test
    # to save memory
    if include_json_string == False :
        line = ""
        
    vertex_data_map[id_int] = (label, properties_identifier, properties_name, line)
    
    vertex_type_count = vertex_type_map.get(label, 0)
    vertex_type_map[label] = vertex_type_count + 1

def load_data(vertex_data_json_filename) :

	vertex_type_map.clear()
	vertex_data_map.clear()

	with open(vertex_data_json_filename, "r") as read_file :
		line_count = 0
		for line in read_file :
			line = line.strip()
			try :
			    json_data = json.loads(line)
			    parse_json_data_2(json_data, line, include_json_string)
			except Exception as e :
			    print(e)
			    continue
			line_count = line_count + 1
		 #print(str(line_count))

	print(len(vertex_data_map))
	print(len(vertex_type_map))

def get_vertex_type_count_list(vertex_type_map) :
    tuple_list = []
    for k, v in vertex_type_map.items() :
        tuple_list.append((str(k) + " (" + str(v) + ")", k))
    return tuple_list
    #print(tuple_list)

def vertices_filtered_by_type(vertex_type, begin = 0, end = 25) :
    #filtered_veretx_data_map = dict()
    filtered_veretx_data_list = []
    for k, v in vertex_data_map.items() :
        if v[0] != vertex_type :
            continue
        else :
            #filtered_veretx_data_map[k] = v
            #filtered_veretx_data_map[(v, k)] = v
            
            if v[2] != "" :
                #filtered_veretx_data_map[k] = v[2]
                filtered_veretx_data_list.append((str(k) + " (" + str(v[2]) + ")", k))
            elif v[1] != "" :
                #filtered_veretx_data_map[k] = v[1]
                filtered_veretx_data_list.append((str(k) + " (" + str(v[1]) + ")", k))
            else :
                #filtered_veretx_data_map[k] = str(k) #v[0]
                filtered_veretx_data_list.append((str(k) + " (" + str(v[2]) + ")", k))
                
    #return filtered_veretx_data_map
    #return filtered_veretx_data_list

    if (end > len(filtered_veretx_data_list)) :
        end = len(filtered_veretx_data_list)
    return filtered_veretx_data_list[begin : end], end, len(filtered_veretx_data_list)

def kw_search_filtered_by_veretx_type(vertex_type_set, keyword, result_list) :
    for k, v in vertex_data_map.items() :        
        if v[0] not in vertex_type_set :
            continue
        match_found = False
        
        complete_json_string = v[3]
        
        match_found = re.search(keyword, complete_json_string, re.IGNORECASE)
        
        if match_found :
            #print(v[3])
            if v[2] != "" :                
                result_list.append((str(k) + " (" + str(v[2]) + ")", k))
            elif v[1] != "" :                
                result_list.append((str(k) + " (" + str(v[1]) + ")", k))
            else :                
                result_list.append((str(k) + " (" + str(v[2]) + ")", k))               
    
    return len(result_list)

def load_more_from_list_by_range(item_list, begin = 0, end = 25) :
    if (end > len(item_list)) :
        end = len(item_list)
    return item_list[begin : end], end, len(item_list)

# # # # # #

# input parameters

s_t_vertices = []
filter_vertex_labels = []

def setup_input_parameters(vertex_data_filename, vertex_type_a_set, vertex_type_b_set, vertex_type_c_set) :
    vertex_type_c_map = dict()

    with open (vertex_data_filename, "r") as read_file :
        for line in read_file :
            if len(vertex_type_c_set) < 1 :
                break
            line = line.strip()
            tokens = line.split(" ")
            vertex = int(tokens[0].strip(' \t\n\r'))
            vertex_data_hash = int(tokens[1].strip(' \t\n\r'))
            vertex_data_string = str(tokens[2].strip(' \t\n\r'))
            if vertex_data_string in vertex_type_c_map :
                continue
            if vertex_data_string in vertex_type_c_set :
                vertex_type_c_map[vertex_data_string] = vertex_data_hash
            if len(vertex_type_c_map) == len(vertex_type_c_set) :
                break

    if len(vertex_type_c_map) != len(vertex_type_c_set) :
        print("Error: items not found.")
        return
    print(vertex_type_c_map)

    for i in vertex_type_a_set :
        s_t_vertices.append(i)

    for i in vertex_type_b_set :
        s_t_vertices.append(i)

    for k, v in vertex_type_c_map.items() :
        filter_vertex_labels.append(v)

    print(s_t_vertices)
    print(filter_vertex_labels)

