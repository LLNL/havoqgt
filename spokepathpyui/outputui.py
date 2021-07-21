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

from spokepathpyui import basicgraph
from spokepathpyui import inpututil

out = Output()

print("outputui.py")

# output

def basic_graph_to_cytoscape_graph (graph) :
    graph_json = {"nodes" : [], "edges" : []}
   
    for v in graph.vertices :
            #print(v)
            v_dict = {"data" : {}}
            v_dict["data"]["id"] = str(v)
            v_dict["data"]["idInt"] = v
            #v_dict["data"]["label"] = str(v)
            v_dict["data"]["type"] = str(inpututil.vertex_data_map[v][0])
            v_dict["data"]["label"] = str(inpututil.vertex_data_map[v][2])
            v_dict["data"]["description"] = str(inpututil.vertex_data_map[v][3])
            graph_json["nodes"].append(v_dict)
        
    for s, t in graph.edges :
            #print(s, end = " ")
            #print(t)
            e_dict = {"data" : {}}
            #e_dict["data"]["id"] = ""
            e_dict["data"]["source"] = str(s)
            e_dict["data"]["target"] = str(t)
            graph_json["edges"].append(e_dict)

    #print(graph_json)
    return graph_json

# # # # # #

node_info_widget = widgets.Textarea(value="", placeholder="Empty", disabled=True, 
                                   layout=widgets.Layout(width='50%', height='200px', align_content='center'))

def reset_node_info_widget() :
	node_info_widget.value = "Empty"

def update_node_info_widget(node) :
    #node_info_widget.value = node["data"]["label"] 
    vertex = int(node["data"]["id"])
    #vertex_data = dict(sorted(vertex_data_map[vertex].items()))
    if str(inpututil.vertex_data_map[vertex][3]) != "" :
        vertex_data = inpututil.vertex_data_map[vertex][3]
    else :
        vertex_data = inpututil.vertex_data_map[vertex]
    node_info_widget.value = str(vertex_data)
    #output_string = ""
    #for k, v in vertex_data.items() :
    #    output_string = output_string + str(k) + " : " + str(v) + "\n"
    #node_info_widget.value = output_string

def cytoscape_node_callback_click(node) :
	update_node_info_widget(node)

def cytoscape_graph_set_style (cytoscape_object) :
    cytoscape_object.set_style([
    {
      "selector": "edge.highlighted",
      "css": {
        "line-color": "red"
      }
    },
    {
      "selector": "node.highlighted",
      "css": {
        "background-color": "yellow",
        "border-color": "black",
        "border-width": "5px",
        "border-opacity": "0.3",  
      },
    },
    {
      "selector": "node.Anatomy",
      "css": {
        "background-color": "rgb(179, 222, 105)"
      },
    },    
    {
      "selector": "node.BiologicalProcess",
      "css": {
        "background-color": "rgb(253, 180, 98)"
      },
    },
    {
      "selector": "node.CellularComponent",
      "css": {
        "background-color": "rgb(255, 255, 179)"
      },
    },        
    {
      "selector": "node.Compound",
      "css": {
        "background-color": "rgb(188, 128, 189)"
      },
    },
    {
      "selector": "node.Disease",
      "css": {
        "background-color": "rgb(251, 128, 114)"
      },
    },       
    {
      "selector": "node.EC",
      "css": {
        "background-color": "rgb(247, 151, 103)"
      },
    },  
    {
      "selector": "node.Food",
      "css": {
        "background-color": "rgb(199, 183, 143)"
      },
    },      
    {
      "selector": "node.Gene",
      "css": {
        "background-color": "rgb(128, 177, 211)"
      },
    },
    {
      "selector": "node.MolecularFunction",
      "css": {
        "background-color": "rgb(255, 237, 111)"
      },
    },    
    {
      "selector": "node.Organism",
      "css": {
        "background-color": "rgb(217, 200, 174)"
      },
    },    
    {
      "selector": "node.Pathway",
      "css": {
        "background-color": "rgb(87, 199, 227)"
      },
    },
    {
      "selector": "node.PharmacologicClass",
      "css": {
        "background-color": "rgb(190, 186, 218)"
      },
    }, 
    {
      "selector": "node.Protein",
      "css": {
        "background-color": "rgb(141, 211, 199)"
      },
    },     
    {
      "selector": "node.Reaction",
      "css": {
        "background-color": "rgb(241, 102, 103)"
      },
    },
    {
      "selector": "node.SideEffect",
      "css": {
        "background-color": "rgb(204, 235, 197)"
      },
    },
    {
      "selector": "node.Symptom",
      "css": {
        "background-color": "rgb(252, 205, 229)"
      },
    },
    {
      "selector": "node.CellType",
      "css": {
        "background-color": "green"
      },
    },
    {
      "selector": "node.Nutrient",
      "css": {
        "background-color": "cyan"
      },
    },
    {
      "selector": "node.MolecularFunction",
      "css": {
        "background-color": "darkcyan"
      },
    },    
    {
      "selector": "node.SARSCov2",
      "css": {
        "background-color": "crimson"
      },
    },
    {
      "selector": "node.DatabaseTimestamp",
      "css": {
        "background-color": "darkgrey"
      },
    },
    {
      "selector": "node.Version",
      "css": {
        "background-color": "darkgrey"
      },
    },    
    {"selector": "node", "style": {"content": "data(label)", "text-valign": "center", "text-halign": "center", "font-size": "8px"}}, 
        {"selector": "node[classes]", "style": {"label": "data(label)"}}])
    
#"background-color": "grey", "content" ...

def set_node_color_by_type (cytoscape_object, target_node, node_type) :
    for node in cytoscape_object.graph.nodes:
        if int(node.data["id"]) == int(target_node["data"]["id"]) :
            classes = set(node.classes.split(" "))    
            classes.add(node_type)
            node.classes = " ".join(classes)
            break

def set_node_color (cytoscape_object) :
    for node in cytoscape_object.graph.nodes :        
        #print(node.data["id"])
        #print(node.data["type"])
        set_node_color_by_type(cytoscape_object, {"data" : {"id": node.data["id"]}}, node.data["type"])
    #set_node_color_by_type(cytoscape_object, {"data" : {"id": "2276338"}}, "compound")
    #set_node_color_by_type(cytoscape_object, {"data" : {"id": "2293366"}}, "reaction")

def unhighlight_nodes (cytoscape_object) :
    for node in cytoscape_object.graph.nodes:
        classes = set(node.classes.split(" "))
        if "highlighted" in classes :
            classes.remove("highlighted")
            node.classes = " ".join(classes)            

def unhighlight_edges (cytoscape_object) :
    for edge in cytoscape_object.graph.edges:
        classes = set(edge.classes.split(" "))
        if "highlighted" in classes :
            classes.remove("highlighted")
            edge.classes = " ".join(classes)                

def highlight_node (cytoscape_object, target_node) :           
    for node in cytoscape_object.graph.nodes:
        if int(node.data["id"]) == int(target_node["data"]["id"]) :
            classes = set(node.classes.split(" "))    
            classes.add("highlighted")
            node.classes = " ".join(classes)
            break
            
def highlight_edges (cytoscape_object, edge_filter) :
    for edge in cytoscape_object.graph.edges:
        #print(edge.data)
        #print(edge_filter)
    
        s = int(edge.data["source"]) 
        t = int(edge.data["target"])             
    
        if (s, t) in edge_filter :
            classes = set(edge.classes.split(" "))
            classes.add("highlighted")
            #classes.remove("highlighted")            
            edge.classes = " ".join(classes)
            continue
        
        #classes = set(edge.classes.split(" "))
        #classes.add("highlighted")
        #edge.classes = " ".join(classes)
    
        #print(edge)

def draw_cytoscape_graph (edge_list, s_t_vertices = list(), cytoscape_graph_layout_name = "cola") :

	basic_graph = basicgraph.BasicGraph()
	basic_graph.read_edge_list(edge_list)	

	cytoscape_graph_json = basic_graph_to_cytoscape_graph(basic_graph)
	cytoscape_object = ipycytoscape.CytoscapeWidget()
	cytoscape_object.graph.add_graph_from_json(cytoscape_graph_json)
	
	cytoscape_object.set_layout(name=cytoscape_graph_layout_name)
	#cytoscapeobj.set_layout(nodeSpacing=100)
	#cytoscapeobj.get_layout()
	#cytoscape_object.set_tooltip_source("label")
	cytoscape_object.set_tooltip_source("description")

	cytoscape_object.on('node', 'click', cytoscape_node_callback_click)

	set_node_color(cytoscape_object)

	for i in range(len(s_t_vertices)) :
		highlight_node(cytoscape_object, {"data" : {"id": s_t_vertices[i]}})

	reset_node_info_widget()

	display(cytoscape_object)
	display(node_info_widget)

