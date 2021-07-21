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

from spokepathpyui import inpututil

out = Output()

print("inputui.py")

def print_string(str) :
	print("UI " + str)

# input widgets
#vertex_type_count_list = inpututil.get_vertex_type_count_list(inpututil.vertex_type_map)

vertex_type_selected_a = set()
vertex_type_selected_b = set()

vertex_type_a_set = set() # TODO: change name, the set contains vertex IDs
vertex_type_b_set = set() # TODO: cange name
vertex_type_c_set = set()

sub_drdn_btn_pair_map = defaultdict(list)    

hbox_vertex_data_a = widgets.HBox([])
hbox_vertex_data_b = widgets.HBox([])

label_selected_vertices_a = widgets.Label("Empty")
label_selected_vertices_b = widgets.Label("Empty")
label_selected_vertex_types = widgets.Label("Empty")

txar_vertex_info = widgets.Textarea(value="", description="Vertex info", placeholder="Empty", disabled=True, 
                                    layout=widgets.Layout(width='30%', height='200px', align_content='center'))

# search
load_more_count = 50

search_result_a_list = list()
search_result_b_list = list()

text_kw_search_a = widgets.Text(value="", placeholder="Empty")
text_kw_search_b = widgets.Text(value="", placeholder="Empty")

label_kw_search_result_count_a = widgets.Label("")
label_kw_search_result_count_b = widgets.Label("")

# # # # # #

def txar_vertex_info_change(vertex) :
    if str(inpututil.vertex_data_map[int(vertex)][3]) != "" :
        txar_vertex_info.value = str(inpututil.vertex_data_map[int(vertex)][3]) # json string
    else :
        txar_vertex_info.value = str(inpututil.vertex_data_map[int(vertex)])        
    #print(inpututil.vertex_data_map[int(vertex)])

def btn_reset_all_on_click(*args) :
    vertex_type_selected_a.clear()
    vertex_type_selected_b.clear()
    
    vertex_type_a_set.clear()
    vertex_type_b_set.clear()
    vertex_type_c_set.clear()
  
    hbox_vertex_data_a.children = tuple()
    hbox_vertex_data_b.children = tuple()

    txar_vertex_info.value = "Empty"
    
    label_selected_vertices_a.value = "Empty"
    label_selected_vertices_b.value = "Empty"
    label_selected_vertex_types.value = "Empty"
    
    sub_drdn_btn_pair_map.clear()
    
    text_kw_search_a.value = ""
    text_kw_search_b.value = ""
    
    label_kw_search_result_count_a.value = ""
    label_kw_search_result_count_b.value = ""
    
    search_result_a_list.clear()
    search_result_b_list.clear()
    
def btn_selected_vertices_a_on_click(*args) :
    # args[0] is the object itself
    vertex_type_a_set.clear()
    label_selected_vertices_a.value = "Empty"    
    txar_vertex_info.value = "Empty"

def btn_selected_vertices_b_on_click(*args) :
    vertex_type_b_set.clear()
    label_selected_vertices_b.value = "Empty"    
    txar_vertex_info.value = "Empty"

def btn_selected_vertex_types_on_click(*args) :
    vertex_type_c_set.clear()
    label_selected_vertex_types.value = "Empty"    
        
def on_sub_drdn_value_change_a(change) :
    #print("a: " + str(change["type"]))
    #print("a: " + str(change["new"]))
    #print("a: " + str(len(change["new"])))
    
    if len(change["new"]) < 1 :
        txar_vertex_info.value = "Empty"
        return
    
    for i in change["new"] :
        vertex_type_a_set.add(i)
    #print(vertex_type_a_set)    
    label_selected_vertices_a.value = str(sorted(vertex_type_a_set))
    txar_vertex_info_change(change["new"][0])
    
def on_sub_drdn_value_change_b(change) :
    #print("b: " + str(change["new"]))
    
    if len(change["new"]) < 1 :
        txar_vertex_info.value = "Empty"
        return
    
    for i in change["new"] :
        vertex_type_b_set.add(i)
    #print(vertex_type_b_set)
    label_selected_vertices_b.value = str(sorted(vertex_type_b_set))
    txar_vertex_info_change(change["new"][0])

def btn_search_result_load_more_on_click(*args) :
    btn = args[0]
    #print(id(btn))
    #print(sub_drdn_btn_pair_map)
    
    if sub_drdn_btn_pair_map[id(btn)][1] >= sub_drdn_btn_pair_map[id(btn)][2] :
        return
    
    drdn_found = False
    
    for vbox in hbox_vertex_data_a.children :
        for child in vbox.children : 
            if sub_drdn_btn_pair_map[id(btn)][0] == id(child) :
                #print(sub_drdn_btn_pair_map)
        
                vertex_type = sub_drdn_btn_pair_map[id(btn)][3]
                begin = int(sub_drdn_btn_pair_map[id(btn)][1]) # starting point for this batch
                end = begin + load_more_count # (begin : end]
                
                #print(vertex_type)
                #print(begin)
                #print(str(end - 1)) # (begin : end-1)
                #print(len(search_result_a_list))
                
                sub_drdn_options_list, sub_drdn_options_list_end, sub_drdn_options_list_len = \
                    inpututil.load_more_from_list_by_range(search_result_a_list, begin, end)
                
                #print(len(sub_drdn_options_list))
                #print(sub_drdn_options_list_end)
                #print(sub_drdn_options_list_len)
                
                child.options = sub_drdn_options_list # replace the previous entries
                
                sub_drdn_btn_pair_map[id(btn)][1] = sub_drdn_options_list_end # starting point for the next batch                
                                
                btn.description = "Load more " + str(begin + 1) + "-" + str(sub_drdn_options_list_end)            
             
                drdn_found = True
                break
            if drdn_found :
                break
                
    for vbox in hbox_vertex_data_b.children :
        for child in vbox.children : 
            if sub_drdn_btn_pair_map[id(btn)][0] == id(child) :
                        
                vertex_type = sub_drdn_btn_pair_map[id(btn)][3]
                begin = int(sub_drdn_btn_pair_map[id(btn)][1]) # starting point for this batch
                end = begin + load_more_count # (begin : end]
                                                
                sub_drdn_options_list, sub_drdn_options_list_end, sub_drdn_options_list_len = \
                    inpututil.load_more_from_list_by_range(search_result_b_list, begin, end)
                                                
                child.options = sub_drdn_options_list # replace the previous entries
                sub_drdn_btn_pair_map[id(btn)][1] = sub_drdn_options_list_end # starting point for the next batch                                
                btn.description = "Load more " + str(begin + 1) + "-" + str(sub_drdn_options_list_end)            
             
                drdn_found = True
                break
            if drdn_found :
                break

def btn_kw_search_a_on_click(*args) :
    hbox_vertex_data_a.children = tuple()
    label_kw_search_result_count_a.value = ""
    label_selected_vertices_a.value = "Empty"
    txar_vertex_info.value = "Empty"    
    search_result_a_list.clear()
    
    if text_kw_search_a.value == "" :
        return
    vertex_type = "dummy" # TODO: remove
    #print(text_kw_search_a.value)    
    #print(vertex_type_selected_a)
    inpututil.kw_search_filtered_by_veretx_type(vertex_type_selected_a, text_kw_search_a.value, search_result_a_list)
    label_kw_search_result_count_a.value = str(len(search_result_a_list))
    #print(len(search_result_a_list))
    
    sub_drdn_options_list, sub_drdn_options_list_end, sub_drdn_options_list_len = \
        inpututil.load_more_from_list_by_range(search_result_a_list, 0, load_more_count)
    
    sub_drdn = widgets.SelectMultiple(options=sub_drdn_options_list, description="")        
    sub_drdn.observe(on_sub_drdn_value_change_a, names='value')
    
    sub_btn = widgets.Button(description="Load more 1-" + str(sub_drdn_options_list_end), tooltip="Tooltip", disabled=False)
    sub_btn.on_click(callback=btn_search_result_load_more_on_click)
               
    vbox_sub_drdn_btn = widgets.VBox([sub_drdn, sub_btn])
    
    #hbox_vertex_data_a.children = hbox_vertex_data_a.children + (sub_drdn,)
    hbox_vertex_data_a.children = hbox_vertex_data_a.children + (vbox_sub_drdn_btn,)
    
    sub_drdn_btn_pair_map[id(sub_btn)] = [id(sub_drdn), sub_drdn_options_list_end, sub_drdn_options_list_len, vertex_type]

def btn_kw_search_b_on_click(*args) :
    hbox_vertex_data_b.children = tuple()
    label_kw_search_result_count_b.value = ""
    label_selected_vertices_b.value = "Empty"
    txar_vertex_info.value = "Empty"    
    search_result_b_list.clear()
    
    if text_kw_search_b.value == "" :
        return
    vertex_type = "dummy" # TODO: remove
    #print(text_kw_search_b.value)    
    #print(vertex_type_selected_b)
    inpututil.kw_search_filtered_by_veretx_type(vertex_type_selected_b, text_kw_search_b.value, search_result_b_list)
    label_kw_search_result_count_b.value = str(len(search_result_b_list))
    #print(len(search_result_b_list))
    
    sub_drdn_options_list, sub_drdn_options_list_end, sub_drdn_options_list_len = \
        inpututil.load_more_from_list_by_range(search_result_b_list, 0, load_more_count)
    
    sub_drdn = widgets.SelectMultiple(options=sub_drdn_options_list, description="")        
    sub_drdn.observe(on_sub_drdn_value_change_b, names='value')
    
    sub_btn = widgets.Button(description="Load more 1-" + str(sub_drdn_options_list_end), tooltip="Tooltip", disabled=False)
    sub_btn.on_click(callback=btn_search_result_load_more_on_click)
               
    vbox_sub_drdn_btn = widgets.VBox([sub_drdn, sub_btn])
    
    #hbox_vertex_data_b.children = hbox_vertex_data_b.children + (sub_drdn,)
    hbox_vertex_data_b.children = hbox_vertex_data_b.children + (vbox_sub_drdn_btn,)
    
    sub_drdn_btn_pair_map[id(sub_btn)] = [id(sub_drdn), sub_drdn_options_list_end, sub_drdn_options_list_len, vertex_type]
    
def on_drdn_value_change_a(change) :
    #print(change["new"])    
    #print(change)
    hbox_vertex_data_a.children = tuple()
    text_kw_search_a.value = ""
    label_kw_search_result_count_a.value = ""
    label_selected_vertices_a.value = "Empty"
    txar_vertex_info.value = "Empty"    
    vertex_type_selected_a.clear()
    vertex_type_a_set.clear()
    
    for i in change["new"] :
        vertex_type = str(i)        
        vertex_type_selected_a.add(vertex_type)              
    
def on_drdn_value_change_b(change) :
    #print(change["new"])    
    #print(change)
    hbox_vertex_data_b.children = tuple()    
    text_kw_search_b.value = ""
    label_kw_search_result_count_b.value = ""
    label_selected_vertices_b.value = "Empty"
    txar_vertex_info.value = "Empty"    
    vertex_type_selected_b.clear()
    vertex_type_b_set.clear()
    
    for i in change["new"] :
        vertex_type = str(i)
        vertex_type_selected_b.add(vertex_type)       
        
def on_drdn_value_change_c(change) :
        #print(change["new"])  
        for i in change["new"] :
            vertex_type_c_set.add(i)
        label_selected_vertex_types.value = str(sorted(vertex_type_c_set))

# # # # # #

def init_inputui() :

	vertex_type_count_list = inpututil.get_vertex_type_count_list(inpututil.vertex_type_map)
        
	#drdn_vertex_type_a = widgets.SelectMultiple(options=[*vertex_type_map], description="Source types")
	drdn_vertex_type_a = widgets.SelectMultiple(options=vertex_type_count_list, description="Source types")
	drdn_vertex_type_a.observe(on_drdn_value_change_a, names='value')

	#drdn_vertex_type_b = widgets.SelectMultiple(options=[*vertex_type_map], description="Target types")
	drdn_vertex_type_b = widgets.SelectMultiple(options=vertex_type_count_list, description="Target types")
	drdn_vertex_type_b.observe(on_drdn_value_change_b, names='value')

	btn_kw_search_a = widgets.Button(description="Search")
	btn_kw_search_a.on_click(callback=btn_kw_search_a_on_click)

	btn_kw_search_b = widgets.Button(description="Search")
	btn_kw_search_b.on_click(callback=btn_kw_search_b_on_click)

	btn_selected_vertices_a = widgets.Button(description="Reset", tooltip="Tooltip", disabled=False)
	btn_selected_vertices_a.on_click(callback=btn_selected_vertices_a_on_click)

	btn_selected_vertices_b = widgets.Button(description="Reset", tooltip="Tooltip", disabled=False)
	btn_selected_vertices_b.on_click(callback=btn_selected_vertices_b_on_click)

	btn_reset_all = widgets.Button(description="Reset", disabled=False)
	btn_reset_all.on_click(callback=btn_reset_all_on_click)

	hbox_vertex_type_ab = widgets.HBox([drdn_vertex_type_a, drdn_vertex_type_b, txar_vertex_info, btn_reset_all])

	# filter by vertex type
	drdn_vertex_type_c = widgets.SelectMultiple(options=vertex_type_count_list, description="Filter by vertex type")
	drdn_vertex_type_c.observe(on_drdn_value_change_c, names='value')

	btn_selected_vertex_types = widgets.Button(description="Reset", disabled=False)
	btn_selected_vertex_types.on_click(callback=btn_selected_vertex_types_on_click)

	#print(widgets.Widget.observe.__doc__)
	print(len(inpututil.vertex_data_map))	
	print(len(vertex_type_count_list))

	# display input widgets
	display(hbox_vertex_type_ab)

	display(widgets.Label("Source vertices"))
	display(widgets.HBox([btn_kw_search_a, text_kw_search_a, label_kw_search_result_count_a]))
	display(hbox_vertex_data_a)

	display(widgets.HBox([btn_selected_vertices_a, widgets.Label("Selected source vertices"), label_selected_vertices_a]))

	display(widgets.Label("Target vertices"))
	display(widgets.HBox([btn_kw_search_b, text_kw_search_b, label_kw_search_result_count_b]))
	display(hbox_vertex_data_b)
	display(widgets.HBox([btn_selected_vertices_b, widgets.Label("Selected target vertices"), label_selected_vertices_b]))

	display(widgets.HBox([drdn_vertex_type_c, btn_selected_vertex_types, widgets.Label("Selected vertex types"), label_selected_vertex_types]))	
