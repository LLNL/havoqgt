import glob
import ijson
import json
import os
import sys
#import thread
from collections import defaultdict
from collections import OrderedDict
#from urlparse import urlparse

maxvid_json_data = 0
maxeid_json_data = 0 

def parse_json_data_2 (json_data, json_string, vertex_data, edge_data) :
	if json_data["type"] == "node" :
		vid = int(json_data["id"])
		labels = json_data["labels"][0]
		#vertex_data[vid] = (labels, json_string)
		vertex_data[vid] = (json_string, labels)

		global maxvid_json_data
		if vid > maxvid_json_data :
			maxvid_json_data = vid
			
	elif json_data["type"] == "relationship" :
		eid = int(json_data["id"])
		#label = json_data["label"]
		elabel= json_data["label"]
		startid = json_data["start"]["id"]
		starlabels = json_data["start"]["labels"][0] 
		endid = json_data["end"]["id"]	
		endlabels = json_data["end"]["labels"][0]
		#edge_data[eid] = (startid, endid, label, json_string)
		edge_data[eid] = (startid, endid, starlabels, endlabels, elabel, json_string)

		global maxeid_json_data
		if eid > maxeid_json_data :
			maxeid_json_data = eid		

def filter_edges (vertex_data, edge_data, required_vertex_data_map, adjacency_list, 
	cmp_rec_adjacency_list, rec_cmp_adjacency_list, rec_active_list) :
	for key, value in edge_data.items() :
		#edge_data[eid] = (startid, endid, starlabels, endlabels, elabel, json_string)
		startid, endid, starlabels, endlabels, elabel, json_string = value #edge_data[key]

		s = int(startid)
		s_data = starlabels
		t = int(endid)
		t_data = endlabels				

		if s_data in required_vertex_data_map and t_data in required_vertex_data_map :
			
			if s_data == "Organism" :
				if t_data == "EC" :
					adjacency_list[s].append(t)
			
			if s_data == "EC" :
				if t_data == "Reaction" :
					adjacency_list[s].append(t)

			if t_data == "Compound" :
				if s_data == "Reaction" :
					adjacency_list[t].append(s)
					cmp_rec_adjacency_list[t].add(s) # Compound Reaction edge
					rec_cmp_adjacency_list[s].add(t) # Reaction Compound edge
					rec_active_list[s] = 1
				  	
def rec_org_edges (vertex_data_map, adjacency_list, rec_org_adjacency_list, rec_active_list, org_rec_adjacency_list) :

	for k, v in adjacency_list.items() :
		if vertex_data_map[k][1] == "Organism" :
			for p in v :
				if  vertex_data_map[p][1] == "EC" and p in adjacency_list :
					for q in adjacency_list[p] :
						if vertex_data_map[q][1] == "Reaction" :
							org_rec_adjacency_list[k].add(q)
							rec_org_adjacency_list[q].add(k) # TODO: don't add if not in rec_cmp_adjacency_list
							rec_active_list[q] = rec_active_list.get(q, 0) + 1

def rec_rec_edges (vertex_data_map, cmp_rec_adjacency_list, rec_active_list, rec_rec_edge_list, max_edge_ID) :

	new_edge_ID = max_edge_ID + 1

	for k, v in cmp_rec_adjacency_list.items() :
		if vertex_data_map[k][1] == "Compound" :
			for p in v :
				if  vertex_data_map[p][1] == "Reaction" and rec_active_list[p] >= 2 :
					for q in v :
						if  vertex_data_map[q][1] == "Reaction" and q != p and q > p and rec_active_list[q] >= 2 :
							if rec_rec_edge_list.get((p, q)) is None :
								rec_rec_edge_list[(p, q)].add(k)
								new_edge_ID = new_edge_ID + 1
								
								rec_active_list[p] = rec_active_list.get(p, 0) + 1
								rec_active_list[q] = rec_active_list.get(q, 0) + 1
							else :
								rec_rec_edge_list[(p, q)].add(k)
							
							# TODO: common Organism		
	
	print(str(new_edge_ID - 1))

def json_output (vertex_data_map, rec_org_adjacency_list, rec_active_list, cmp_rec_adjacency_list, 
	org_rec_adjacency_list, rec_rec_edge_list, max_edge_ID, output_filename) :
	
	new_edge_ID = max_edge_ID + 1
		
	output_file = open(output_filename, "w")

	# node 

	for k, v in rec_active_list.items() :
		if v >= 3 :
			json_obj = json.loads(vertex_data_map[k][0])
			json_obj["Organisms"] = list(map(str, rec_org_adjacency_list[k]))
			json_data = json.dumps(json_obj)
			output_file.write(json_data + "\n")

	for k, v in org_rec_adjacency_list.items() :
		json_obj = json.loads(vertex_data_map[k][0])
		json_data = json.dumps(json_obj)
		output_file.write(json_data + "\n")

	for k, v in cmp_rec_adjacency_list.items() :
		json_obj = json.loads(vertex_data_map[k][0])
		json_data = json.dumps(json_obj)
		output_file.write(json_data + "\n")

	# relationship

	for k, v in rec_rec_edge_list.items() :
		s, t = k
		json_obj = dict()
		json_obj["type"] = "relationship"
		json_obj["id"] = str(new_edge_ID)
		json_obj["label"] = "NULL"
		json_obj["start"] = {"id" : str(s), "labels" : [str(vertex_data_map[s][1])]}
		json_obj["end"] = {"id" : str(t), "labels" : [str(vertex_data_map[t][1])]}
		json_obj["Compounds"] = list(map(str, v))

		json_data = json.dumps(json_obj)
		output_file.write(json_data + "\n")				

		new_edge_ID = new_edge_ID + 1

	output_file.close()	

################################################################################

# main

input_filename_1 = sys.argv[1]
output_filename_1 = sys.argv[2]

vertex_data = dict()
edge_data = dict()
	
with open (input_filename_1, 'r') as read_file : 
	line_count = 0	
	for line in read_file :
		#print(str(line_count) + " : " + line)
		try :
			json_data = json.loads(line)
			parse_json_data_2(json_data, line.strip(), vertex_data, edge_data)
			
		except Exception as e :
			print(e)
			continue					

		line_count = line_count + 1
		#if line_count > 15 :
		#	break

	print("Read " + str(line_count) + " lines")

#for k, v in vertex_data.items() :	
#	print(v[1])

#for k, v in edge_data.items() :
#	print(v[3])

print(len(vertex_data))
print(len(edge_data))
print(maxvid_json_data)
print(maxeid_json_data)

# # #

required_vertex_data_map = {"Reaction" : "13996733750220222916", "Compound" : "1546094233688003844", 
	"Organism" : "11838663228819696432", "EC" : "8121734289560993151"}

adjacency_list = defaultdict(list)

cmp_rec_adjacency_list = defaultdict(set)
rec_cmp_adjacency_list = defaultdict(set)
org_rec_adjacency_list = defaultdict(set)
rec_org_adjacency_list = defaultdict(set)

rec_active_list = dict()

rec_rec_edge_list = defaultdict(set)

# # #

print("filter_edges")

filter_edges(vertex_data, edge_data, required_vertex_data_map, 
	adjacency_list, cmp_rec_adjacency_list, rec_cmp_adjacency_list, rec_active_list)

print(len(adjacency_list))
print(len(cmp_rec_adjacency_list))
print(len(rec_cmp_adjacency_list))
print(len(rec_active_list))

#for k, v in adjacency_list.items() :
#	print(k, end = " ")	
#	print(v)

# # #

print("rec_org_edges")

rec_org_edges(vertex_data, adjacency_list, rec_org_adjacency_list, rec_active_list, org_rec_adjacency_list)

print(len(org_rec_adjacency_list))
print(len(rec_org_adjacency_list))
print(len(rec_active_list))

tmp_count = 0
for k, v in rec_active_list.items() :
	if v >= 2 : # has both Compound and Organism
		tmp_count = tmp_count + 1
print(tmp_count)

# # #

print("rec_rec_edges")

rec_rec_edges(vertex_data, cmp_rec_adjacency_list, rec_active_list, rec_rec_edge_list, maxeid_json_data)

print(len(rec_rec_edge_list))

tmp_count = 0
for k, v in rec_active_list.items() :
	if v >= 3 : # has both Compound and Organism; in the Reaction Reaction graph
		tmp_count = tmp_count + 1
print(tmp_count)		

# # #

print("json_output")

json_output(vertex_data, rec_org_adjacency_list, rec_active_list, cmp_rec_adjacency_list, 
	org_rec_adjacency_list, rec_rec_edge_list, maxeid_json_data, output_filename_1) 

# # #

print("Done.")

sys.exit()

# # # 
