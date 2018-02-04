#!/usr/bin/env python3

"""
PyMol visualization of weighted contact network and communication pathways. 

Inputs:
- topology file
- contact frequency file 
- edge list for communication pathway or visualizing portions of full network 

"""

import sys
from pymol import *

def fix_amino_acid_names(key):
	key = key.replace("HSD", "HIS")
	key = key.replace("HSE", "HIS")
	key = key.replace("HSP", "HIS")
	key = key.replace("HIE", "HIS")
	key = key.replace("HIP", "HIS")
	key = key.replace("HID", "HIS")
	key = key.replace("GLH", "GLU")
	key = key.replace("ASH", "ASP")
	key = key.replace("CYP", "CYS")
	key = key.replace("CYX", "CYS")
	return key

def get_edge_to_weight_map(edge_weights):
	"""
	Map the edge between a pair of residues to weight 

	Parameters
	----------
	edge_weights: string
		Path to table associated edge between a pair of amino acids to an edge weight

	Returns
	-------
	edge_weight_map: dictionary
		{("B:SER:60", "B:TYR:99"): 0.64032, ...}
	"""

	edge_weight_map = {}
	f = open(edge_weights, 'r')

	for line in f:
		linfo = line.strip().split("\t")
		aa1, aa2, freq = linfo[0], linfo[1], float(linfo[2])
		aa1 = fix_amino_acid_names(aa1)
		aa2 = fix_amino_acid_names(aa2)
		edge_weight_map[(aa1, aa2)] = freq 

	return edge_weight_map

def gen_atom_selection_str(top_name, node):
	"""
	Generate atom selection string 

	Parameters
	----------
	top_name: string
		Name of the protein structure 
	node: string
		":" delimited atom or residue string

	Return 
	------
	selection_str: string 
		Selection string 
	"""

	node_info = node.split(":")

	# Residue level interactions
	if(len(node_info) == 3):
		chain, resn, resi = node_info
		selection_str = "%s and chain %s and resn %s and resi %s and name CA" % (top_name, chain, resn, resi)

	# Atomic level interactions
	elif(len(node_info) == 4):
		chain, resn, resi, name = node_info
		selection_str = "%s and chain %s and resn %s and resi %s and name %s" % (top_name, chain, resn, resi, name)

	return selection_str

def get_line_thickness(weight):
	"""
	Determine radius of cylinder based on weight of contact edge 

	Parameters
	----------
	weight: float 
		Edge weight

	Return
	------
	radius: float
		Radius of sphere
	"""
	MIN_WIDTH = 0.01
	MAX_WIDTH = 0.5
	radius = float(weight)*MAX_WIDTH
	if(radius < MIN_WIDTH):
		radius = MIN_WIDTH
	return radius

def draw_edge(top_name, node1, node2, weight):
	"""
	Draw an edge between two nodes with specified weight

	Parameters
	----------
	top_name: string
		Name of the protein structure 
	node1: string 
		Represents first residue or atom
	node2: string
		Represents second residue or atom
	weight: float 
		Edge between between node1 and node2
	"""

	label = "%s_contact_network" % (top_name)
	sel1 = gen_atom_selection_str(top_name, node1)
	sel2 = gen_atom_selection_str(top_name, node2)

	### Draw cylinder 
	xyz1 = cmd.get_coords(sel1, 1)[0]
	xyz2 = cmd.get_coords(sel2, 1)[0]
	x1, y1, z1 = xyz1
	x2, y2, z2 = xyz2
	r1, g1, b1 = 0.5, 0.5, 0.5 # color (red)
	r2, g2, b2 = 0.5, 0.5, 0.5 # color (red)
	radius = get_line_thickness(weight)
	cmd.load_cgo( [ 9.0, x1, y1, z1, x2, y2, z2, radius, r1, g1, b1, r2, g2, b2 ], "cylinder_%s_%s" % (node1, node2))

	cmd.do("show spheres, %s" % (sel1))
	cmd.do("show spheres, %s" % (sel2))

def visualize_protein_network(structure, edge_weights, sub_network=[], cutoff=0.0):
	"""
	PyMol rendering of protein structure network with edge weights derived from MD simulation

	Parameters
	----------
	structure: string
		Path to structure to superimpose network on
	edge_weights: string
		Path to table associated edge between a pair of amino acids to an edge weight
	sub_network: list of tuples 
		Default is empty array. User can use this argument to draw a subnetwork or shortest path
		between a pair of nodes
	cutoff: float
		Default 0.0, cutoff value for frequency of contact 
	"""

	### Load in background protein structure
	top_name = top.split("/")[-1].split(".")[0]
	cmd.bg_color("white")
	cmd.load(top)
	cmd.hide()
	cmd.show("cartoon")
	cmd.cartoon("loop")
	cmd.color("white", top_name)
	cmd.do("set cartoon_transparency, 0.7")
	cmd.do("set sphere_scale, 0.35")
	cmd.do("set sphere_transparency, 0.5")
	cmd.do("set sphere_color, deepteal")


	### Draw protein network with edge thickness representing the frequency of a contact in simulation
	edge_weight_map = get_edge_to_weight_map(edge_weights)
	
	if(sub_network == []): # If sub_network argument empty then draw entire network
		for edge in edge_weight_map:
			aa1, aa2 = edge 
			weight = edge_weight_map[edge]
			if(aa1 == aa2): continue
			if(weight < cutoff): continue 
			draw_edge(top_name, aa1, aa2, weight)
	else:
		for aa1, aa2 in sub_network:
			aa1 = fix_amino_acid_names(aa1)
			aa2 = fix_amino_acid_names(aa2)
			k1, k2 = (aa1, aa2), (aa2, aa1)
			if(k1 in edge_weight_map):
				weight = edge_weight_map[k1]
			elif(k2 in edge_weight_map):
				weight = edge_weight_map[k2]
			if(aa1 == aa2): continue
			if(weight < cutoff): continue
			draw_edge(top_name, aa1, aa2, weight)


# Example 
top = "top.pdb"
contact_freq="wb_freq.tsv"
sub_network = [('C:HSD:62', 'C:PRO:1'), ('A:TYR:95', 'B:PRO:1'), ('A:TYR:99', 'A:HSD:62'), ('C:TYR:99', 'C:TYR:95'), ('B:TYR:99', 'A:HSD:62'), ('B:TYR:99', 'B:TYR:95'), ('A:TYR:99', 'A:TYR:95'), ('C:TYR:99', 'C:GLY:65'), ('C:GLY:65', 'C:PRO:1'), ('B:TYR:99', 'B:SER:63'), ('C:TYR:99', 'A:ASN:97'), ('A:ASN:97', 'B:PRO:1'), ('C:TYR:95', 'A:PRO:1'), ('B:SER:63', 'B:PRO:1'), ('A:TYR:99', 'C:HSD:62'), ('B:TYR:95', 'C:PRO:1'), ('A:HSD:62', 'A:PRO:1')]
visualize_protein_network(top, contact_freq, sub_network)


