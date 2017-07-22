##############################################################################
# MDContactNetworks: A Python Library for computing non-covalent contacts
#                    throughout Molecular Dynamics Trajectories. 
# Copyright 2016-2017 Stanford University and the Authors
#
# Authors: Anthony Kai Kwang Ma
# Email: anthony.ma@yale.edu, anthonyma27@gmail.com, akma327@stanford.edu
##############################################################################

##############################################################################
# Imports
##############################################################################

from vmd import *
import molecule 
import itertools
from contact_utils import *

__all__ = ["compute_vanderwaals"]

##############################################################################
# Globals
##############################################################################

ALPHA_CARBON_DIST_CUTOFF = 10.0 # Angstroms
SOFT_VDW_CUTOFF = 5.0 # Angstroms
VDW_EPSILON = 0.5 # Angstroms
ATOM_RADIUS = {'H': 1.20,
				'C': 1.70,
				'N': 1.55,
				'O': 1.52,
				'S': 1.80}

##############################################################################
# Functions
##############################################################################

def compute_vanderwaals(traj_frag_molid, frame_idx, index_to_label, chain_id):
	"""
	Compute all vanderwaals interactions in a frame of simulation

	Parameters
	----------
	traj_frag_molid: int
		Identifier to simulation fragment in VMD
	frame_idx: int
		Frame number to query
	index_to_label: dict 
		Maps VMD atom index to label "chain:resname:resid:name:index"
		{11205: "A:ASP:114:CA:11205, ...}
	chain_id: string, default = None
		Specify chain of protein to perform computation on 

	Returns
	-------
	vanderwaals: list of tuples, [(frame_idx, atom1_label, atom2_label, itype), ...]
		itype = "vdw"
	"""
	vanderwaals = []
	if(chain_id == None):
		evaltcl("set full_protein [atomselect %s \" noh and protein \" frame %s]" % (traj_frag_molid, frame_idx))
	else:
		evaltcl("set full_protein [atomselect %s \" noh and protein and chain %s \" frame %s]" % (traj_frag_molid, chain_id, frame_idx))

	contacts = evaltcl("measure contacts %s $full_protein" % (SOFT_VDW_CUTOFF))
	contact_index_pairs = parse_contacts(contacts)
	for atom1_index, atom2_index in contact_index_pairs:
		atom1_label, atom2_label = index_to_label[atom1_index], index_to_label[atom2_index]
		element1 = atom1_label.split(":")[3][0]
		element2 = atom2_label.split(":")[3][0]

		distance = compute_distance(traj_frag_molid, frame_idx, atom1_label, atom2_label)
		vanderwaal_cutoff = ATOM_RADIUS[element1] + ATOM_RADIUS[element2] + VDW_EPSILON
		if(distance < vanderwaal_cutoff):
			vanderwaals.append([frame_idx, atom1_label, atom2_label, "vdw"])

	return vanderwaals



# def filter_alpha_carbon(traj_frag_molid, frame_idx, chain_id):
# 	"""
# 	Consider all pairs of alpha carbons throughout the protein and 
# 	compute distance between them. Return list of residue pairs
# 	that are sufficiently close enough for vanderwaals interactions to form

# 	Returns
# 	-------
# 	candidate_residue_pairs: list of tuples
# 		[("ASP:115", "GLU:204"), ...]
# 	"""

# 	if(chain_id == None):
# 		evaltcl("set alpha_carbons [atomselect %s \" name CA \" frame %s]" % (traj_frag_molid, frame_idx))
# 	else:
# 		evaltcl("set alpha_carbons [atomselect %s \" (name CA) and (chain %s) \" frame %s]" % (traj_frag_molid, chain_id, frame_idx))

# 	alpha_carbons = get_atom_selection_labels("alpha_carbons")
# 	evaltcl("$alpha_carbons delete")

# 	candidate_residue_pairs = []
# 	for ca_1, ca_2 in itertools.combinations(alpha_carbons, 2):
# 		distance = compute_distance(traj_frag_molid, frame_idx, ca_1, ca_2)
# 		if(distance < ALPHA_CARBON_DIST_CUTOFF):
# 			res1 = ":".join(ca_1.split(":")[1:3])
# 			res2 = ":".join(ca_2.split(":")[1:3])
# 			candidate_residue_pairs.append((res1, res2))

# 	### USE CONTACTS
# 	return candidate_residue_pairs


# def compute_vanderwaals(traj_frag_molid, frame_idx, index_to_label, chain_id):
# 	"""
# 	Compute all vanderwaals interactions in a frame of simulation

# 	Parameters
# 	----------
# 	traj_frag_molid: int
# 		Identifier to simulation fragment in VMD
# 	frame_idx: int
# 		Frame number to query
# 	index_to_label: dict 
# 		Maps VMD atom index to label "chain:resname:resid:name:index"
# 		{11205: "A:ASP:114:CA:11205, ...}
# 	chain_id: string, default = None
# 		Specify chain of protein to perform computation on 


# 	Returns
# 	-------
# 	vanderwaals: list of tuples, [(frame_idx, atom1_label, atom2_label, itype), ...]
# 		itype = "vdw"
# 	"""
# 	vanderwaals = []

# 	candidate_residue_pairs = filter_alpha_carbon(traj_frag_molid, frame_idx, chain_id)
# 	for res1, res2 in candidate_residue_pairs:
# 		resname1, resid1 = res1.split(":")
# 		resname2, resid2 = res2.split(":")

# 		evaltcl("set res1_atoms [atomselect %s \" resname %s and resid %s \" frame %s]" % (traj_frag_molid, resname1, resid1, frame_idx))
# 		evaltcl("set res2_atoms [atomselect %s \" resname %s and resid %s \" frame %s]" % (traj_frag_molid, resname2, resid2, frame_idx))
# 		contacts = evaltcl("measure contacts %s $res1_atoms $res2_atoms" % (SOFT_VDW_CUTOFF))

# 		contact_index_pairs = parse_contacts(contacts)
# 		for atom1_index, atom2_index in contact_index_pairs:
# 			atom1_label, atom2_label = index_to_label[atom1_index], index_to_label[atom2_index]
# 			distance = compute_distance(traj_frag_molid, frame_idx, atom1_label, atom2_label)
# 			element1 = get_element(traj_frag_molid, frame_idx, atom1_index)
# 			element2 = get_element(traj_frag_molid, frame_idx, atom2_index)
# 			vanderwaal_cutoff = ATOM_RADIUS[element1] + ATOM_RADIUS[element2] + VDW_EPSILON
# 			if(distance < vanderwaal_cutoff):
# 				print(atom1_label, atom2_label, distance)
# 				vanderwaals.append([frame_idx, atom1_label, atom2_label, 'vdw'])

# 		evaltcl("$res1_atoms delete")
# 		evaltcl("$res2_atoms delete")

# 	return vanderwaals





