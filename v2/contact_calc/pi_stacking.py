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
from contact_utils import *

__all__ = ['compute_pi_stacking']

##############################################################################
# Globals
##############################################################################
DISTANCE_CUTOFF = 7.0 # Angstrom
SOFT_DISTANCE_CUTOFF = 10.0 # Angstroms
PI_STACKING_ANGLE_CUTOFF = 30
PSI_ANGLE_CUTOFF = 45

##############################################################################
# Functions
##############################################################################


def get_aromatic_triplet(traj_frag_molid, frame_idx, aromatic_residue_label):
	"""
	Given an aromatic residue label return triplet of atoms on the ring 

	Parameters
	----------
	traj_frag_molid: int
		Identifier to simulation fragment in VMD
	frame_idx: int
		Frame number to query
	aromatic_residue_label: string
		ie "A:PHE:222" 

	Returns
	-------
	aromatic_atom_triplet: list of strings 
		ie ["A:PHE:329:CE1:55228", "A:PHE:329:CE2:55234", "A:PHE:329:CG:55225"]
	"""

	residue_to_atom_names = {"PHE": "CG CE1 CE2", "TRP": "CD2 CZ2 CZ3", "TYR": "CG CE1 CE2"}

	chain, resname, resid = aromatic_residue_label.split(":")
	evaltcl("set aromatic_atoms [atomselect %s \" (chain %s) and (resname %s) and (resid %s) and (name %s)\" frame %s]" %(traj_frag_molid, chain, resname, resid, residue_to_atom_names[resname], frame_idx))
	aromatic_atom_triplet = get_atom_selection_labels("aromatic_atoms")	
	return aromatic_atom_triplet


def compute_pi_stacking(traj_frag_molid, frame_idx, index_to_label, chain_id):
	"""
	Compute pi-stacking interactions in a frame of simulation

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
	pi_stacking = list of tuples, [(frame_index, atom1_label, atom2_label, itype), ...]
		itype = "ps"
	"""
	pi_stacking = []

	if(chain_id == None):
		aromatic_atom_sel = "set aromatic_atoms [atomselect %s \" ((resname PHE) and (name CG CE1 CE2)) or ((resname TRP) and (name CD2 CZ2 CZ3)) or ((resname TYR) and (name CG CE1 CE2)) \" frame %s]" % (traj_frag_molid, frame_idx)
	else:
		aromatic_atom_sel = "set aromatic_atoms [atomselect %s \" ((resname PHE) and (name CG CE1 CE2)) or ((resname TRP) and (name CD2 CZ2 CZ3)) or ((resname TYR) and (name CG CE1 CE2)) and chain %s\" frame %s]" % (traj_frag_molid, chain_id, frame_idx)

	evaltcl(aromatic_atom_sel)
	contacts = evaltcl("measure contacts %s $aromatic_atoms" % (SOFT_DISTANCE_CUTOFF))
	evaltcl("$aromatic_atoms delete")

	### Calculate set of distinct aromatic candidate pairs that may have pi-stacking
	contact_index_pairs = parse_contacts(contacts)
	res_pairs = set()
	residue_to_atom_labels = {}
	for aromatic1_index, aromatic2_index in contact_index_pairs:
		aromatic1_label = index_to_label[aromatic1_index]
		aromatic2_label = index_to_label[aromatic2_index]

		### Check if the two atoms belong to same aromatic group
		aromatic1_res = ":".join(aromatic1_label.split(":")[0:3])
		aromatic2_res = ":".join(aromatic2_label.split(":")[0:3])
		if(aromatic1_res == aromatic2_res): continue

		if(aromatic1_res not in residue_to_atom_labels):
			residue_to_atom_labels[aromatic1_res] = get_aromatic_triplet(traj_frag_molid, frame_idx, aromatic1_res)
		if(aromatic2_res not in residue_to_atom_labels):
			residue_to_atom_labels[aromatic2_res] = get_aromatic_triplet(traj_frag_molid, frame_idx, aromatic2_res)

		k1 = (aromatic1_res, aromatic2_res)
		k2 = (aromatic2_res, aromatic1_res)
		if(k1 not in res_pairs and k2 not in res_pairs):
			res_pairs.add((aromatic1_res, aromatic2_res))

	### Perform strict geometric criterion on candidate aromatic pairs
	for aromatic1_res, aromatic2_res in res_pairs:
		# print(frame_idx, aromatic1_res, aromatic2_res)
		aromatic1_atom_labels = residue_to_atom_labels[aromatic1_res]
		aromatic2_atom_labels = residue_to_atom_labels[aromatic2_res]


		### Distance between two aromatic centers must be below DISTANCE_CUTOFF
		arom1_atom1_coord = get_coord(traj_frag_molid, frame_idx, aromatic1_atom_labels[0])
		arom1_atom2_coord = get_coord(traj_frag_molid, frame_idx, aromatic1_atom_labels[1])
		arom1_atom3_coord = get_coord(traj_frag_molid, frame_idx, aromatic1_atom_labels[2])

		arom2_atom1_coord = get_coord(traj_frag_molid, frame_idx, aromatic2_atom_labels[0])
		arom2_atom2_coord = get_coord(traj_frag_molid, frame_idx, aromatic2_atom_labels[1])
		arom2_atom3_coord = get_coord(traj_frag_molid, frame_idx, aromatic2_atom_labels[2])

		aromatic1_centroid = calc_geom_centroid(arom1_atom1_coord, arom1_atom2_coord, arom1_atom3_coord)
		aromatic2_centroid = calc_geom_centroid(arom2_atom1_coord, arom2_atom2_coord, arom2_atom3_coord)
		aromatic_centers_distance = calc_geom_distance(aromatic1_centroid, aromatic2_centroid)
		if(aromatic_centers_distance > DISTANCE_CUTOFF): continue

		### Angle between vectors normal to each aromatic plane must be below cutoff 
		aromatic1_normal_vector = calc_geom_normal_vector(arom1_atom1_coord, arom1_atom2_coord, arom1_atom3_coord)
		aromatic2_normal_vector = calc_geom_normal_vector(arom2_atom1_coord, arom2_atom2_coord, arom2_atom3_coord)
		aromatic_plane_alignment_angle = calc_angle_between_vectors(aromatic1_normal_vector, aromatic2_normal_vector)
		aromatic_plane_alignment_angle = min(math.fabs(aromatic_plane_alignment_angle - 0), math.fabs(aromatic_plane_alignment_angle - 180))
		if(aromatic_plane_alignment_angle > PI_STACKING_ANGLE_CUTOFF): continue 

		print(frame_idx, aromatic1_res, aromatic2_res)
		### Psi Angle cutoff 
		psi_angle1 = calc_geom_psi_angle(aromatic1_centroid, aromatic2_centroid, aromatic1_normal_vector)
		psi_angle2 = calc_geom_psi_angle(aromatic2_centroid, aromatic1_centroid, aromatic2_normal_vector)
		psi_angle = min(psi_angle1, psi_angle2)
		if(psi_angle > PSI_ANGLE_CUTOFF): continue 

		

		# for arom1_atom_label in aromatic1_atom_labels:
		# 	for arom2_atom_label in aromatic2_atom_labels:
		# 		pi_stacking.append([frame_idx, arom1_atom_label, arom2_atom_label, "ps"])

	return pi_stacking




