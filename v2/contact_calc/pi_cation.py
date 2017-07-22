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

__all__ = ['prep_pi_cation_computation', 'compute_pi_cation']

##############################################################################
# Globals
##############################################################################
DISTANCE_CUTOFF = 60.0 # Angstrom
ANGLE_CUTOFF = 60 # Degree


##############################################################################
# Functions
##############################################################################

def prep_pi_cation_computation(traj_frag_molid, frame_idx, chain_id):
	"""
	Compute all possible cation atom and aromatic residues in simulation

	Returns
	-------
	cation_list: list of strings
		List of atom labels for atoms in LYS, ARG, HIS that
		can form salt bridges
	aromatic_atom_triplet_list: list of strings
		List of TYR, TRP, PHE residues that can form pi-cation interactions.
		strings formatted as "chain:resname:resid" 
	"""

	cation_list = get_cation_atoms(traj_frag_molid, frame_idx, chain_id)
	aromatic_atom_triplet_list = get_aromatic_atom_triplets(traj_frag_molid, frame_idx, chain_id)

	return cation_list, aromatic_atom_triplet_list

def compute_pi_cation(traj_frag_molid, frame_idx, cation_list, aromatic_atom_triplet_list):
	"""
	Compute pi-cation interactions in a frame of simulation

	Parameters
	----------
	traj_frag_molid: int
		Identifier to simulation fragment in VMD
	frame_idx: int
		Frame number to query
	cation_list: list of strings
		List of atom labels for atoms in LYS, ARG, HIS that
		can form salt bridges
	aromatic_atom_triplet_list: list of strings
		List of TYR, TRP, PHE residues that can form pi-cation interactions.
		strings formatted as "chain:resname:resid" 

	Returns
	-------
	pi_cations = list of tuples, [(frame_index, atom1_label, atom2_label, itype), ...]
		itype = "pc"
	"""

	pi_cations = []

	for cidx, cation_atom_label in enumerate(cation_list):
		cation_coord = get_coord(traj_frag_molid, frame_idx, cation_atom_label)
		for arom_atom1_label, arom_atom2_label, arom_atom3_label in aromatic_atom_triplet_list:
			arom_atom1_coord = get_coord(traj_frag_molid, frame_idx, arom_atom1_label)
			arom_atom2_coord = get_coord(traj_frag_molid, frame_idx, arom_atom2_label)
			arom_atom3_coord = get_coord(traj_frag_molid, frame_idx, arom_atom3_label)

			### Perform distance criterion
			aromatic_centroid = calc_geom_centroid(arom_atom1_coord, arom_atom2_coord, arom_atom3_coord)
			cation_to_centroid_distance = calc_geom_distance(cation_coord, aromatic_centroid)
			if(cation_to_centroid_distance > DISTANCE_CUTOFF): continue 

			### Perform angle criterion
			aromatic_plane_norm_vec = calc_geom_normal_vector(arom_atom1_coord, arom_atom2_coord, arom_atom3_coord)
			aromatic_center_to_cation_vec = points_to_vector(aromatic_centroid, cation_coord)
			cation_norm_offset_angle = calc_angle_between_vectors(aromatic_plane_norm_vec, aromatic_center_to_cation_vec)
			cation_norm_offset_angle = min(math.fabs(cation_norm_offset_angle - 0), math.fabs(cation_norm_offset_angle - 180))
			if(cation_norm_offset_angle > ANGLE_CUTOFF): continue

			print(cation_to_centroid_distance, cation_norm_offset_angle)

			### Append three of the aromatic atoms
			pi_cations.append([frame_idx, cation_atom_label, arom_atom1_label, "pc"])
			pi_cations.append([frame_idx, cation_atom_label, arom_atom2_label, "pc"])
			pi_cations.append([frame_idx, cation_atom_label, arom_atom2_label, "pc"])


	return pi_cations


