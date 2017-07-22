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

	for cidx, cation_atom in enumerate(cation_list):
		if(cidx > 0): break
		print(cation_atom, get_coord(traj_frag_molid, frame_idx, cation_atom))


	return pi_cations



