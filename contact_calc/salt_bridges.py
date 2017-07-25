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

__all__ = ['prep_salt_bridge_computation', 'compute_salt_bridges']

##############################################################################
# Functions
##############################################################################

def prep_salt_bridge_computation(traj_frag_molid, frame_idx, sele_id):
	"""
	Compute all possible anion and cation atoms from first frame of simulation

	Returns
	-------
	anion_list: list of strings
		List of atom labels for atoms in ASP or GLU that
		can form salt bridges
	cation_list: list of strings
		List of atom labels for atoms in LYS, ARG, HIS that
		can form salt bridges
	"""
	anion_list = get_anion_atoms(traj_frag_molid, frame_idx, sele_id)
	cation_list = get_cation_atoms(traj_frag_molid, frame_idx, sele_id)
	return anion_list, cation_list


def compute_salt_bridges(traj_frag_molid, frame_idx, sele_id, SALT_BRIDGE_CUTOFF_DISTANCE=4.0):
	"""
	Compute salt bridges in a frame of simulation

	Parameters
	----------
	traj_frag_molid: int
		Identifier to simulation fragment in VMD
	frame_idx: int
		Frame number to query
	sele_id: string, default = None
		Compute contacts on subset of atom selection based on VMD query
	SALT_BRIDGE_CUTOFF_DISTANCE: float, default = 4.0 angstroms
		cutoff for distance between anion and cation atoms 

	Returns
	-------
	salt_bridges: list of tuples, [(frame_index, atom1_label, atom2_label, itype), ...]
		itype = "sb"
	"""
	anion_list, cation_list = prep_salt_bridge_computation(traj_frag_molid, frame_idx, sele_id)
	salt_bridges = []
	for anion_atom in anion_list:
		for cation_atom in cation_list:
			dist = compute_distance(traj_frag_molid, frame_idx, anion_atom, cation_atom)
			if(dist < SALT_BRIDGE_CUTOFF_DISTANCE):
				salt_bridges.append([frame_idx, anion_atom, cation_atom, "sb"])

	return salt_bridges

