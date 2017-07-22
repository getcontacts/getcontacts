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
import datetime

__all__ = ['prep_salt_bridge_computation', 'compute_salt_bridges']

##############################################################################
# Globals
##############################################################################
CUTOFF_DISTANCE = 4.0

##############################################################################
# Functions
##############################################################################

def get_anion_atoms(traj_frag_molid, frame_idx, chain_id):
	"""
	Get list of anion atoms that can form salt bridges

	Returns
	-------
	anion_list: list of strings
		List of atom labels for atoms in ASP or GLU that
		can form salt bridges
	"""
	anion_list = []

	if(chain_id == None):
		evaltcl("set ASP [atomselect %s \" (resname ASP) and (name OD1 OD2) \" frame %s]" % (traj_frag_molid, frame_idx))
		evaltcl("set GLU [atomselect %s \" (resname GLU) and (name OE1 OE2) \" frame %s]" % (traj_frag_molid, frame_idx))
	else:
		evaltcl("set ASP [atomselect %s \" (resname ASP) and (name OD1 OD2) and (chain %s) \" frame %s]" % (traj_frag_molid, chain_id, frame_idx))
		evaltcl("set GLU [atomselect %s \" (resname GLU) and (name OE1 OE2) and (chain %s) \" frame %s]" % (traj_frag_molid, chain_id, frame_idx))

	anion_list += get_atom_selection_labels("ASP")
	anion_list += get_atom_selection_labels("GLU")

	evaltcl('$ASP delete')
	evaltcl('$GLU delete')

	return anion_list

def get_cation_atoms(traj_frag_molid, frame_idx, chain_id):
	"""
	Get list of cation atoms that can form salt bridges

	Returns
	-------
	cation_list: list of strings
		List of atom labels for atoms in LYS, ARG, HIS that
		can form salt bridges
	"""
	cation_list = []
	if(chain_id == None):
		evaltcl("set LYS [atomselect %s \" (resname LYS) and (name NZ) \" frame %s]" % (traj_frag_molid, frame_idx))
		evaltcl("set ARG [atomselect %s \" (resname ARG) and (name NH1 NH2) \" frame %s]" % (traj_frag_molid, frame_idx))
		evaltcl("set HIS [atomselect %s \" (resname HIS HSD HSE HSP HIE HIP HID) and (name ND1 NE2) \" frame %s]" % (traj_frag_molid, frame_idx))
	else:
		evaltcl("set LYS [atomselect %s \" (resname LYS) and (name NZ) and (chain %s) \" frame %s]" % (traj_frag_molid, chain_id, frame_idx))
		evaltcl("set ARG [atomselect %s \" (resname ARG) and (name NH1 NH2) and (chain %s) \" frame %s]" % (traj_frag_molid, chain_id, frame_idx))
		evaltcl("set HIS [atomselect %s \" (resname HIS HSD HSE HSP HIE HIP HID) and (name ND1 NE2) and (chain %s) \" frame %s]" % (traj_frag_molid, chain_id, frame_idx))

	cation_list += get_atom_selection_labels("LYS")
	cation_list += get_atom_selection_labels("ARG")
	cation_list += get_atom_selection_labels("HIS")

	evaltcl('$LYS delete')
	evaltcl('$ARG delete')
	evaltcl('$HIS delete')

	return cation_list


def prep_salt_bridge_computation(traj_frag_molid, frame_idx, chain_id):
	"""
	Compute all possible anion and cation atoms from first frame of simulation
	"""
	anion_list = get_anion_atoms(traj_frag_molid, frame_idx, chain_id)
	cation_list = get_cation_atoms(traj_frag_molid, frame_idx, chain_id)
	return anion_list, cation_list


def compute_salt_bridges(traj_frag_molid, frame_idx, index_to_label, anion_list, cation_list):
	salt_bridges = []
	for anion_atom in anion_list:
		for cation_atom in cation_list:
			dist = compute_dist(traj_frag_molid, frame_idx, anion_atom, cation_atom)
			if(dist < CUTOFF_DISTANCE):
				salt_bridges.append([frame_idx, anion_atom, cation_atom, "sb"])
				
	return salt_bridges





