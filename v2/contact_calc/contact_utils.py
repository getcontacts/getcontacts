##############################################################################
# MDContactNetworks: A Python Library for computing non-covalent contacts
#					throughout Molecular Dynamics Trajectories. 
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

def get_file_type(file_name):
	"""
	Determine file type by extracting suffix of file_name
	"""
	file_type = file_name.split(".")[-1].strip()
	if(file_type == "nc"): file_type = 'netcdf'
	return file_type


def load_traj(TOP, TRAJ, beg_frame, end_frame, stride):
	"""
	Loads in topology and trajectory into VMD

	Parameters
	----------
	TOP: MD Topology
	TRAJ: MD Trajectory
	beg_frame: int 
	end_frame: int 
	stride: int 

	Returns
	-------
	trajid: int
	simulation molid object 
	"""

	top_file_type = get_file_type(TOP)
	traj_file_type = get_file_type(TRAJ)
	trajid = molecule.load(top_file_type, TOP)
	molecule.read(trajid, traj_file_type, TRAJ, beg=beg_frame, end=end_frame, skip=stride, waitfor=-1)
	return trajid


def gen_index_to_atom_label(TOP, TRAJ):
	"""
	Read in first frame of simulation and generate mapping from 
	VMD index to atom labels

	Parameters
	----------
	TOP: MD Topology
	TRAJ: MD Trajectory

	Returns
	-------
	index_to_label: dict mapping int to string 
		Maps VMD atom index to label "chain:resname:resid:name:index"
		{11205: "A:ASP:114:CA:11205, ...}

	"""
	### Select all atoms from first frame of trajectory
	trajid = load_traj(TOP, TRAJ, 1, 2, 1)
	all_atom_sel = "set all_atoms [atomselect %s \" all \" frame %s]" % (trajid, 0)
	all_atoms = evaltcl(all_atom_sel)

	chains = map(str, evaltcl("$all_atoms get chain").split(" "))
	resnames = map(str, evaltcl("$all_atoms get resname").split(" "))
	resids = map(str, evaltcl("$all_atoms get resid").split(" "))
	names = map(str, evaltcl("$all_atoms get name").split(" "))
	indices = map(str, evaltcl("$all_atoms get index").split(" "))

	### Generate mapping
	index_to_label = {}

	for idx, index in enumerate(indices):
		chain = chains[idx]
		resname = resnames[idx]
		resid = resids[idx]
		name = names[idx]
		atom_label = "%s:%s:%s:%s:%s" % (chain, resname, resid, name, index)
		index_key = int(index)
		index_to_label[index_key] = atom_label

	return index_to_label






	molecule.delete(trajid)
