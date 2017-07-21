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

	molecule.delete(trajid)
	return index_to_label



def calc_water_to_residues_map(water_hbonds, solvent_resn):
	"""
	Returns
	-------
	frame_idx: int
		Specify frame index with respect to the smaller trajectory fragment
	water_to_residues: dict mapping string to list of strings
		Map each water molecule to the list of residues it forms
		contacts with (ie {"W:TIP3:8719:OH2:29279" : ["A:PHE:159:N:52441", ...]})
	solvent_bridges: list
		List of hbond interactions between two water molecules
		[("W:TIP3:757:OH2:2312", "W:TIP3:8719:OH2:29279"), ...]
	"""
	water_to_residues = {}
	_solvent_bridges = []
	for frame_idx, atom1_label, atom2_label, itype in water_hbonds:
		if(solvent_resn in atom1_label and solvent_resn in atom2_label): 
			# print atom1_label, atom2_label
			_solvent_bridges.append((atom1_label, atom2_label))
			continue
		elif(solvent_resn in atom1_label and solvent_resn not in atom2_label):
			water = atom1_label
			protein = atom2_label
		elif(solvent_resn not in atom1_label and solvent_resn in atom2_label):
			water = atom2_label
			protein = atom1_label

		if(water not in water_to_residues):
			water_to_residues[water] = set()
		water_to_residues[water].add(protein)

	### Remove duplicate solvent bridges (w1--w2 and w2--w1 are the same)
	solvent_bridges = set()
	for water1, water2 in _solvent_bridges:
		key1 = (water1, water2)
		key2 = (water2, water1)
		if(key1 not in solvent_bridges and key2 not in solvent_bridges):
			solvent_bridges.add(key1)
	solvent_bridges = sorted(list(solvent_bridges))

	return frame_idx, water_to_residues, solvent_bridges