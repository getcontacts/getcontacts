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
import time

__all__ = ['compute_hydrogen_bonds']

##############################################################################
# Globals
##############################################################################
WATER_TO_PROTEIN_DIST = 5
WATER_TO_LIGAND_DIST = 12


##############################################################################
# Functions
##############################################################################


def filter_duplicates(donors, acceptors):
	"""
	Filter out duplicate donor acceptor atom pairs
	"""

	pairs = sorted(list(set([(d, acceptors[idx]) for idx, d in enumerate(donors)])))

	new_donors, new_acceptors = [], []
	for d, a in pairs:
		new_donors.append(d)
		new_acceptors.append(a)

	return new_donors, new_acceptors


def calc_donor_acceptor_pairs(traj_frag_molid, frame_idx, solvent_resn, chain_id, ligand, distance_cutoff, angle_cutoff):
	"""
	Compute donor and acceptor atom pairs for hydrogen bonds in terms of numeric VMD indices
	"""

	### Measure Hbonds command

	if(ligand == None):
		if(chain_id == None):
			measure_hbonds_command = "measure hbonds %s %s [atomselect %s \"(resname %s and within %s of protein) or protein and not lipid and not carbon and not sulfur\" frame %s]" % (distance_cutoff, angle_cutoff, traj_frag_molid, solvent_resn, WATER_TO_PROTEIN_DIST, frame_idx)
		else:
			measure_hbonds_command = "measure hbonds %s %s [atomselect %s \"(resname %s and within %s of protein and chain %s) or protein and chain %s and not lipid and not carbon and not sulfur\" frame %s]" % (distance_cutoff, angle_cutoff, traj_frag_molid, solvent_resn, WATER_TO_PROTEIN_DIST, chain_id, chain_id, frame_idx)
	else:
		if(chain_id == None):
			measure_hbonds_command = "measure hbonds %s %s [atomselect %s \"(resname %s and within %s of resname %s) or (not carbon and not sulfur and protein within %s of resname %s) or (not carbon and not sulfur and resname %s) and (not lipid)\" frame %s]" % (distance_cutoff, angle_cutoff, traj_frag_molid, solvent_resn, WATER_TO_LIGAND_DIST, ligand, WATER_TO_LIGAND_DIST, ligand, ligand, frame_idx)
		else:
			measure_hbonds_command = "measure hbonds %s %s [atomselect %s \"(resname %s and within %s of resname %s) or (not carbon and not sulfur and protein and chain %s and within %s of resname %s) or (not carbon and not sulfur and resname %s) and (not lipid)\" frame %s]" % (distance_cutoff, angle_cutoff, traj_frag_molid, solvent_resn, WATER_TO_LIGAND_DIST, ligand, chain_id, WATER_TO_LIGAND_DIST, ligand, ligand, frame_idx)


	donor_acceptor_indices = evaltcl(measure_hbonds_command)

	### Parse atom indices
	donor_acceptor_lists = donor_acceptor_indices.split("}")
	donor_list = donor_acceptor_lists[0].split("{")[1].split(" ")
	donors = [int(d) for d in donor_list]

	acceptor_list = donor_acceptor_lists[1].split("{")[1].split(" ")
	acceptors = [int(a) for a in acceptor_list]

	return donors, acceptors


def compute_hydrogen_bonds(traj_frag_molid, frame_idx, index_to_label, solvent_resn, chain_id, ligand=None, distance_cutoff=3.5, angle_cutoff=70):
	"""
	Compute hydrogen bonds involving protein for a single frame of simulation

	Parameters
	----------
	traj_frag_molid: int
		Specifies which trajectory fragment in VMD to perform computations upon
	frame_idx: int
		Specify frame index with respect to the smaller trajectory fragment
	solvent_resn: string, default = TIP3
		Denotes the resname of solvent in simulation
	chain_id: string, default = None
		Specify chain of protein to perform computation on 
	distance_cutoff: float, default = 3.5 Angstroms
	angle_cutoff: float, default = 70 degrees

	Return
	------
	hbonds: list of tuples, [(frame_idx, atom1_label, atom2_label, itype), ...]
	"""
	itype = "hb"
	if(ligand != None): itype = "lhb"

	hbonds = []
	donors, acceptors = calc_donor_acceptor_pairs(traj_frag_molid, frame_idx, solvent_resn, chain_id, ligand, distance_cutoff, angle_cutoff)
	donors, acceptors = filter_duplicates(donors, acceptors)

	for idx, donor in enumerate(donors):
		acceptor = acceptors[idx]
		donor_label, acceptor_label = index_to_label[donor], index_to_label[acceptor]

		### If computing ligand contacts then interaction must involve ligand molecule
		if(itype == "lhb" and ligand not in donor_label and ligand not in acceptor_label): continue
		hbonds.append([frame_idx, donor_label, acceptor_label, itype])

	return hbonds

