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

__all__ = ["stratify_hbond_subtypes"]


##############################################################################
# Functions
##############################################################################


def residue_vs_water_hbonds(hbonds, solvent_resn):
	"""
	Split hbonds into those involving residues only and those involving water 
	"""
	residue_hbonds = []
	water_hbonds = []
	for hbond in hbonds:
		frame_idx, atom1_label, atom2_label, itype = hbond
		if(solvent_resn in atom1_label or solvent_resn in atom2_label):
			water_hbonds.append(hbond)
		else:
			residue_hbonds.append(hbond)

	return residue_hbonds, water_hbonds

def stratify_residue_hbonds(residue_hbonds):
	"""
	Stratify residue to residue hbonds into those between sidechain-sidechain,
	sidechain-backbone, and backbone-backbone
	"""

	backbone_atoms = ['N', 'O']
	hbss, hbsb, hbbb = [], [], []

	### Iterate through each residue hbond and bin into appropriate subtype
	for frame_idx, atom1_label, atom2_label, itype in residue_hbonds:
		atom1 = atom1_label.split(":")[3]
		atom2 = atom2_label.split(":")[3]

		# hbss
		if(atom1 not in backbone_atoms and atom2 not in backbone_atoms): 
			hbss.append([frame_idx, atom1_label, atom2_label, "hbss"])

		# hbsb
		if((atom1 not in backbone_atoms and atom2 in backbone_atoms) or (atom1 in backbone_atoms and atom2 not in backbone_atoms)):
			hbsb.append([frame_idx, atom1_label, atom2_label, "hbsb"])

		# hbbb
		if(atom1 in backbone_atoms and atom2 in backbone_atoms):
			hbbb.append([frame_idx, atom1_label, atom2_label, "hbbb"])

	return hbss, hbsb, hbbb



def stratify_hbond_subtypes(hbonds, solvent_resn):
	"""
	Stratify the full hbonds list into the following subtypes: sidechain-sidechain,
	sidechain-backbone, backbone-backbone, water-bridge, and extended water-bridge

	Parameters
	----------
	hbonds: list, [[frame_idx, atom1_label, atom2_label, itype], ...]
		List of all hydrogen bond contacts in a single frame. itype = "hb"

	solvent_resn: string, default = TIP3
		Denotes the resname of solvent in simulation

	Returns
	-------
	hbond_subtypes: list, [[frame_idx, atom1_label, atom2_label, itype], ...]
		List of all hydrogen contacts with itype = "hbss", "hbsb", "hbbb", "wb", or "wb2"
	"""

	residue_hbonds, water_hbonds = residue_vs_water_hbonds(hbonds, solvent_resn)
	hbss, hbsb, hbbb = stratify_residue_hbonds(residue_hbonds)


	hbonds = hbss + hbsb + hbbb
	return hbonds
