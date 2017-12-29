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
#import molecule
import itertools
from contact_utils import *

__all__ = ["stratify_hbond_subtypes"]


##############################################################################
# Functions
##############################################################################


def residue_vs_water_hbonds(hbonds, solvent_resn):
    """
    Split hbonds into those involving residues only and those mediated by water.
    """
    residue_hbonds, water_hbonds = [], []
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

        if(atom1 not in backbone_atoms and atom2 not in backbone_atoms):
            hbss.append([frame_idx, atom1_label, atom2_label, "hbss"])

        if((atom1 not in backbone_atoms and atom2 in backbone_atoms) or (atom1 in backbone_atoms and atom2 not in backbone_atoms)):
            hbsb.append([frame_idx, atom1_label, atom2_label, "hbsb"])

        if(atom1 in backbone_atoms and atom2 in backbone_atoms):
            hbbb.append([frame_idx, atom1_label, atom2_label, "hbbb"])

    return hbss, hbsb, hbbb

def stratify_water_bridge(water_hbonds, solvent_resn):
    """
    Infer direct water bridges between residues that both have hbond
    with the same water (ie res1 -- water -- res2)
    """
    frame_idx, water_to_residues, _ = calc_water_to_residues_map(water_hbonds, solvent_resn)
    water_bridges = set()
    ### Infer direct water bridges
    for water in water_to_residues:
        protein_atoms = sorted(list(water_to_residues[water]))
        for res_atom_pair in itertools.combinations(protein_atoms, 2):
            res_atom1, res_atom2 = res_atom_pair
            if(res_atom1 != res_atom2):
                water_bridges.add((frame_idx, res_atom1, water, res_atom2, "wb"))

    wb = sorted([list(entry) for entry in water_bridges])
    return wb


def stratify_extended_water_bridge(water_hbonds, solvent_resn):
    """
    Infer extended water bridges between residues that form hbond with
    water molecules that also have hbond between them.
    (ie res1 -- water1 -- water2 -- res2)
    """
    frame_idx, water_to_residues, solvent_bridges = calc_water_to_residues_map(water_hbonds, solvent_resn)
    extended_water_bridges = set()
    for water1, water2 in solvent_bridges:
        if(water1 not in water_to_residues or water2 not in water_to_residues): continue
        res_atom1_list, res_atom2_list = water_to_residues[water1], water_to_residues[water2]

        for atom1 in res_atom1_list:
            for atom2 in res_atom2_list:
                extended_water_bridges.add((frame_idx, atom1, water1, water2, atom2, "wb2"))

    extended_water_bridges = sorted(list(extended_water_bridges))

    wb2 = []
    for frame_idx, atom1, water1, water2, atom2, itype in extended_water_bridges:
        wb2.append([frame_idx, atom1, water1, water2, atom2, itype])

    return wb2


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
        corresponding to sidechain-sidechain, sidechain-backbone, backbone-backbone,
        water bridge and extended water bridge respectively.
    """
    residue_hbonds, water_hbonds = residue_vs_water_hbonds(hbonds, solvent_resn)
    hbss, hbsb, hbbb = stratify_residue_hbonds(residue_hbonds)
    wb = stratify_water_bridge(water_hbonds, solvent_resn)
    wb2 = stratify_extended_water_bridge(water_hbonds, solvent_resn)
    hbonds = hbss + hbsb + hbbb + wb + wb2

    return hbonds
