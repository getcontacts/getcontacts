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
import time
from stratify_hbonds import *
from stratify_ligand_hbonds import *

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

def calc_donor_acceptor_pairs(traj_frag_molid, frame_idx, solvent_resn, sele_id, ligand, HBOND_CUTOFF_DISTANCE, HBOND_CUTOFF_ANGLE):
    """
    Compute donor and acceptor atom pairs for hydrogen bonds in terms of numeric VMD indices
    """
    ### Measure Hbonds command
    if(ligand == None):
        if(sele_id == None):
            measure_hbonds_command = "measure hbonds %s %s [atomselect %s \"(resname %s and within %s of protein) or protein and not lipid and not carbon and not sulfur\" frame %s]" % (HBOND_CUTOFF_DISTANCE, HBOND_CUTOFF_ANGLE, traj_frag_molid, solvent_resn, WATER_TO_PROTEIN_DIST, frame_idx)
        else:
            measure_hbonds_command = "measure hbonds %s %s [atomselect %s \"(resname %s and within %s of (protein and (%s))) or protein and (%s) and not lipid and not carbon and not sulfur\" frame %s]" % (HBOND_CUTOFF_DISTANCE, HBOND_CUTOFF_ANGLE, traj_frag_molid, solvent_resn, WATER_TO_PROTEIN_DIST, sele_id, sele_id, frame_idx)
    else:
        if(sele_id == None):
            measure_hbonds_command = "measure hbonds %s %s [atomselect %s \"(resname %s and within %s of resname %s) or (not carbon and not sulfur and protein within %s of resname %s) or (not carbon and not sulfur and resname %s) and (not lipid)\" frame %s]" % (HBOND_CUTOFF_DISTANCE, HBOND_CUTOFF_ANGLE, traj_frag_molid, solvent_resn, WATER_TO_LIGAND_DIST, ligand, WATER_TO_LIGAND_DIST, ligand, ligand, frame_idx)
        else:
            measure_hbonds_command = "measure hbonds %s %s [atomselect %s \"(resname %s and within %s of resname %s) or ((not carbon and not sulfur and protein and (%s)) and within %s of resname %s) or (not carbon and not sulfur and resname %s) and (not lipid)\" frame %s]" % (HBOND_CUTOFF_DISTANCE, HBOND_CUTOFF_ANGLE, traj_frag_molid, solvent_resn, WATER_TO_LIGAND_DIST, ligand, sele_id, WATER_TO_LIGAND_DIST, ligand, ligand, frame_idx)

    donor_acceptor_indices = evaltcl(measure_hbonds_command)

    ### Parse atom indices
    donor_acceptor_lists = donor_acceptor_indices.split("}")
    donor_list = donor_acceptor_lists[0].split("{")[1].split(" ")
    acceptor_list = donor_acceptor_lists[1].split("{")[1].split(" ")

    donors, acceptors = [], []
    for idx, d in enumerate(donor_list):
        a = acceptor_list[idx]
        if(d == "" or a == ""):continue
        donors.append(int(d))
        acceptors.append(int(a))

    return donors, acceptors

def compute_hydrogen_bonds(traj_frag_molid, frame_idx, index_to_label, solvent_resn, sele_id, ligand=None, HBOND_CUTOFF_DISTANCE=3.5, HBOND_CUTOFF_ANGLE=70):
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
    sele_id: string, default = None
        Compute contacts on subset of atom selection based on VMD query
    HBOND_CUTOFF_DISTANCE: float, default = 3.5 Angstroms
    HBOND_CUTOFF_ANGLE: float, default = 70 degrees

    Return
    ------
    hbonds: list of tuples, [(frame_idx, atom1_label, atom2_label, itype), ...]
    """
    itype = "hb"
    if(ligand != None): itype = "lhb"

    hbonds = []
    donors, acceptors = calc_donor_acceptor_pairs(traj_frag_molid, frame_idx, solvent_resn, sele_id, ligand, HBOND_CUTOFF_DISTANCE, HBOND_CUTOFF_ANGLE)
    donors, acceptors = filter_duplicates(donors, acceptors)

    for idx, donor in enumerate(donors):
        acceptor = acceptors[idx]
        donor_label, acceptor_label = index_to_label[donor], index_to_label[acceptor]
        hbonds.append([frame_idx, donor_label, acceptor_label, itype])


    ### Perform post processing on hbonds list to stratify into different subtypes
    if(itype == "hb"):
        hbond_subtypes = stratify_hbond_subtypes(hbonds, solvent_resn)
    elif(itype == "lhb"):
        hbond_subtypes = stratify_ligand_hbond_subtypes(hbonds, solvent_resn, ligand)

    return hbond_subtypes
