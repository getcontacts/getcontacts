############################################################################
# Copyright 2018 Anthony Ma & Stanford University                          #
#                                                                          #
# Licensed under the Apache License, Version 2.0 (the "License");          #
# you may not use this file except in compliance with the License.         #
# You may obtain a copy of the License at                                  #
#                                                                          #
#     http://www.apache.org/licenses/LICENSE-2.0                           #
#                                                                          #
# Unless required by applicable law or agreed to in writing, software      #
# distributed under the License is distributed on an "AS IS" BASIS,        #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. #
# See the License for the specific language governing permissions and      #
# limitations under the License.                                           #
############################################################################

############################################################################
# Imports
############################################################################
from __future__ import print_function

from vmd import *
import numpy as np
import math
import re
import sys
import os
from contextlib import contextmanager
from numpy.linalg import norm
from .atom import Atom

int_pattern = re.compile(r'\d+')
float_pattern = re.compile(r'-?\d+\.\d+')


def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    return [atoi(c) for c in re.split(r'(\d+)', text)]


def get_file_type(file_name):
    """
    Determine file type by extracting suffix of file_name
    """
    if file_name is None:
        return None

    file_type = file_name.split(".")[-1].strip()
    if file_type == "nc":
        file_type = "netcdf"
    if file_type == "prmtop":
        file_type = "parm7"
    if file_type == "cms":
        file_type = "mae"
    if file_type == "cif":
        file_type = "pdbx"

    return file_type


@contextmanager
def suppress_stdout():
    """
    Temporarily suppresses stdout.

    From: https://stackoverflow.com/questions/4178614/suppressing-output-of-module-calling-outside-library

    Example
    =======
        with suppress_stdout():
            print "You cannot see this"
    """
    with open('/dev/null', "w") as devnull:
        old_stdout = os.dup(sys.stdout.fileno())
        os.dup2(devnull.fileno(), 1)
        try:
            yield
        finally:
            os.dup2(old_stdout, 1)


def load_traj(top, traj, beg_frame, end_frame, stride):
    """
    Loads in topology and trajectory into VMD

    Parameters
    ----------
    top: MD Topology
    traj: MD Trajectory
    beg_frame: int
    end_frame: int
    stride: int

    Returns
    -------
    trajid: int
        simulation molid object
    """
    with suppress_stdout():
        top_file_type = get_file_type(top)
        traj_file_type = get_file_type(traj)
        trajid = molecule.load(top_file_type, top)
        molecule.delframe(trajid)  # Ensure topology doesn't count as a frame
        if traj is not None:
            molecule.read(trajid, traj_file_type, traj, beg=beg_frame, end=end_frame, skip=stride, waitfor=-1)
        else:
            molecule.read(trajid, top_file_type, top, beg=beg_frame, end=end_frame, skip=stride, waitfor=-1)

    return trajid


def simulation_length(top, traj):
    """
    Computes the simulation length efficiently
    """

    trajid = load_traj(top, traj, 0, -1, 100)
    num_frags = molecule.numframes(trajid)

    # There are between (num_frags-1)*100 and num_frags*100 frames.
    # Read all frames of the last fragment to determine the exact amount
    trajid = load_traj(top, traj, (num_frags-1)*100, -1, 1)
    last_frag_frames = molecule.numframes(trajid)
    return (num_frags - 1) * 100 + last_frag_frames


# def get_atom_selection_labels(selection_id):
#     """
#     Returns set of atom labels for each atom in selection_id
#     """
#     chains, resnames, resids, names, indices = get_atom_selection_properties(selection_id)
#     atom_labels = set()
#     for idx in range(len(chains)):
#         chain, resname, resid, name, index = chains[idx], resnames[idx], resids[idx], names[idx], indices[idx]
#         label = "%s:%s:%s:%s:%s" % (chain, resname, resid, name, index)
#         atom_labels.add(label)

#     return atom_labels

def get_atom_selection_indices(selection_id):
    """
    *** Major optimization to be implemented ***
    GetContacts should perform all VMD based computation with just atom indices until needing to convert to labels with index_to_label
    """

    return set(map(int, int_pattern.findall(evaltcl("$" + selection_id + " get index"))))

    # Old old
    # chains, resnames, resids, names, elements, indices = get_atom_selection_properties(selection_id)

    # Old
    # indices = list(map(str, evaltcl("$%s get index" % selection_id).split(" ")))
    # if indices == ['']:
    #     indices = []
    #
    # return set([int(idx) for idx in indices])

    # # Old old
    # atom_indices = set()
    # for idx in range(len(indices)):
    #     index = int(indices[idx])
    #     atom_indices.add(index)
    #
    # return atom_indices
    

def get_atom_selection_properties(selection_id):
    """
    After executing an evaltcl atom selection command, this function
    is called to retrieve the chain, resname, resid, name, and index
    of each atom in the selection

    Parameters
    ----------
    selection_id: string
        Denotes the atom selection identifier

    Returns
    -------
    chains: list of strings
        Chain identifier for each atom ["A", "B", "A", ...]
    resnames: list of strings
        Resname for each atom ["ASP", "GLU", "TYR", ...]
    resids: list of strings
        Residue index for each atom ["12", "45", ...]
    names: list of strings
        Atom name of each atom ["NZ", "CA", "N", ...]
    indices: list of strings
        VMD based index for each atom ["12304", "1231", ...]
    icodes: list of strings
        insertion code for each atom where present, and '' otherwise
        ['','','','','A','A','A','','','',...]
        """
    chains   = safely_parsed_evaltcl("$%s get chain" % selection_id)
    resnames = safely_parsed_evaltcl("$%s get resname" % selection_id)
    resids   = safely_parsed_evaltcl("$%s get resid" % selection_id)
    names    = safely_parsed_evaltcl("$%s get name" % selection_id)
    elements = safely_parsed_evaltcl("$%s get element" % selection_id)
    indices  = safely_parsed_evaltcl("$%s get index" % selection_id)
    icodes   = safely_parsed_evaltcl("$%s get insertion" % selection_id)
    if icodes == [''] or chains == [''] or resnames == [''] or resids == [''] or names == [''] or indices == ['']:
        return [], [], [], [], [], []
    return chains, resnames, resids, names, elements, indices, icodes

def safely_parsed_evaltcl(tclstr):
    ''' When a value is not found, VMD outputs '{ }', which can throw off a split function.
    Therefore we define a function here which takes this into account and correctly returns values 
    from evaltcl. '''
    nonchar = '{ }'
    working_string = evaltcl(tclstr).replace(nonchar, '{}')
    working_list = working_string.split(' ')
    working_list = [str(i) if i != '{}' else '' for i in working_list]
    return working_list

def gen_index_to_atom(top, traj):
    """
    Read in first frame of simulation and generate mapping from VMD index to atom labels

    Parameters
    ----------
    top: MD Topology
    traj: MD Trajectory

    Returns
    -------
    index_to_atom: dict
        Maps VMD atom index to Atom

    """
    # Select all atoms from first frame of trajectory
    trajid = load_traj(top, traj, 0, 1, 1)
    all_atom_sel = "set all_atoms [atomselect %s \" all \" frame %s]" % (trajid, 0)
    evaltcl(all_atom_sel)
    chains, resnames, resids, names, elements, indices, icodes = get_atom_selection_properties("all_atoms")
    evaltcl('$all_atoms delete')
    # Generate mapping
    index_to_atom = {}
    for idx, index in enumerate(indices):
        chain = chains[idx]
        resname = resnames[idx]
        resid = resids[idx]
        name = names[idx]
        element = elements[idx]
        icode = icodes[idx]
        # atom_label = "%s:%s:%s:%s:%s" % (chain, resname, resid, name, index)
        index_key = int(index)
        # index_to_atom[index_key] = atom_label
        index_to_atom[index_key] = Atom(int(index_key), chain, resname, int(resid), name, element, icode=icode)

    molecule.delete(trajid)
    return index_to_atom


# def get_selection_atoms(traj_frag_molid, frame_idx, selection_id):
#     """
#     Get list of atoms in protein within VMD selection query 

#     Returns
#     -------
#     protein_selection_list: list of strings 
#         List of atom labels for atoms in protein selection 
#     """
#     # evaltcl("set selection_id [atomselect %s \" (protein) and (%s) \" frame %s]" % (traj_frag_molid, selection_id, frame_idx))
#     evaltcl("set selection_id [atomselect %s \" (%s) \" frame %s]" % (traj_frag_molid, selection_id, frame_idx))
#     protein_selection_set = get_atom_selection_labels("selection_id")
#     evaltcl('$selection_id delete')
#     return protein_selection_set

def get_selection_indices(traj_frag_molid, frame_idx, selection_id):
    """
    Get list of atoms in protein within VMD selection query 

    Returns
    -------
    protein_selection_indices: set of int
        Set of atom indices for atoms in protein selection 
    """
    # evaltcl("set selection_id [atomselect %s \" (protein) and (%s) \" frame %s]" % (traj_frag_molid, selection_id, frame_idx))
    evaltcl("set selection_id [atomselect %s \" (%s) \" frame %s]" % (traj_frag_molid, selection_id, frame_idx))
    protein_selection_indices = get_atom_selection_indices("selection_id")
    evaltcl('$selection_id delete')
    return protein_selection_indices


# def get_aromatic_atom_triplets(traj_frag_molid, frame_idx, chain_id):
#     """
#     Get list of aromatic atom triplets

#     Returns
#     -------
#     aromatic_atom_triplet_list: list of tuples corresponding to three equally spaced points
#     on the 6-membered rings of TYR, TRP, or PHE residues.
#         [(A:PHE:72:CG:51049, A:PHE:72:CE1:51052, A:PHE:72:CE2:51058), ...]
#     """
#     aromatic_atom_list = []
#     if chain_id is None:
#         evaltcl("set PHE [atomselect %s \" (resname PHE) and (name CG CE1 CE2) \" frame %s]" %
#                 (traj_frag_molid, frame_idx))
#         evaltcl("set TRP [atomselect %s \" (resname TRP) and (name CD2 CZ2 CZ3) \" frame %s]" %
#                 (traj_frag_molid, frame_idx))
#         evaltcl("set TYR [atomselect %s \" (resname TYR) and (name CG CE1 CE2) \" frame %s]" %
#                 (traj_frag_molid, frame_idx))
#     else:
#         evaltcl("set PHE [atomselect %s \" (resname PHE) and (name CG CE1 CE2) and (chain %s)\" frame %s]" %
#                 (traj_frag_molid, frame_idx, chain_id))
#         evaltcl("set TRP [atomselect %s \" (resname TRP) and (name CD2 CZ2 CZ3) and (chain %s)\" frame %s]" %
#                 (traj_frag_molid, frame_idx, chain_id))
#         evaltcl("set TYR [atomselect %s \" (resname TYR) and (name CG CE1 CE2) and (chain %s)\" frame %s]" %
#                 (traj_frag_molid, frame_idx, chain_id))

#     aromatic_atom_list += get_atom_selection_labels("PHE")
#     aromatic_atom_list += get_atom_selection_labels("TRP")
#     aromatic_atom_list += get_atom_selection_labels("TYR")

#     evaltcl("$PHE delete")
#     evaltcl("$TRP delete")
#     evaltcl("$TYR delete")

#     aromatic_atom_triplet_list = []

#     # Generate triplets of the three equidistant atoms on an aromatic ring
#     for i in range(0, len(aromatic_atom_list), 3):
#         aromatic_atom_triplet_list.append(aromatic_atom_list[i:i+3])

#     return aromatic_atom_triplet_list


def convert_to_single_atom_aromatic_string(aromatic_atom_labels):
    """
    Checks the triple of aromatic atom labels and returns the one with lower atom name

    Parameters
    ----------
    aromatic_atom_labels: list of str
        triple of atom labels
    Returns
    -------
    str: A label of a single atom representative of the aromatic ring (ie C:PHE:49:CG)
    """
    return min(aromatic_atom_labels, key=lambda l: l.split(":")[3])
    # aromatic_CG_atom_label = ":".join(aromatic_atom_label.split(":")[0:3]) + ":CG:vmd_idx"
    # return aromatic_CG_atom_label


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
    frame_idx = 0
    water_to_residues = {}
    _solvent_bridges = []
    for frame_idx, atom1_label, atom2_label, itype in water_hbonds:
        if solvent_resn in atom1_label and solvent_resn in atom2_label:
            _solvent_bridges.append((atom1_label, atom2_label))
            continue
        elif solvent_resn in atom1_label and solvent_resn not in atom2_label:
            water = atom1_label
            protein = atom2_label
        elif solvent_resn not in atom1_label and solvent_resn in atom2_label:
            water = atom2_label
            protein = atom1_label
        else:
            raise ValueError("Solvent residue name can't be resolved")

        if water not in water_to_residues:
            water_to_residues[water] = set()
        water_to_residues[water].add(protein)

    # Remove duplicate solvent bridges (w1--w2 and w2--w1 are the same)
    solvent_bridges = set()
    for water1, water2 in _solvent_bridges:
        key1 = (water1, water2)
        key2 = (water2, water1)
        if key1 not in solvent_bridges and key2 not in solvent_bridges:
            solvent_bridges.add(key1)
    solvent_bridges = sorted(list(solvent_bridges))

    return frame_idx, water_to_residues, solvent_bridges

def find_disulfide(top, traj):
    """
    Find Cysteines making disulfide bridges and return their resids
    """
    disulfide_pairs = list()
    molid = load_traj(top, traj, 0, 1, 1)
    evaltcl("set cys_all [atomselect %s \" resname CYS  \" frame 0]" % molid)
    cys_all = set(evaltcl("$cys_all get resid").split())
    for cys in cys_all:
        evaltcl("set bridge_cys [atomselect %s \" name SG and within 2.1 of (resid %s and name SG)  \" frame 0]" % (molid, cys))
        disulfide_pair = set(map(int,evaltcl("$bridge_cys get resid").split()))
        if len(disulfide_pair) == 2 and disulfide_pair not in disulfide_pairs:
            disulfide_pairs.append(disulfide_pair)
    return(disulfide_pairs)

def configure_solv(top, traj, solvent_sele):
    """
    Detects the solvent atoms and creates a corresponding VMD selection macro called 'solv'. Will print a warning
    if no solvent molecule could be located and returns an empty set.

    Parameters
    ----------
    top: Topology
        In .pdb or .mae format
    traj: Trajectory
        In .nc or .dcd format
    solvent_sele: str
        The solvent selection. If empty or None, attempt to determine solvent resname

    Returns
    -------
    set of int:
        A set of vmd-indices of atoms that are in the solvent selection
    """
    # Set up custom water selection
    if solvent_sele:
        evaltcl("atomselect macro solv \" " + solvent_sele + " \"")

        # Determine resnames in selection
        molid = load_traj(top, traj, 0, 1, 1)
        evaltcl("set all_solv [atomselect %s \" solv \" frame 0]" % molid)
        solv_resn = set(evaltcl("$all_solv get resname").split())
        solv_ids = get_atom_selection_indices("all_solv")
        evaltcl("$all_solv delete")
        molecule.delete(molid)
        if solv_ids:
            print("Detected %d solvent atoms (resnames %s) matching --solv '%s'" %
                  (len(solv_ids), " ".join(set(solv_resn)), solvent_sele))
        else:
            print("Detected no solvent atoms matching --solv '%s'" % solvent_sele)
        return solv_ids
    else:
        solv_resnames = set("H2O HH0 OHH HOH OH2 SOL WAT TIP TIP2 TIP3 TIP4 T3P IP3".split())
        molid = load_traj(top, traj, 0, 1, 1)
        evaltcl("set all_atoms [atomselect %s \" all \" frame 0]" % molid)
        all_resn = set(evaltcl("$all_atoms get resname").split())
        solv_resnames = solv_resnames & all_resn
        evaltcl("$all_atoms delete")

        if not solv_resnames:
            print("Detected no solvent atoms (specify manually using --solv)")
            evaltcl("atomselect macro solv \" none \"")
            molecule.delete(molid)
            return set([])
        else:
            solvent_resn = " ".join(solv_resnames)
            evaltcl("atomselect macro solv \" (resname " + solvent_resn + ") \"")
            solv_ids = get_selection_indices(molid, 0, "solv")
            print("Detected %d solvent atoms (resname %s)" % (len(solv_ids), solvent_resn))
            molecule.delete(molid)
            return solv_ids

def configure_lipid(top, traj, lipid_sele):
    """
    Detects the lipid residue name and creates a corresponding VMD selection macro called 'lipid'.

    Parameters
    ----------
    top: Topology
        In .pdb or .mae format
    traj: Trajectory
        In .nc or .dcd format
    lipid_sele: str
        VMD-selection for lipids. If empty or None, attempt to determine lipid resname
    """
    if lipid_sele:
        evaltcl("atomselect macro lipid \" " + lipid_sele + " \"")

        # Determine resnames in selection
        molid = load_traj(top, traj, 0, 1, 1)
        evaltcl("set all_lipid [atomselect %s \" lipid \" frame 0]" % molid)
        lipid_resn = set(evaltcl("$all_lipid get resname").split())
        lipid_ids = get_atom_selection_indices("all_lipid")
        evaltcl("$all_lipid delete")
        molecule.delete(molid)
        if lipid_ids:
            print("Detected %d lipid atoms (resnames %s) matching --lipid '%s'" %
                  (len(lipid_ids), " ".join(set(lipid_resn)), lipid_sele))
        else:
            print("Detected no lipid atoms matching --lipid '%s'" % lipid_sele)
        return lipid_ids
    else:
        # If we need to expand default definition of lipids uncomment this line and add to this list
        # evaltcl("atomselect macro lipid \" (resname DLPE DMPC DPPC GPC LPPC PALM PC PGCL POPC POPE) \"")
        # return set("DLPE DMPC DPPC GPC LPPC PALM PC PGCL POPC POPE".split())

        lipid_resnames = set("DLPE DMPC DPPC GPC LPPC PALM PC PGCL POPC POPE".split())
        molid = load_traj(top, traj, 0, 1, 1)
        evaltcl("set all_atoms [atomselect %s \" all \" frame 0]" % molid)
        all_resn = set(evaltcl("$all_atoms get resname").split())
        lipid_resnames = lipid_resnames & all_resn
        evaltcl("$all_atoms delete")

        if not lipid_resnames:
            print("Detected no lipid atoms (specify manually using --lipid)")
            evaltcl("atomselect macro lipid \" none \"")
            molecule.delete(molid)
            return set([])
        else:
            lipid_resn = " ".join(lipid_resnames)
            evaltcl("atomselect macro lipid \" (resname " + lipid_resn + ") \"")
            lipid_ids = get_selection_indices(molid, 0, "lipid")
            print("Detected %d lipid atoms (resname %s)" % (len(lipid_ids), lipid_resn))
            molecule.delete(molid)
            return lipid_ids


def configure_ligand(top, traj, ligand_sele, sele1, sele2):
    """
    Detects the ligands if not specified in ligand_sele and creates a VMD selection macro called 'ligand'

    Parameters
    ----------
    top: Topology
        In .pdb or .mae format
    traj: Trajectory
        In .nc or .dcd format
    ligand_sele: str | None
        Ligand selection
    sele1: str
        Selection 1
    sele2: str
        Selection 2
    """
    molid = load_traj(top, traj, 0, 1, 1)
    if ligand_sele:
        try:
            evaltcl("atomselect macro ligand \" (" + ligand_sele + ") \"")
        except ValueError as err:
            print("Error:", err, "- not a valid VMD-selection.", file=sys.stderr)
            sys.exit(-1)

        evaltcl("set ligatoms [atomselect %s \" ligand \" frame 0]" % molid)
        ligatoms = get_selection_indices(molid, 0, "ligand")

        if ligatoms:
            print("Detected %d ligand atoms matching --ligand '%s'" % (len(ligatoms), ligand_sele))

            # Warn if there is no overlap between ligand and sele1 or sele2
            ligatoms = get_selection_indices(molid, 0, "ligand and (%s or %s)" % (sele1, sele2))
            if not ligatoms:
                print("Warning: None of the selected ligand atoms are in --sele or --sele2")
        else:
            print("Detected no ligand atoms matching --ligand '%s'" % ligand_sele)

    else:
        evaltcl("atomselect macro ligand \"not (lipid or solv or protein or nucleic)\"")
        evaltcl("set ligatoms [atomselect %s \" ligand \" frame 0]" % molid)
        ligatoms = get_atom_selection_indices("ligatoms")
        ligand_resnames = set(evaltcl("$ligatoms get resname").split())
        if ligatoms:
            ligand_resn = " ".join(ligand_resnames)
            print("Detected %d ligand atoms (resnames: %s)" % (len(ligatoms), ligand_resn))

            # Warn if there is no overlap between ligand and sele1 or sele2
            ligatoms = get_selection_indices(molid, 0, "ligand and (%s or %s)" % (sele1, sele2))
            if not ligatoms:
                print("Warning: No ligand atoms are in --sele or --sele2")
        else:
            print("Detected no ligand atoms (specify manually using e.g. --ligand 'resname GTP')")
            evaltcl("atomselect macro ligand \" none \"")

    evaltcl("$ligatoms delete")
    molecule.delete(molid)

def is_sp3(molid, index_to_atom, atom1, atom2, atom3):
    atom1 = index_to_atom[atom1].get_label()
    atom2 = index_to_atom[atom2].get_label()
    atom3 = index_to_atom[atom3].get_label()
    
    angle = compute_angle(molid, 0, atom1, atom2, atom3)
    return (109.5 - 5) < angle and angle < (109.5 + 5)

def is_sp2(molid, index_to_atom, atom1, atom2, atom3):
    atom1 = index_to_atom[atom1].get_label()
    atom2 = index_to_atom[atom2].get_label()
    atom3 = index_to_atom[atom3].get_label()

    angle = compute_angle(molid, 0, atom1, atom2, atom3)
    return (120. - 5) < angle and angle < (120. + 5)

def is_sp(molid, atom1, atom2, atom3):
    atom1 = index_to_atom[atom1].get_label()
    atom2 = index_to_atom[atom2].get_label()
    atom3 = index_to_atom[atom3].get_label()
    
    angle = compute_angle(molid, 0, atom1, atom2, atom3)
    return (180. - 5) < angle and angle < (180. + 5)

def extract_ligand_features(top, traj, index_to_atom):
    """
    Extracts lists of cationc and anionic atoms identified in the ligand

    Parameters
    ----------
    top: Topology
        In .pdb or .mae format
    traj: Trajectory
        In .nc or .dcd format
    index_to_atom: dict
        Maps VMD atom index to Atom

    Returns
    -------
    ligand_anions/ligand_cations: list
        VMD indices of anions/cations, respectively, detected among the ligand atoms
    """
    molid = load_traj(top, traj, 0, 1, 1)
    ligand_indices = get_selection_indices(molid, 0, "ligand")

    ligand_anions = []
    ligand_cations = []

    metal_cations = ["MG", "MN", "RH", "ZN", "FE", "BI", "AS", "AG"]

    ''' Find all neighbors '''
    index_to_neighbors = {}
    for atom_idx in ligand_indices:
        # print("Atom idx {}, element {}".format(atom_idx, index_to_atom[atom_idx].element))
        evaltcl("set neighbors [atomselect %s \"within 1.95 of (index %d)\" frame %s]" % (molid, atom_idx, 0))
        neighbor_indices = [idx for idx in get_atom_selection_indices("neighbors") if idx != atom_idx]
        evaltcl("$neighbors delete")
        index_to_neighbors[atom_idx] = neighbor_indices

    ''' Identify ligand cations/anions '''
    for atom_idx in ligand_indices:
        neighbors = index_to_neighbors[atom_idx] if atom_idx in index_to_neighbors else []
        ''' Check if the atom is a metal cation. I.E. one of the metal_cations names appears in its label. '''
        if any([cation in index_to_atom[atom_idx].get_label() for cation in metal_cations]):
            ligand_cations += [atom_idx]
            continue

        ''' Detect carboxylates in ligand. There are two restrictions:
        - sp2 carbon,
        - attached to 2 O's and 1 C
        '''
        if index_to_atom[atom_idx].element == 'C' and len(neighbors) == 3:
            # Check hybridization
            if not is_sp2(molid, index_to_atom, neighbors[0], atom_idx, neighbors[1]): continue
            # Check bonded atoms
            neighbor_elements = [index_to_atom[n_idx].element for n_idx in neighbors]
            from collections import Counter
            neighbor_elem_counts = Counter(neighbor_elements)
            if not (neighbor_elem_counts['C'] == 1 and neighbor_elem_counts['O'] == 2):
                continue
            # It's a carboxylate (probably)! Add both O's to ligand_anions
            ligand_anions += [n_idx for n_idx in neighbors if index_to_atom[n_idx].element == 'O']

    molecule.delete(molid)
    return ligand_anions, ligand_cations



# def compute_distance(molid, frame_idx, atom1_label, atom2_label):
#     """
#     Compute distance between two atoms in a specified frame of simulation

#     Parameters
#     ----------
#     molid: int
#         Denotes trajectory id in VMD
#     frame_idx: int
#         Frame of simulation in trajectory fragment
#     atom1_label: string
#         Atom label (ie "A:GLU:323:OE2:55124")
#     atom2_label: string
#         Atom label (ie "A:ARG:239:NH1:53746")

#     Returns
#     -------
#     distance: float
#     """
#     # atom_index1 = atom1_label.split(":")[-1]
#     # atom_index2 = atom2_label.split(":")[-1]
#     atom_index1 = atom1_label[atom1_label.rfind(":")+1:]
#     atom_index2 = atom2_label[atom2_label.rfind(":")+1:]
#     distance = float(evaltcl("measure bond {%s %s} molid %s frame %s" % (atom_index1, atom_index2, molid, frame_idx)))
#     return distance

def compute_distance(molid, frame_idx, atom1_index, atom2_index):
    """
    Compute distance between two atoms in a specified frame of simulation

    Parameters
    ----------
    molid: int
        Denotes trajectory id in VMD
    frame_idx: int
        Frame of simulation in trajectory fragment
    atom1_index: int
        Atom index 1
    atom2_index: int
        Atom index 2

    Returns
    -------
    distance: float
    """
    distance = float(evaltcl("measure bond {%s %s} molid %s frame %s" % (atom1_index, atom2_index, molid, frame_idx)))
    return distance


def compute_angle(molid, frame_idx, atom1, atom2, atom3):
    """
    Compute distance between two atoms in a specified frame of simulation

    Parameters
    ----------
    molid: int
        Denotes trajectory id in VMD
    frame_idx: int
        Frame of simulation in trajectory fragment
    atom1: string
        Atom label (ie "A:GLU:323:OE2:55124")
    atom2: string
        Atom label (ie "A:ARG:239:NH1:53746")
    atom3: string
        Atom label (ie "A:GLU:118:OE1:51792")

    Returns
    -------
    angle: float
        Expressed in degrees
    """
    atom_index1 = atom1.split(":")[-1]
    atom_index2 = atom2.split(":")[-1]
    atom_index3 = atom3.split(":")[-1]

    angle = float(evaltcl("measure angle {%s %s %s} molid %s frame %s" %
                          (atom_index1, atom_index2, atom_index3, molid, frame_idx)))
    return angle


def get_chain(traj_frag_molid, frame_idx, index):
    """
    Parse atom label and return element

    Parameters
    ----------
    traj_frag_molid: int
        Denotes trajectory id in VMD
    frame_idx: int
        Frame of simulation in trajectory fragment
    index: string
        VMD atom index

    Returns
    -------
    chain: string (ie "A", "B")
    """
    evaltcl("set sel [atomselect %s \" index %s \" frame %s]" % (traj_frag_molid, index, frame_idx))
    chain = evaltcl("$sel get chain")
    evaltcl("$sel delete")
    return chain


def get_resname(traj_frag_molid, frame_idx, index):
    """
    Parse atom label and return element

    Parameters
    ----------
    traj_frag_molid: int
        Denotes trajectory id in VMD
    frame_idx: int
        Frame of simulation in trajectory fragment
    index: string
        VMD atom index

    Returns
    -------
    resname: string (ie "ASP", "GLU")
    """
    evaltcl("set sel [atomselect %s \" index %s \" frame %s]" % (traj_frag_molid, index, frame_idx))
    resname = evaltcl("$sel get resname")
    evaltcl("$sel delete")
    return resname


def get_resid(traj_frag_molid, frame_idx, index):
    """
    Parse atom label and return element

    Parameters
    ----------
    traj_frag_molid: int
        Denotes trajectory id in VMD
    frame_idx: int
        Frame of simulation in trajectory fragment
    index: string
        VMD atom index

    Returns
    -------
    resid: string (ie "115", "117")
    """
    evaltcl("set sel [atomselect %s \" index %s \" frame %s]" % (traj_frag_molid, index, frame_idx))
    resid = evaltcl("$sel get resid")
    evaltcl("$sel delete")
    return resid


def get_name(traj_frag_molid, frame_idx, index):
    """
    Parse atom label and return element

    Parameters
    ----------
    traj_frag_molid: int
        Denotes trajectory id in VMD
    frame_idx: int
        Frame of simulation in trajectory fragment
    index: string
        VMD atom index

    Returns
    -------
    name: string (ie "CA", "NZ", )
    """
    evaltcl("set sel [atomselect %s \" index %s \" frame %s]" % (traj_frag_molid, index, frame_idx))
    name = evaltcl("$sel get name")
    evaltcl("$sel delete")
    return name


def get_element(traj_frag_molid, frame_idx, index):
    """
    Parse atom label and return element

    Parameters
    ----------
    traj_frag_molid: int
        Denotes trajectory id in VMD
    frame_idx: int
        Frame of simulation in trajectory fragment
    index: string
        VMD atom index

    Returns
    -------
    element: string (ie "C", "H", "O", "N", "S")
    """
    evaltcl("set sel [atomselect %s \" index %s \" frame %s]" % (traj_frag_molid, index, frame_idx))
    element = evaltcl("$sel get element")
    evaltcl("$sel delete")
    return element


def get_atom_label(traj_frag_molid, frame_idx, index):
    chain = get_chain(traj_frag_molid, frame_idx, index)
    resname = get_resname(traj_frag_molid, frame_idx, index)
    resid = get_resid(traj_frag_molid, frame_idx, index)
    name = get_name(traj_frag_molid, frame_idx, index)

    atom_label = "%s:%s:%s:%s:%s" % (chain, resname, resid, name, index)
    return atom_label


# def parse_contacts(contact_string):
#     """
#     Parameters
#     ----------
#     contact_string: string
#         Output from measure contacts function {indices} {indices}
#
#     Returns
#     -------
#     contact_index_pairs: list of tuples
#         List of index int pairs
#     """
#     contact_index_pairs = []
#
#     # Handle case where only one pair of atoms form contacts
#     if "} {" not in contact_string:
#         # print("CASE1: "+contact_string)
#         atom1_index, atom2_index = map(int, contact_string.split(" "))
#         contact_index_pairs.append((atom1_index, atom2_index))
#     else:
#         # print("CASE2: "+contact_string)
#         contacts_list = contact_string.split("} {")
#         atom1_list_str = contacts_list[0].strip("{}")
#         atom2_list_str = contacts_list[1].strip("{}")
#         if atom1_list_str == "" or atom2_list_str == "":
#             return []
#
#         atom1_list = list(map(int, atom1_list_str.split(" ")))
#         atom2_list = list(map(int, atom2_list_str.split(" ")))
#
#         for idx in range(len(atom1_list)):
#             atom1_index = atom1_list[idx]
#             atom2_index = atom2_list[idx]
#             contact_index_pairs.append((atom1_index, atom2_index))
#
#     return contact_index_pairs


def parse_contacts(contact_string):
    """
    Parameters
    ----------
    contact_string: string
        Output from measure contacts function {indices} {indices}

    Returns
    -------
    contact_index_pairs: list of tuples
        List of index int pairs
    """
    atom_indices = [int(index) for index in int_pattern.findall(contact_string)]
    half = len(atom_indices) // 2
    return zip(atom_indices[0:half], atom_indices[half:])


# Geometry Tools
def get_coord(traj_frag_molid, frame_idx, atom_label):
    """
    Get x, y, z coordinate of an atom specified by its label

    Parameters
    ----------
    traj_frag_molid: int
        Denotes trajectory id in VMD
    frame_idx: int
        Frame of simulation in trajectory fragment
    atom_label: string
        Atom label (ie "A:GLU:323:OE2:55124")

    Returns
    -------
    coord: np.array[x, y, z]
    """
    index = atom_label.split(":")[-1]
    evaltcl("set sel [atomselect %s \" index %s \" frame %s]" % (traj_frag_molid, index, frame_idx))
    # coord = np.array(list(map(float, float_pattern.findall(evaltcl("$sel get {x y z}")))))
    x = float(evaltcl("$sel get x"))
    y = float(evaltcl("$sel get y"))
    z = float(evaltcl("$sel get z"))
    evaltcl("$sel delete")
    coord = np.array([x, y, z])
    return coord


def points_to_vector(point1, point2):
    """
    Return vector from point1 to point2
    """
    return point2 - point1


def calc_vector_length(vector):
    """
    Compute length of vector
    """
    return norm(vector)
    # vector_length = math.sqrt(np.dot(vector, vector))
    # return vector_length


def calc_angle_between_vectors(vector1, vector2):
    """
    Returns
    -------
    angle_between_vectors: float
        Degrees between two vectors
    """
    radians_between_vectors = math.acos(np.dot(vector1, vector2) /
                                        (calc_vector_length(vector1) * calc_vector_length(vector2)))
    angle_between_vectors = math.degrees(radians_between_vectors)
    return angle_between_vectors


def calc_geom_distance(point1, point2):
    """
    Compute distance between two points

    Parameters
    ----------
    point1: array_like of number
        First point
    point2: array_like of number
        Second point

    Returns
    -------
    distance: float
    """
    distance = np.linalg.norm(point1 - point2)
    return distance


def calc_geom_centroid(point1, point2, point3):
    """
    Compute centroid between three points

    Parameters
    ----------
    point1: array_like of number
        First point
    point2: array_like of number
        Second point
    point3: array_like of number
        Second point

    Returns
    -------
    centroid: np.array[x, y, z]
    """
    centroid = (point1 + point2 + point3)/3
    return centroid


def calc_geom_normal_vector(point1, point2, point3):
    """
    Compute normal vector to the plane constructed by three points

    Returns
    -------
    normal_vector: np.array[x, y, z]

    """
    v1 = point3 - point1
    v2 = point2 - point1
    normal_vector = np.cross(v1, v2)
    return normal_vector


def calc_geom_psi_angle(center1, center2, normal_vector):
    """
    Parameters
    ----------
    center1: np.array of floats
        Coordinates of the center of aromatic plane 1
    center2: np.array of floats
        Coordinates of the center of aromatic plane 2
    normal_vector: np.array of floats
        Vector coordinate pointing from aromatic plane 1

    Returns
    -------
    psi_angle: float in degrees
        Angle between normal vector and vector(center2 - center1)
    """
    center_to_center_vector = points_to_vector(center2, center1)
    psi_angle = calc_angle_between_vectors(normal_vector, center_to_center_vector)
    psi_angle = min(math.fabs(psi_angle - 0), math.fabs(psi_angle - 180))
    return psi_angle
