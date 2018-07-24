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

##############################################################################
# Imports
##############################################################################

from .contact_utils import *
from .atom import *

__all__ = ["compute_vanderwaals"]

##############################################################################
# Globals
##############################################################################

ALPHA_CARBON_DIST_CUTOFF = 10.0  # Angstroms
SOFT_VDW_CUTOFF = 5.0  # Angstroms

# # http://periodictable.com/Properties/A/VanDerWaalsRadius.v.html
# ATOM_RADIUS = {'H': 1.20,
#                'C': 1.70,
#                'N': 1.55,
#                'O': 1.52,
#                'S': 1.80,
#                'P': 1.80,
#                'K': 2.75,
#                'NA': 2.27,
#                'MG': 1.73,
#                'LI': 1.82,
#                'F': 1.47,
#                'CL': 1.75}

##############################################################################
# Functions
##############################################################################

# def infer_element(atom_name):
#     """
#     Infer the element of an atom based on atom name
#
#     Parameters
#     ----------
#     atom_name: string
#         Atom name as specified in topology
#
#     Returns
#     -------
#     element: string
#         Inferred element identity from the atom name
#     """
#
#     # Single letter atom_name identifies element
#     if(len(atom_name) == 1):
#         return atom_name
#
#     # Otherwise default element is carbon
#     element = "C"
#
#     # Consider both single and double letter element names
#     cand_elem1 = atom_name[0:1]
#     cand_elem2 = atom_name[0:2]
#
#     # If single letter matches element dictionary
#     if(cand_elem1 in ATOM_RADIUS and cand_elem2 not in ATOM_RADIUS):
#         element = cand_elem1
#     if(cand_elem2 in ATOM_RADIUS):
#         element = cand_elem2
#
#     return element


def compute_vanderwaals(traj_frag_molid, frame_idx, index_to_atom, sele_id, sele_id2, ligands,
                        VDW_EPSILON, VDW_RES_DIFF):
    """
    Compute all vanderwaals interactions in a frame of simulation

    Parameters
    ----------
    traj_frag_molid: int
        Identifier to simulation fragment in VMD
    frame_idx: int
        Frame number to query
    index_to_atom: dict
        Maps VMD atom index to Atom
    sele_id: string, default = None
        Compute contacts on subset of atom selection based on VMD query
    sele_id2: string, default = None
        If second VMD query is specified, then compute contacts between atom selection 1 and 2
    sele1_atoms: list 
        List of atom label strings for all atoms in selection 1
    sele2_atoms: list 
        List of atom label strings for all atoms in selection 2
    ligands: list of string
        Residue names of ligands
    VDW_EPSILON: float, default = 0.5 angstroms
        amount of padding for calculating vanderwaals contacts
    VDW_RES_DIFF: int, default = 2
        minimum residue distance for which to consider computing vdw interactions

    Returns
    -------
    vanderwaals: list of tuples, [(frame_idx, itype, atom1_label, atom2_label), ...]
        itype = "vdw"
    """
    
    sel1 = "" if sele_id is None else "and (%s) " % (sele_id)
    sel2 = "" if sele_id2 is None else "and (%s) " % (sele_id2)
    custom_lig = "" if not ligands else "or (resname " + (" ".join(ligands)) + ") "
    evaltcl("set vdw_atoms1 [atomselect %s \" noh and ( protein %s) %s\" frame %s]" %
            (traj_frag_molid, custom_lig, sel1, frame_idx))
    evaltcl("set vdw_atoms2 [atomselect %s \" noh and ( protein %s) %s\" frame %s]" %
            (traj_frag_molid, custom_lig, sel2, frame_idx))

    if sel2 == "":
        contacts = evaltcl("measure contacts %s $vdw_atoms1" % (SOFT_VDW_CUTOFF))
    else:  
        contacts = evaltcl("measure contacts %s $vdw_atoms1 $vdw_atoms2" % (SOFT_VDW_CUTOFF))
    contact_index_pairs = parse_contacts(contacts)
    evaltcl("$vdw_atoms1 delete")
    evaltcl("$vdw_atoms2 delete")
    
    vanderwaals = []
    for atom1_index, atom2_index in contact_index_pairs:

        # Convert to atom label
        atom1, atom2 = index_to_atom[atom1_index], index_to_atom[atom2_index]

        if atom1.chain == atom2.chain and abs(atom1.resid - atom2.resid) < VDW_RES_DIFF:
            continue

        # Perform distance cutoff with atom indices
        distance = compute_distance(traj_frag_molid, frame_idx, atom1_index, atom2_index)

        # vanderwaal_cutoff = ATOM_RADIUS[element1] + ATOM_RADIUS[element2] + VDW_EPSILON
        vanderwaal_cutoff = atom1.vdwradius + atom2.vdwradius + VDW_EPSILON
        if distance < vanderwaal_cutoff:
            vanderwaals.append([frame_idx, "vdw", atom1.get_label(), atom2.get_label()])

    return vanderwaals
