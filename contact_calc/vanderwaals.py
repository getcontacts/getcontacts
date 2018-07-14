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

__all__ = ["compute_vanderwaals"]

##############################################################################
# Globals
##############################################################################

ALPHA_CARBON_DIST_CUTOFF = 10.0  # Angstroms
SOFT_VDW_CUTOFF = 5.0  # Angstroms
ATOM_RADIUS = {'H': 1.20,
               'C': 1.70,
               'N': 1.55,
               'O': 1.52,
               'S': 1.80,
               'P': 1.80}

##############################################################################
# Functions
##############################################################################


def filter_dual_selection_vdw(sele1_atoms, sele2_atoms, atom1_index, atom2_index):
    """
    Filter out van der waals interactions that are not between selection 1 and selection 2

    Parameters
    ----------
    sele1_atoms: list 
        List of atom label strings for all atoms in selection 1
    sele2_atoms: list 
        List of atom label strings for all atoms in selection 2
    atom1_index: int
        Index for atom1 participating in vdw interaction
    atom2_index: int
        Index for atom2 participating in vdw interaction

    Returns
    -------
    filter_bool: bool 
        True to filter out interaction
    """

    dual_sel1 = (atom1_index in sele1_atoms) and (atom2_index in sele2_atoms)
    if(dual_sel1):
        return False

    dual_sel2 = (atom1_index in sele2_atoms) and (atom2_index in sele1_atoms)
    if(dual_sel2):
        return False

    return True


def compute_vanderwaals(traj_frag_molid, frame_idx, index_to_label, sele_id, sele_id2,
                        sele1_atoms, sele2_atoms, ligands, VDW_EPSILON, VDW_RES_DIFF):
    """
    Compute all vanderwaals interactions in a frame of simulation

    Parameters
    ----------
    traj_frag_molid: int
        Identifier to simulation fragment in VMD
    frame_idx: int
        Frame number to query
    index_to_label: dict
        Maps VMD atom index to label "chain:resname:resid:name:index"
        {11205: "A:ASP:114:CA:11205, ...}
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
    if sele_id is None and sele_id2 is None:
        custom_sele = ""
    elif sele_id is not None and sele_id2 is None:
        custom_sele = "and (%s) " % sele_id
    elif sele_id is not None and sele_id2 is not None:
        sele_union = "(%s) or (%s)" % (sele_id, sele_id2)
        custom_sele = "and (%s) " % sele_union
    # custom_sele = "" if sele_id is None else "and (" + sele_id + ") "
    custom_lig = "" if not ligands else "or (resname " + (" ".join(ligands)) + ") "

    evaltcl("set full_protein [atomselect %s \" noh and ( protein %s) %s\" frame %s]" %
            (traj_frag_molid, custom_lig, custom_sele, frame_idx))

    # if sele_id is None:
    #     evaltcl("set full_protein [atomselect %s \" noh and protein \" frame %s]" %
    #             (traj_frag_molid, frame_idx))
    # else:
    #     evaltcl("set full_protein [atomselect %s \" noh and protein and (%s) \" frame %s]" %
    #             (traj_frag_molid, sele_id, frame_idx))

    contacts = evaltcl("measure contacts %s $full_protein" % SOFT_VDW_CUTOFF)
    contact_index_pairs = parse_contacts(contacts)
    evaltcl('$full_protein delete')

    vanderwaals = []
    for atom1_index, atom2_index in contact_index_pairs:
        # Perform dual selection with atom indices
        if sele1_atoms is not None and sele2_atoms is not None:
            if filter_dual_selection_vdw(sele1_atoms, sele2_atoms, atom1_index, atom2_index):
                continue

        # Convert to atom label
        atom1_label, atom2_label = index_to_label[atom1_index], index_to_label[atom2_index]

        # Another optimization is only convert to labels after passing initial checks. 
        # Use numeric comparison for sele1_atoms, sele2_atoms
        atom1_label_split = atom1_label.split(":")
        atom2_label_split = atom2_label.split(":")

        chain1 = atom1_label_split[0]
        chain2 = atom2_label_split[0]
        resi1 = int(atom1_label_split[2])
        resi2 = int(atom2_label_split[2])
        if chain1 == chain2 and abs(resi1-resi2) < VDW_RES_DIFF:
            continue

        element1 = atom1_label_split[3][0]
        element2 = atom2_label_split[3][0]

        if element1 not in ATOM_RADIUS:
            element1 = 'C'
        if element2 not in ATOM_RADIUS:
            element2 = 'C'

        # Perform distance cutoff with atom indices
        distance = compute_distance(traj_frag_molid, frame_idx, atom1_index, atom2_index)
        vanderwaal_cutoff = ATOM_RADIUS[element1] + ATOM_RADIUS[element2] + VDW_EPSILON
        if distance < vanderwaal_cutoff:
            vanderwaals.append([frame_idx, "vdw", atom1_label, atom2_label])

    return vanderwaals
