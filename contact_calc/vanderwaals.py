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

# http://periodictable.com/Properties/A/VanDerWaalsRadius.v.html
# https://www.webelements.com/iron/atom_sizes.html (for metals)
ATOM_RADIUS = {'H': 1.20,
               'C': 1.70,
               'N': 1.55,
               'O': 1.52,
               'F': 1.47,
               'LI': 1.82,
               'NA': 2.27,
               'MG': 1.73,
               'AL': 1.84,
               'P': 1.80,
               'S': 1.80,
               'CL': 1.75,
               'K': 2.75,
               'CA': 2.31,
               'CR': 2.45,
               'MN': 2.05,
               'FE': 2.44,
               'CO': 2.40,
               'NI': 2.40,
               'CU': 2.38,
               'ZN': 2.39,
               'BR': 1.85, 
               'I': 1.98
               }

##############################################################################
# Functions
##############################################################################

def infer_element(resname, atom_name):
    """
    Infer the element of an atom based on atom name 

    Parameters
    ----------
    resname: string
        Residue name as specified in topology
    atom_name: string 
        Atom name as specified in topology

    Returns
    -------
    element: string
        Inferred element identity from the atom name
    """

    # Single letter atom_name identifies element
    if(len(atom_name) == 1):
        return atom_name

    # Consider both single and double letter element names
    cand_elem1 = atom_name[0:1]
    cand_elem2 = atom_name[0:2]

    # If single letter matches element dictionary
    if(cand_elem1 in ATOM_RADIUS and cand_elem2 not in ATOM_RADIUS):
        return cand_elem1
    if(cand_elem2 in ATOM_RADIUS):
        # Special ambiguous case where CA can be either calcium or alpha carbon
        if(cand_elem2 == "CA"):
            if(resname[0:2] != "CA"):
                return "C"
        return cand_elem2

    # Otherwise default element is carbon
    return "C"


def compute_vanderwaals(traj_frag_molid, frame_idx, index_to_label, sele_id, sele_id2,
                        ligands, VDW_EPSILON, VDW_RES_DIFF):
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
    
    sel1 = "" if sele_id is None else "and (%s) " % (sele_id)
    sel2 = "" if sele_id2 is None else "and (%s) " % (sele_id2)
    custom_lig = "" if not ligands else "or (resname " + (" ".join(ligands)) + ") "
    evaltcl("set vdw_atoms1 [atomselect %s \" noh and ( protein or (hetero and not water and not lipid) %s) %s\" frame %s]" %
            (traj_frag_molid, custom_lig, sel1, frame_idx))
    evaltcl("set vdw_atoms2 [atomselect %s \" noh and ( protein or (hetero and not water and not lipid) %s) %s\" frame %s]" %
            (traj_frag_molid, custom_lig, sel2, frame_idx))

    if(sel2 == ""):
        contacts = evaltcl("measure contacts %s $vdw_atoms1" % (SOFT_VDW_CUTOFF))
    else:  
        contacts = evaltcl("measure contacts %s $vdw_atoms1 $vdw_atoms2" % (SOFT_VDW_CUTOFF))
    contact_index_pairs = parse_contacts(contacts)
    evaltcl("$vdw_atoms1 delete")
    evaltcl("$vdw_atoms2 delete")
    
    vanderwaals = []
    for atom1_index, atom2_index in contact_index_pairs:

        # Convert to atom label
        atom1_label, atom2_label = index_to_label[atom1_index], index_to_label[atom2_index]

        # Another optimization is only convert to labels after passing initial checks. 
        # Use numeric comparison for sele1_atoms, sele2_atoms
        atom1_label_split = atom1_label.split(":")
        atom2_label_split = atom2_label.split(":")

        chain1 = atom1_label_split[0]
        chain2 = atom2_label_split[0]
        resname1 = atom1_label_split[1]
        resname2 = atom2_label_split[1]
        resi1 = int(atom1_label_split[2])
        resi2 = int(atom2_label_split[2])
        if chain1 == chain2 and abs(resi1-resi2) < VDW_RES_DIFF:
            continue

        # Perform distance cutoff with atom indices
        distance = compute_distance(traj_frag_molid, frame_idx, atom1_index, atom2_index)

        # Infer element
        atom_name1 = atom1_label_split[3]
        atom_name2 = atom2_label_split[3]
        element1 = infer_element(resname1, atom_name1)
        element2 = infer_element(resname2, atom_name2)

        vanderwaal_cutoff = ATOM_RADIUS[element1] + ATOM_RADIUS[element2] + VDW_EPSILON
        if distance < vanderwaal_cutoff:
            vanderwaals.append([frame_idx, "vdw", atom1_label, atom2_label])

    return vanderwaals
