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

from .contact_utils import *

__all__ = ["compute_vanderwaals"]

SOFT_VDW_CUTOFF = 0


def update_soft_cutoff(traj_frag_molid, index_to_atom, sele1, sele2, epsilon):
    global SOFT_VDW_CUTOFF
    all_ids = get_selection_indices(traj_frag_molid, 0, "(%s) or (%s)" % (sele1, sele2))
    max_vdw = max([index_to_atom[i].vdwradius for i in all_ids])
    SOFT_VDW_CUTOFF = max_vdw * 2 + epsilon


def compute_vanderwaals(traj_frag_molid, frame, index_to_atom, sele1, sele2, geom_criteria):
    """
    Compute all vanderwaals interactions in a frame of simulation

    Parameters
    ----------
    traj_frag_molid: int
        Identifier to simulation fragment in VMD
    frame: int
        Frame number to query
    index_to_atom: dict
        Maps VMD atom index to Atom
    sele1: string, default = None
        Compute contacts on subset of atom selection based on VMD query
    sele2: string, default = None
        If second VMD query is specified, then compute contacts between atom selection 1 and 2
    geom_criteria: dict
        Container for geometric criteria

    Returns
    -------
    vanderwaals: list of tuples, [(frame_idx, itype, atom1_label, atom2_label), ...]
        itype = "vdw"
    """
    epsilon = geom_criteria['VDW_EPSILON']
    res_diff = geom_criteria['VDW_RES_DIFF']
    if SOFT_VDW_CUTOFF == 0:
        update_soft_cutoff(traj_frag_molid, index_to_atom, sele1, sele2, epsilon)

    evaltcl("set vdw_atoms1 [atomselect %s \" noh and (%s)\" frame %s]" % (traj_frag_molid, sele1, frame))
    evaltcl("set vdw_atoms2 [atomselect %s \" noh and (%s)\" frame %s]" % (traj_frag_molid, sele2, frame))

    contacts = evaltcl("measure contacts %s $vdw_atoms1 $vdw_atoms2" % SOFT_VDW_CUTOFF)

    contact_index_pairs = parse_contacts(contacts)
    evaltcl("$vdw_atoms1 delete")
    evaltcl("$vdw_atoms2 delete")
    
    vanderwaals = []
    for atom1_index, atom2_index in contact_index_pairs:

        # Convert to atom label
        atom1, atom2 = index_to_atom[atom1_index], index_to_atom[atom2_index]

        if atom1.chain == atom2.chain and abs(atom1.resid - atom2.resid) < res_diff:
            continue

        # Perform distance cutoff with atom indices
        distance = compute_distance(traj_frag_molid, frame, atom1_index, atom2_index)

        # vanderwaal_cutoff = ATOM_RADIUS[element1] + ATOM_RADIUS[element2] + VDW_EPSILON
        vanderwaal_cutoff = atom1.vdwradius + atom2.vdwradius + epsilon
        if distance < vanderwaal_cutoff:
            vanderwaals.append([frame, "vdw", atom1.get_label(), atom2.get_label()])

    return vanderwaals
