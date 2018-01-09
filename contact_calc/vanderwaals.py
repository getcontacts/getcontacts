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


def compute_vanderwaals(traj_frag_molid, frame_idx, index_to_label, sele_id, VDW_EPSILON):
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
    VDW_EPSILON: float, default = 0.5 angstroms
        amount of padding for calculating vanderwaals contacts

    Returns
    -------
    vanderwaals: list of tuples, [(frame_idx, itype, atom1_label, atom2_label), ...]
        itype = "vdw"
    """
    vanderwaals = []
    if sele_id is None:
        evaltcl("set full_protein [atomselect %s \" noh and protein \" frame %s]" % (traj_frag_molid, frame_idx))
    else:
        evaltcl("set full_protein [atomselect %s \" noh and protein and (%s) \" frame %s]" %
                (traj_frag_molid, sele_id, frame_idx))

    contacts = evaltcl("measure contacts %s $full_protein" % (SOFT_VDW_CUTOFF))
    contact_index_pairs = parse_contacts(contacts)
    evaltcl('$full_protein delete')
    for atom1_index, atom2_index in contact_index_pairs:
        atom1_label, atom2_label = index_to_label[atom1_index], index_to_label[atom2_index]
        element1 = atom1_label.split(":")[3][0]
        element2 = atom2_label.split(":")[3][0]

        distance = compute_distance(traj_frag_molid, frame_idx, atom1_label, atom2_label)
        vanderwaal_cutoff = ATOM_RADIUS[element1] + ATOM_RADIUS[element2] + VDW_EPSILON
        if distance < vanderwaal_cutoff:
            vanderwaals.append([frame_idx, "vdw", atom1_label, atom2_label])

    return vanderwaals
