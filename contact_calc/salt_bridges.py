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

__all__ = ['compute_salt_bridges']

acidic_asp = "((resname ASP) and (name OD1 OD2))"
acidic_glu = "((resname GLU) and (name OE1 OE2))"
acidic_nucl = "((resname C U G A DC DT DG DA) and (name OP1 OP2))"
basic_his = "((resname HIS HSD HSE HSP HIE HIP HID) and (name ND1 NE2))"
basic_lys = "((resname LYS) and (name NZ))"
basic_arg = "((resname ARG) and (name NH1 NH2))"


def compute_salt_bridges(traj_frag_molid, frame, index_to_atom, sele1, sele2, geom_criteria, ligand_anions, ligand_cations):
    """
    Compute salt bridges in a frame of simulation

    Parameters
    ----------
    traj_frag_molid: int
        Identifier to simulation fragment in VMD
    frame: int
        Frame number to query
    index_to_atom: dict
        Maps VMD atom index to label "chain:resname:resid:name:index"
        {11205: "A:ASP:114:CA:11205, ...}
    sele1: string, default = None
        Compute contacts on subset of atom selection based on VMD query
    sele2: string, default = None
        If second VMD query is specified, then compute contacts between atom selection 1 and 2
    geom_criteria: dict
        Container for geometric criteria

    Returns
    -------
    salt_bridges: list of tuples, [(frame_index, itype, atom1_label, atom2_label), ...]
        itype = "sb"
    """
    cutoff_dist = geom_criteria['SALT_BRIDGE_CUTOFF_DISTANCE']

    ligand_indices = get_selection_indices(traj_frag_molid, frame, "ligand")
    ligand_anions_string, ligand_cations_string = "", ""
    if len(ligand_anions) > 0:
        ligand_anions_string = "or (index {})".format(' '.join([str(idx) for idx in ligand_anions]))
    if len(ligand_cations) > 0:
        ligand_cations_string = "or (index {})".format(' '.join([str(idx) for idx in ligand_cations]))

    s1_anions = "(%s or %s or %s %s) and (%s)" % (acidic_asp, acidic_glu, acidic_nucl, ligand_anions_string, sele1)
    s2_anions = "(%s or %s or %s %s) and (%s)" % (acidic_asp, acidic_glu, acidic_nucl, ligand_anions_string, sele2)
    s1_cations = "(%s or %s or %s %s) and (%s)" % (basic_his, basic_lys, basic_arg, ligand_cations_string, sele1)
    s2_cations = "(%s or %s or %s %s) and (%s)" % (basic_his, basic_lys, basic_arg, ligand_cations_string, sele2)

    evaltcl("set s1anions [atomselect %s \" %s \" frame %s]" % (traj_frag_molid, s1_anions, frame))
    evaltcl("set s2anions [atomselect %s \" %s \" frame %s]" % (traj_frag_molid, s2_anions, frame))
    evaltcl("set s1cations [atomselect %s \" %s \" frame %s]" % (traj_frag_molid, s1_cations, frame))
    evaltcl("set s2cations [atomselect %s \" %s \" frame %s]" % (traj_frag_molid, s2_cations, frame))
    contacts_12 = set(parse_contacts(evaltcl("measure contacts %f $s1anions $s2cations" % cutoff_dist)))
    if sele1 == sele2:
        contacts_21 = set([])
    else:
        contacts_21 = set(parse_contacts(evaltcl("measure contacts %f $s2cations $s1anions" % cutoff_dist)))
    evaltcl("$s1anions delete")
    evaltcl("$s2anions delete")
    evaltcl("$s1cations delete")
    evaltcl("$s2cations delete")

    acpairs = contacts_12 | contacts_21

    salt_bridges = []
    for (a, c) in acpairs:
        if a in ligand_indices and c not in ligand_indices:
            salt_bridges.append([frame, "sblp", index_to_atom[a].get_label(), index_to_atom[c].get_label()])
        if a not in ligand_indices and c in ligand_indices:
            salt_bridges.append([frame, "sbpl", index_to_atom[a].get_label(), index_to_atom[c].get_label()])
        if a in ligand_indices and c in ligand_indices:
            salt_bridges.append([frame, "sbll", index_to_atom[a].get_label(), index_to_atom[c].get_label()])
        else:
            salt_bridges.append([frame, "sb", index_to_atom[a].get_label(), index_to_atom[c].get_label()])
    # return [[frame, "sb", index_to_atom[a].get_label(), index_to_atom[c].get_label()] for (a, c) in acpairs]
    return salt_bridges

