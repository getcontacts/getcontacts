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

__all__ = ['prep_salt_bridge_computation', 'compute_salt_bridges']

##############################################################################
# Functions
##############################################################################


def filter_dual_selection_salt_bridges(sele1_atoms, sele2_atoms, anion_atom, cation_atom):
    """
    Filter out salt bridge interaction that is not between selection 1 and selection 2

    Parameters
    ----------
    sele1_atoms: list 
        List of atom label strings for all atoms in selection 1
    sele2_atoms: list 
        List of atom label strings for all atoms in selection 2
    anion_atom: string 
        Atom label for anion
    cation_atom: string 
        Atom label for cation

    """
    dual_sel1 = (anion_atom in sele1_atoms and cation_atom in sele2_atoms)
    if(dual_sel1):
        return False

    dual_sel2 = (anion_atom in sele2_atoms and cation_atom in sele1_atoms)
    if(dual_sel2):
        return False
    return True


def prep_salt_bridge_computation(traj_frag_molid, frame_idx, sele_id, sele_id2):
    """
    Compute all possible anion and cation atoms from first frame of simulation

    Returns
    -------
    anion_list: list of strings
        List of atom labels for atoms in ASP or GLU that
        can form salt bridges
    cation_list: list of strings
        List of atom labels for atoms in LYS, ARG, HIS that
        can form salt bridges
    """
    anion_set = get_anion_atoms(traj_frag_molid, frame_idx, sele_id, sele_id2)
    cation_set = get_cation_atoms(traj_frag_molid, frame_idx, sele_id, sele_id2)
    return anion_set, cation_set


def compute_salt_bridges(traj_frag_molid, frame, index_to_atom, sele1, sele2, geom_criteria):
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
    sele1_atoms: list 
        List of atom label indices for all atoms in selection 1
    sele2_atoms: list 
        List of atom label indices for all atoms in selection 2
    geom_criteria: dict
        Container for geometric criteria

    Returns
    -------
    salt_bridges: list of tuples, [(frame_index, itype, atom1_label, atom2_label), ...]
        itype = "sb"
    """
    cutoff_dist = geom_criteria['SALT_BRIDGE_CUTOFF_DISTANCE']

    acidic_asp = "((resname ASP) and (name OD1 OD2))"
    acidic_glu = "((resname GLU) and (name OE1 OE2))"
    s1_anions = "(%s or %s) and %s" % (acidic_asp, acidic_glu, sele1)
    s2_anions = "(%s or %s) and %s" % (acidic_asp, acidic_glu, sele2)

    basic_his = "((resname HIS HSD HSE HSP HIE HIP HID) and (name ND1 NE2))"
    basic_lys = "((resname LYS) and (name NZ))"
    basic_arg = "((resname ARG) and (name NH1 NH2))"
    s1_cations = "( %s or %s or %s) and %s" % (basic_his, basic_lys, basic_arg, sele1)
    s2_cations = "( %s or %s or %s) and %s" % (basic_his, basic_lys, basic_arg, sele2)

    evaltcl("set s1anions [atomselect %s \" %s \" frame %s]" % (traj_frag_molid, s1_anions, frame))
    evaltcl("set s2anions [atomselect %s \" %s \" frame %s]" % (traj_frag_molid, s2_anions, frame))
    evaltcl("set s1cations [atomselect %s \" %s \" frame %s]" % (traj_frag_molid, s1_cations, frame))
    evaltcl("set s2cations [atomselect %s \" %s \" frame %s]" % (traj_frag_molid, s2_cations, frame))
    contacts_12 = set(parse_contacts(evaltcl("measure contacts %f $s1anions $s2cations" % cutoff_dist)))
    contacts_21 = set(parse_contacts(evaltcl("measure contacts %f $s2anions $s1cations" % cutoff_dist)))
    evaltcl("$s1anions delete")
    evaltcl("$s2anions delete")
    evaltcl("$s1cations delete")
    evaltcl("$s2cations delete")

    acpairs = contacts_12 | contacts_21
    return [[frame, "sb", index_to_atom[a].get_label(), index_to_atom[c].get_label()] for (a, c) in acpairs]

