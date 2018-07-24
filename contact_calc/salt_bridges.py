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


def compute_salt_bridges(traj_frag_molid, frame_idx, index_to_atom, sele_id, sele_id2, sele1_atoms, sele2_atoms,
                         SALT_BRIDGE_CUTOFF_DISTANCE=4.0):
    """
    Compute salt bridges in a frame of simulation

    Parameters
    ----------
    traj_frag_molid: int
        Identifier to simulation fragment in VMD
    frame_idx: int
        Frame number to query
    index_to_atom: dict
        Maps VMD atom index to label "chain:resname:resid:name:index"
        {11205: "A:ASP:114:CA:11205, ...}
    sele_id: string, default = None
        Compute contacts on subset of atom selection based on VMD query
    sele_id2: string, default = None
        If second VMD query is specified, then compute contacts between atom selection 1 and 2
    sele1_atoms: list 
        List of atom label indices for all atoms in selection 1
    sele2_atoms: list 
        List of atom label indices for all atoms in selection 2
    SALT_BRIDGE_CUTOFF_DISTANCE: float, default = 4.0 angstroms
        cutoff for distance between anion and cation atoms

    Returns
    -------
    salt_bridges: list of tuples, [(frame_index, itype, atom1_label, atom2_label), ...]
        itype = "sb"
    """
    anion_set, cation_set = prep_salt_bridge_computation(traj_frag_molid, frame_idx, sele_id, sele_id2)

    salt_bridges = []
    for anion_atom in anion_set:
        for cation_atom in cation_set:
            # Process dual selection output if user provides two selection queries
            if sele1_atoms is not None and sele2_atoms is not None:
                if filter_dual_selection_salt_bridges(sele1_atoms, sele2_atoms, anion_atom, cation_atom):
                    continue
                    
            dist = compute_distance(traj_frag_molid, frame_idx, anion_atom, cation_atom)
            if dist < SALT_BRIDGE_CUTOFF_DISTANCE:
                anion_label = index_to_atom[anion_atom].get_label()
                cation_label = index_to_atom[cation_atom].get_label()
                salt_bridges.append([frame_idx, "sb", anion_label, cation_label])

    return salt_bridges
