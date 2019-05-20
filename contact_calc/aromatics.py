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

from vmd import *
from .contact_utils import *

__all__ = ['compute_pi_stacking', 'compute_t_stacking']

PI_STACK_SOFT_DISTANCE_CUTOFF = 10.0  # angstroms
T_STACK_SOFT_DISTANCE_CUTOFF = 6.0  # angstroms

aromatic_phe = "resname PHE and name CG CE1 CE2"
aromatic_trp = "resname TRP and name CD2 CZ2 CZ3"
aromatic_tyr = "resname TYR and name CG CE1 CE2"
aromatic_his = "(resname HIS HSD HSE HSP HIE HIP HID) and (name CG CE1 CD2)"
aromatic_nucl = "nucleic and name C4 C2 C6"

############################################################################
# Functions
############################################################################


def filter_dual_selection_aromatic(sele1_atoms, sele2_atoms, aromatic1_index, aromatic2_index):
    """
    Filter out aromatic interactions that are not between selection 1 and selection 2

    Parameters
    ----------
    sele1_atoms: list 
        List of atom label strings for all atoms in selection 1
    sele2_atoms: list 
        List of atom label strings for all atoms in selection 2
    aromatic1_index: int
        Aromatic atom index 1
    aromatic2_index: int
        Aromatic atom index 2

    Returns
    -------
    filter_bool: bool 
        True to filter out interaction

    """

    if (aromatic1_index in sele1_atoms) and (aromatic2_index in sele2_atoms):
        return False

    if (aromatic1_index in sele2_atoms) and (aromatic2_index in sele1_atoms):
        return False

    return True


def get_aromatic_triplet(molid, frame, aatom, index_to_atom):
    """
    Given an aromatic residue label return triplet of atoms on the ring

    Parameters
    ----------
    molid: int
        Identifier to simulation fragment in VMD
    frame: int
        Frame number to query
    aatom: Atom
        Aromatic atom
    index_to_atom: dict
        Maps VMD atom index to Atom

    Returns
    -------
    aromatic_atom_triplet_labels: list of strings
        ie ["A:PHE:329:CE1:55228", "A:PHE:329:CE2:55234", "A:PHE:329:CG:55225"]
    """
    residue_to_atom_names = {
        "PHE": "CG CE1 CE2",
        "TRP": "CD2 CZ2 CZ3",
        "TYR": "CG CE1 CE2",
        "HIS": "CG CE1 CD2 NE2 NE1",
        "HSD": "CG CE1 CD2 NE2 NE1",
        "HSE": "CG CE1 CD2 NE2 NE1",
        "HSP": "CG CE1 CD2 NE2 NE1",
        "HIE": "CG CE1 CD2 NE2 NE1",
        "HIP": "CG CE1 CD2 NE2 NE1", 
        "HID": "CG CE1 CD2 NE2 NE1",
    }
    chain, resname, resid = aatom.chain, aatom.resname, aatom.resid
    # Use the dict above and assume nucleic acid if residue name is not present
    arom_triple = residue_to_atom_names[resname] if resname in residue_to_atom_names else "C2 C4 C6"
    evaltcl("set aromatic_atoms [atomselect %s \"(chain '%s') and (resname '%s') and (resid '%s') and (name %s)\" frame %s]"
            % (molid, chain, resname, resid, arom_triple, frame))
    aromatic_atom_triplet_indices = get_atom_selection_indices("aromatic_atoms")
    evaltcl('$aromatic_atoms delete')

    # Need aromatic_atom_triplet as list 
    aromatic_atom_triplet_labels = [index_to_atom[atom_index].get_label() for atom_index in aromatic_atom_triplet_indices]
    return aromatic_atom_triplet_labels


def compute_aromatics(molid, frame, index_to_atom, sele1, sele2, sele1_atoms, sele2_atoms, itype,
                      soft_distance_cutoff, distance_cutoff, angle_cutoff, psi_angle_cutoff):
    """
    Compute aromatic interactions in a frame of simulation

    Parameters
    ----------
    molid: int
        Identifier to simulation fragment in VMD
    frame: int
        Frame number to query
    index_to_atom: dict
        Maps VMD atom index to Atom
    sele1: string, default = None
        Compute contacts on subset of atom selection based on VMD query
    sele2: string, default = None
        If second VMD query is specified, then compute contacts between atom selection 1 and 2 
    sele1_atoms: list 
        List of atom label strings for all atoms in selection 1
    sele2_atoms: list 
        List of atom label strings for all atoms in selection 2
    itype: string
        Specify which type of aromatics ("ps" or "ts") to compute
    soft_distance_cutoff: float
        Soft distance cutoff to find candidate aromatic pairs
    distance_cutoff: float
        Cutoff distance between aromatic centers
    angle_cutoff: float
        Cutoff angle for the angle between normal vectors of aromatic planes
    psi_angle_cutoff: float
        Cutoff angle for how aligned two aromatic planes are

    Returns
    -------
    aromatics = list of tuples, [(frame_index, atom1_label, atom2_label, itype), ...]
        itype = "ps" or "ts"
    """

    aromatics = []

    sele_union = "(%s) or (%s)" % (sele1, sele2)
    aromatic_atom_sel = "set aromatic_atoms [atomselect %s \" ((%s) or (%s) or (%s) or (%s) or (%s)) and (%s) \" frame %s]" % \
                        (molid, aromatic_phe, aromatic_trp, aromatic_tyr, aromatic_his, aromatic_nucl, sele_union, frame)
    # aromatic_atom_sel = "set aromatic_atoms [atomselect %s \"" \
    #                     "((resname PHE) and (name CG CE1 CE2) and (%s)) or " \
    #                     "((resname TRP) and (name CD2 CZ2 CZ3) and (%s)) or " \
    #                     "((resname TYR) and (name CG CE1 CE2) and (%s)) \" frame %s]" % \
    #                     (molid, sele_union, sele_union, sele_union, frame)

    evaltcl(aromatic_atom_sel)
    contacts = evaltcl("measure contacts %s $aromatic_atoms" % soft_distance_cutoff)
    evaltcl("$aromatic_atoms delete")

    # Calculate set of distinct aromatic candidate pairs that may have pi-stacking
    contact_index_pairs = parse_contacts(contacts)
    res_pairs = set()
    residue_to_atom_labels = {}  # map from aromatic residues that passed soft cutoff to their three atom labels
    for aromatic1_index, aromatic2_index in contact_index_pairs:

        # Perform dual selection with integer indices
        if sele1_atoms is not None and sele2_atoms is not None:
            if filter_dual_selection_aromatic(sele1_atoms, sele2_atoms, aromatic1_index, aromatic2_index):
                continue 

        aromatic1_label = index_to_atom[aromatic1_index].get_label()
        aromatic2_label = index_to_atom[aromatic2_index].get_label()

        # Check if the two atoms belong to same aromatic group
        aromatic1_res = ":".join(aromatic1_label.split(":")[0:3])
        aromatic2_res = ":".join(aromatic2_label.split(":")[0:3])
        if aromatic1_res == aromatic2_res:
            continue

        if aromatic1_res not in residue_to_atom_labels:
            residue_to_atom_labels[aromatic1_res] = get_aromatic_triplet(molid, frame, index_to_atom[aromatic1_index], index_to_atom)
        if aromatic2_res not in residue_to_atom_labels:
            residue_to_atom_labels[aromatic2_res] = get_aromatic_triplet(molid, frame, index_to_atom[aromatic2_index], index_to_atom)

        k1 = (aromatic1_res, aromatic2_res)
        k2 = (aromatic2_res, aromatic1_res)
        if k1 not in res_pairs and k2 not in res_pairs:
            res_pairs.add((aromatic1_res, aromatic2_res))

    # Perform strict geometric criterion on candidate aromatic pairs
    for aromatic1_res, aromatic2_res in res_pairs:
        aromatic1_atom_labels = residue_to_atom_labels[aromatic1_res]
        aromatic2_atom_labels = residue_to_atom_labels[aromatic2_res]

        # Ignore interactions between residues with missing atoms
        if len(aromatic1_atom_labels) != 3 or len(aromatic2_atom_labels) != 3:
            continue

        # Distance between two aromatic centers must be below DISTANCE_CUTOFF
        arom1_atom1_coord = get_coord(molid, frame, aromatic1_atom_labels[0])
        arom1_atom2_coord = get_coord(molid, frame, aromatic1_atom_labels[1])
        arom1_atom3_coord = get_coord(molid, frame, aromatic1_atom_labels[2])

        arom2_atom1_coord = get_coord(molid, frame, aromatic2_atom_labels[0])
        arom2_atom2_coord = get_coord(molid, frame, aromatic2_atom_labels[1])
        arom2_atom3_coord = get_coord(molid, frame, aromatic2_atom_labels[2])

        aromatic1_centroid = calc_geom_centroid(arom1_atom1_coord, arom1_atom2_coord, arom1_atom3_coord)
        aromatic2_centroid = calc_geom_centroid(arom2_atom1_coord, arom2_atom2_coord, arom2_atom3_coord)
        aromatic_centers_distance = calc_geom_distance(aromatic1_centroid, aromatic2_centroid)
        if aromatic_centers_distance > distance_cutoff:
            continue

        # print(aromatic1_atom_labels, aromatic2_atom_labels)
        # Angle between vectors normal to each aromatic plane must be below cutoff
        aromatic1_normal_vector = calc_geom_normal_vector(arom1_atom1_coord, arom1_atom2_coord, arom1_atom3_coord)
        aromatic2_normal_vector = calc_geom_normal_vector(arom2_atom1_coord, arom2_atom2_coord, arom2_atom3_coord)
        aromatic_plane_alignment_angle = calc_angle_between_vectors(aromatic1_normal_vector, aromatic2_normal_vector)
        if itype == "ps":
            aromatic_plane_alignment_angle = min(math.fabs(aromatic_plane_alignment_angle - 0),
                                                 math.fabs(aromatic_plane_alignment_angle - 180))
            if aromatic_plane_alignment_angle > angle_cutoff:
                continue
        elif itype == "ts":
            aromatic_plane_perpendicular_angle = \
                math.fabs(calc_angle_between_vectors(aromatic1_normal_vector, aromatic2_normal_vector) - 90)
            if aromatic_plane_perpendicular_angle > angle_cutoff:
                continue
        # print(aromatic1_atom_labels, aromatic2_atom_labels)
        # Psi Angle cutoff
        psi_angle1 = calc_geom_psi_angle(aromatic1_centroid, aromatic2_centroid, aromatic1_normal_vector)
        psi_angle2 = calc_geom_psi_angle(aromatic2_centroid, aromatic1_centroid, aromatic2_normal_vector)
        psi_angle = min(psi_angle1, psi_angle2)
        if psi_angle > psi_angle_cutoff:
            continue

        # print(aromatic1_atom_labels, aromatic2_atom_labels)
        # Returns a single interaction between the CG atom of each aromatic ring
        arom1_CG_label = convert_to_single_atom_aromatic_string(aromatic1_atom_labels)
        arom2_CG_label = convert_to_single_atom_aromatic_string(aromatic2_atom_labels)
        # print(arom1_CG_label, arom2_CG_label)
        aromatics.append([frame, itype, arom1_CG_label, arom2_CG_label])

        # Returns all 3x3 combinations of aromatic interaction pairs (DEPRECATED)
        # for arom1_atom_label in aromatic1_atom_labels:
        #     for arom2_atom_label in aromatic2_atom_labels:
        #         aromatics.append([frame_idx, itype, arom1_atom_label, arom2_atom_label])
    return aromatics


def compute_pi_stacking(traj_frag_molid, frame_idx, index_to_atom, sele1, sele2, sele1_atoms, sele2_atoms,
                        geom_criteria):
    """
    Compute pi-stacking interactions in a frame of simulation

    Parameters
    ----------
    traj_frag_molid: int
        Identifier to simulation fragment in VMD
    frame_idx: int
        Frame number to query
    index_to_atom: dict
        Maps VMD atom index to Atom
    sele1: string, default = None
        Compute contacts on subset of atom selection based on VMD query
    sele2: string, default = None
        If second VMD query is specified, then compute contacts between atom selection 1 and 2
    sele1_atoms: list 
        List of atom label strings for all atoms in selection 1
    sele2_atoms: list 
        List of atom label strings for all atoms in selection 2
    geom_criteria: dict
        Container for geometric criteria

    Returns
    -------
    pi_stacking = list of tuples, [(frame_index, itype, atom1_label, atom2_label), ...]
        itype = "ps"
    """
    cutoff_distance = geom_criteria['PI_STACK_CUTOFF_DISTANCE']
    cutoff_angle = geom_criteria['PI_STACK_CUTOFF_ANGLE']
    psi_angle = geom_criteria['PI_STACK_PSI_ANGLE']
    pi_stacking = compute_aromatics(traj_frag_molid, frame_idx, index_to_atom,
                                    sele1, sele2,
                                    sele1_atoms, sele2_atoms, "ps",
                                    PI_STACK_SOFT_DISTANCE_CUTOFF, cutoff_distance,
                                    cutoff_angle, psi_angle)
    return pi_stacking


def compute_t_stacking(traj_frag_molid, frame_idx, index_to_atom, sele1, sele2, sele1_atoms, sele2_atoms,
                       geom_criteria):
    """
    Compute t-stacking interactions in a frame of simulation

    Parameters
    ----------
    traj_frag_molid: int
        Identifier to simulation fragment in VMD
    frame_idx: int
        Frame number to query
    index_to_atom: dict
        Maps VMD atom index to Atom
    sele1: string, default = None
        Compute contacts on subset of atom selection based on VMD query
    sele2: string, default = None
        If second VMD query is specified, then compute contacts between atom selection 1 and 2
    sele1_atoms: list 
        List of atom label strings for all atoms in selection 1
    sele2_atoms: list 
        List of atom label strings for all atoms in selection 2
    geom_criteria: dict
        Container for geometric criteria

    Returns
    -------
    t_stacking = list of tuples, [(frame_index, itype, atom1_label, atom2_label), ...]
        itype = "ts"
    """
    cutoff_distance = geom_criteria['T_STACK_CUTOFF_DISTANCE']
    cutoff_angle = geom_criteria['T_STACK_CUTOFF_ANGLE']
    psi_angle = geom_criteria['T_STACK_PSI_ANGLE']

    t_stacking = compute_aromatics(traj_frag_molid, frame_idx, index_to_atom,
                                   sele1, sele2,
                                   sele1_atoms, sele2_atoms, "ts",
                                   T_STACK_SOFT_DISTANCE_CUTOFF, cutoff_distance,
                                   cutoff_angle, psi_angle)
    return t_stacking


