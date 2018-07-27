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

__all__ = ['compute_pi_cation']

SOFT_DISTANCE_CUTOFF = 10.0  # Angstroms

basic_his = "((resname HIS HSD HSE HSP HIE HIP HID) and (name ND1 NE2))"
basic_lys = "((resname LYS) and (name NZ))"
basic_arg = "((resname ARG) and (name NH1 NH2))"
aromatic_phe = "((resname PHE) and (name CG CE1 CE2))"
aromatic_trp = "((resname TRP) and (name CD2 CZ2 CZ3))"
aromatic_tyr = "((resname TYR) and (name CG CE1 CE2))"


def compute_pi_cation(traj_frag_molid, frame, index_to_atom, sele1, sele2, geom_criteria):
    """
    Compute pi-cation interactions in a frame of simulation

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
    sele1_atoms: list 
        List of atom label strings for all atoms in selection 1
    sele2_atoms: list 
        List of atom label strings for all atoms in selection 2
    geom_criteria: dict
        Container for geometric criteria

    Returns
    -------
    pi_cations = list of tuples, [(frame_index, itype, atom1_label, atom2_label), ...]
        itype = "pc"
    """


    # if sele_id is None and sele_id2 is None:
    #     cation_atom_sel = "set cation_atoms [atomselect %s \" " \
    #                       "((resname LYS) and (name NZ)) or " \
    #                       "((resname ARG) and (name NH1 NH2)) or " \
    #                       "((resname HIS HSD HSE HSP HIE HIP HID) and (name ND1 NE2)) \" frame %s]" % \
    #                       (traj_frag_molid, frame_idx)
    #     aromatic_atom_sel = "set aromatic_atoms [atomselect %s \" " \
    #                         "((resname PHE) and (name CG CE1 CE2)) or " \
    #                         "((resname TRP) and (name CD2 CZ2 CZ3)) or " \
    #                         "((resname TYR) and (name CG CE1 CE2)) \" frame %s]" % (traj_frag_molid, frame_idx)
    # elif sele_id is not None and sele_id2 is None:
    #     cation_atom_sel = "set cation_atoms [atomselect %s \" " \
    #                       "((resname LYS) and (name NZ) and (%s)) or " \
    #                       "((resname ARG) and (name NH1 NH2) and (%s)) or " \
    #                       "((resname HIS HSD HSE HSP HIE HIP HID) and (name ND1 NE2) and (%s))\" frame %s]" % \
    #                       (traj_frag_molid, sele_id, sele_id, sele_id, frame_idx)
    #     aromatic_atom_sel = "set aromatic_atoms [atomselect %s \" " \
    #                         "((resname PHE) and (name CG CE1 CE2) and (%s)) or " \
    #                         "((resname TRP) and (name CD2 CZ2 CZ3) and (%s)) or " \
    #                         "((resname TYR) and (name CG CE1 CE2) and (%s))\" frame %s]" % \
    #                         (traj_frag_molid, sele_id, sele_id, sele_id, frame_idx)
    # elif sele_id is not None and sele_id2 is not None:
    #     sele_union = "(%s) or (%s)" % (sele_id, sele_id2)
    #     cation_atom_sel = "set cation_atoms [atomselect %s \" " \
    #                       "((resname LYS) and (name NZ) and (%s)) or " \
    #                       "((resname ARG) and (name NH1 NH2) and (%s)) or " \
    #                       "((resname HIS HSD HSE HSP HIE HIP HID) and (name ND1 NE2) and (%s))\" frame %s]" % \
    #                       (traj_frag_molid, sele_union, sele_union, sele_union, frame_idx)
    #     aromatic_atom_sel = "set aromatic_atoms [atomselect %s \" " \
    #                         "((resname PHE) and (name CG CE1 CE2) and (%s)) or " \
    #                         "((resname TRP) and (name CD2 CZ2 CZ3) and (%s)) or " \
    #                         "((resname TYR) and (name CG CE1 CE2) and (%s))\" frame %s]" % \
    #                         (traj_frag_molid, sele_union, sele_union, sele_union, frame_idx)
    #
    # evaltcl(cation_atom_sel)
    # evaltcl(aromatic_atom_sel)
    # contacts = evaltcl("measure contacts %s $cation_atoas $aromatic_atoms" % SOFT_DISTANCE_CUTOFF)
    # evaltcl("$cation_atoms delete")
    # evaltcl("$aromatic_atoms delete")
    #
    # # Evaluate geometric criterion if all three points of an aromatic
    # # residue are sufficiently close to a cation atom
    # contact_index_pairs = parse_contacts(contacts)

    PI_CATION_CUTOFF_DISTANCE = geom_criteria['PI_CATION_CUTOFF_DISTANCE']
    PI_CATION_CUTOFF_ANGLE = geom_criteria['PI_CATION_CUTOFF_ANGLE']

    s1_aroms = "(%s or %s or %s) and %s" % (aromatic_phe, aromatic_trp, aromatic_tyr, sele1)
    s2_aroms = "(%s or %s or %s) and %s" % (aromatic_phe, aromatic_trp, aromatic_tyr, sele2)
    s1_cations = "( %s or %s or %s) and %s" % (basic_his, basic_lys, basic_arg, sele1)
    s2_cations = "( %s or %s or %s) and %s" % (basic_his, basic_lys, basic_arg, sele2)

    evaltcl("set s1aroms [atomselect %s \" %s \" frame %s]" % (traj_frag_molid, s1_aroms, frame))
    evaltcl("set s2aroms [atomselect %s \" %s \" frame %s]" % (traj_frag_molid, s2_aroms, frame))
    evaltcl("set s1cations [atomselect %s \" %s \" frame %s]" % (traj_frag_molid, s1_cations, frame))
    evaltcl("set s2cations [atomselect %s \" %s \" frame %s]" % (traj_frag_molid, s2_cations, frame))
    contacts_12 = set(parse_contacts(evaltcl("measure contacts %f $s1cations $s2aroms" % SOFT_DISTANCE_CUTOFF)))
    contacts_21 = set(parse_contacts(evaltcl("measure contacts %f $s2cations $s1aroms" % SOFT_DISTANCE_CUTOFF)))
    evaltcl("$s1aroms delete")
    evaltcl("$s2aroms delete")
    evaltcl("$s1cations delete")
    evaltcl("$s2cations delete")

    contact_index_pairs = contacts_12 | contacts_21

    # map every distinct combination of cation atom and aromatic residue to the three atoms on the aromatic atom
    pi_cation_aromatic_grouping = {}
    for cation_index, aromatic_index in contact_index_pairs:
        # if sele1_atoms is not None and sele2_atoms is not None:
        #     if filter_dual_selection_pi_cation(sele1_atoms, sele2_atoms, cation_index, aromatic_index):
        #         continue

        cation_label = index_to_atom[cation_index].get_label()
        aromatic_label = index_to_atom[aromatic_index].get_label()
        pi_cation_aromatic_res_key = cation_label + ":" + ":".join(aromatic_label.split(":")[0:3])
        if pi_cation_aromatic_res_key not in pi_cation_aromatic_grouping:
            pi_cation_aromatic_grouping[pi_cation_aromatic_res_key] = set()
        pi_cation_aromatic_grouping[pi_cation_aromatic_res_key].add(aromatic_label)

    # Apply strict geometric criterion
    pi_cations = []
    for pi_cation_aromatic_res_key in pi_cation_aromatic_grouping:
        cation_atom_label = ":".join(pi_cation_aromatic_res_key.split(":")[0:5])
        aromatic_atom_labels = pi_cation_aromatic_grouping[pi_cation_aromatic_res_key]
        if len(aromatic_atom_labels) != 3:
            continue
        aromatic_atom_labels = sorted(list(aromatic_atom_labels))
        arom_atom1_label, arom_atom2_label, arom_atom3_label = aromatic_atom_labels

        # Compute coordinates of cation and aromatic atoms
        cation_coord = get_coord(traj_frag_molid, frame, cation_atom_label)
        arom_atom1_coord = get_coord(traj_frag_molid, frame, arom_atom1_label)
        arom_atom2_coord = get_coord(traj_frag_molid, frame, arom_atom2_label)
        arom_atom3_coord = get_coord(traj_frag_molid, frame, arom_atom3_label)

        # Perform distance criterion
        aromatic_centroid = calc_geom_centroid(arom_atom1_coord, arom_atom2_coord, arom_atom3_coord)
        cation_to_centroid_distance = calc_geom_distance(cation_coord, aromatic_centroid)
        if cation_to_centroid_distance > PI_CATION_CUTOFF_DISTANCE:
            continue

        # Perform angle criterion
        aromatic_plane_norm_vec = calc_geom_normal_vector(arom_atom1_coord, arom_atom2_coord, arom_atom3_coord)
        aromatic_center_to_cation_vec = points_to_vector(aromatic_centroid, cation_coord)
        cation_norm_offset_angle = calc_angle_between_vectors(aromatic_plane_norm_vec, aromatic_center_to_cation_vec)
        cation_norm_offset_angle = min(math.fabs(cation_norm_offset_angle - 0), math.fabs(cation_norm_offset_angle - 180))
        if cation_norm_offset_angle > PI_CATION_CUTOFF_ANGLE:
            continue

        # Append just the CG atom of the aromatic ring
        single_arom_atom_label = convert_to_single_atom_aromatic_string(arom_atom1_label)
        pi_cations.append([frame, "pc", cation_atom_label, single_arom_atom_label])

        # Append three of the aromatic atoms
        # pi_cations.append([frame_idx, "pc", cation_atom_label, arom_atom1_label])
        # pi_cations.append([frame_idx, "pc", cation_atom_label, arom_atom2_label])
        # pi_cations.append([frame_idx, "pc", cation_atom_label, arom_atom3_label])

    return pi_cations
