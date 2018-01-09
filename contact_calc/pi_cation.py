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

__all__ = ['compute_pi_cation']

##############################################################################
# Globals
##############################################################################
SOFT_DISTANCE_CUTOFF = 10.0  # Angstroms

##############################################################################
# Functions
##############################################################################


def compute_pi_cation(traj_frag_molid, frame_idx, index_to_label, sele_id,
                      PI_CATION_CUTOFF_DISTANCE=6.0, PI_CATION_CUTOFF_ANGLE=60):
    """
    Compute pi-cation interactions in a frame of simulation

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
    PI_CATION_CUTOFF_DISTANCE: float, default = 6.0 angstroms
        cutoff for distance between cation and centroid of aromatic ring
    PI_CATION_CUTOFF_ANGLE: float, default = 60 degrees
        cutoff for angle between normal vector projecting
        from aromatic plane and vector from aromatic center
        to cation atom


    Returns
    -------
    pi_cations = list of tuples, [(frame_index, itype, atom1_label, atom2_label), ...]
        itype = "pc"
    """
    pi_cations = []

    if sele_id is None:
        cation_atom_sel = "set cation_atoms [atomselect %s \" " \
                          "((resname LYS) and (name NZ)) or " \
                          "((resname ARG) and (name NH1 NH2)) or " \
                          "((resname HIS HSD HSE HSP HIE HIP HID) and (name ND1 NE2)) \" frame %s]" % \
                          (traj_frag_molid, frame_idx)
        aromatic_atom_sel = "set aromatic_atoms [atomselect %s \" " \
                            "((resname PHE) and (name CG CE1 CE2)) or " \
                            "((resname TRP) and (name CD2 CZ2 CZ3)) or " \
                            "((resname TYR) and (name CG CE1 CE2)) \" frame %s]" % (traj_frag_molid, frame_idx)
    else:
        cation_atom_sel = "set cation_atoms [atomselect %s \" " \
                          "((resname LYS) and (name NZ) and (%s)) or " \
                          "((resname ARG) and (name NH1 NH2) and (%s)) or " \
                          "((resname HIS HSD HSE HSP HIE HIP HID) and (name ND1 NE2) and (%s))\" frame %s]" % \
                          (traj_frag_molid, sele_id, sele_id, sele_id, frame_idx)
        aromatic_atom_sel = "set aromatic_atoms [atomselect %s \" " \
                            "((resname PHE) and (name CG CE1 CE2) and (%s)) or " \
                            "((resname TRP) and (name CD2 CZ2 CZ3) and (%s)) or " \
                            "((resname TYR) and (name CG CE1 CE2) and (%s))\" frame %s]" % \
                            (traj_frag_molid, sele_id, sele_id, sele_id, frame_idx)

    evaltcl(cation_atom_sel)
    evaltcl(aromatic_atom_sel)
    contacts = evaltcl("measure contacts %s $cation_atoms $aromatic_atoms" % SOFT_DISTANCE_CUTOFF)
    evaltcl("$cation_atoms delete")
    evaltcl("$aromatic_atoms delete")

    # Evaluate geometric criterion if all three points of an aromatic
    # residue are sufficiently close to a cation atom
    contact_index_pairs = parse_contacts(contacts)
    pi_cation_aromatic_grouping = {}
    for cation_index, aromatic_index in contact_index_pairs:
        cation_label = index_to_label[cation_index]
        aromatic_label = index_to_label[aromatic_index]
        pi_cation_aromatic_res_key = cation_label + ":" + ":".join(aromatic_label.split(":")[2])
        if pi_cation_aromatic_res_key not in pi_cation_aromatic_grouping:
            pi_cation_aromatic_grouping[pi_cation_aromatic_res_key] = set()
        pi_cation_aromatic_grouping[pi_cation_aromatic_res_key].add(aromatic_label)

    # Apply strict geometric criterion
    for pi_cation_aromatic_res_key in pi_cation_aromatic_grouping:
        cation_atom_label = ":".join(pi_cation_aromatic_res_key.split(":")[0:5])
        aromatic_atom_labels = pi_cation_aromatic_grouping[pi_cation_aromatic_res_key]
        if len(aromatic_atom_labels) != 3:
            continue
        aromatic_atom_labels = sorted(list(aromatic_atom_labels))
        arom_atom1_label, arom_atom2_label, arom_atom3_label = aromatic_atom_labels

        # Compute coordinates of cation and aromatic atoms
        cation_coord = get_coord(traj_frag_molid, frame_idx, cation_atom_label)
        arom_atom1_coord = get_coord(traj_frag_molid, frame_idx, arom_atom1_label)
        arom_atom2_coord = get_coord(traj_frag_molid, frame_idx, arom_atom2_label)
        arom_atom3_coord = get_coord(traj_frag_molid, frame_idx, arom_atom3_label)

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

        # Append three of the aromatic atoms
        pi_cations.append([frame_idx, "pc", cation_atom_label, arom_atom1_label])
        pi_cations.append([frame_idx, "pc", cation_atom_label, arom_atom2_label])
        pi_cations.append([frame_idx, "pc", cation_atom_label, arom_atom2_label])

    return pi_cations
