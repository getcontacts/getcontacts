##############################################################################
# MDContactNetworks: A Python Library for computing non-covalent contacts
#                    throughout Molecular Dynamics Trajectories. 
# Copyright 2016-2017 Stanford University and the Authors
#
# Authors: Anthony Kai Kwang Ma
# Email: anthony.ma@yale.edu, anthonyma27@gmail.com, akma327@stanford.edu
##############################################################################

##############################################################################
# Imports
##############################################################################

from vmd import *
from contact_utils import *

__all__ = ['compute_pi_stacking', 'compute_t_stacking']

##############################################################################
# Functions
##############################################################################


def get_aromatic_triplet(traj_frag_molid, frame_idx, aromatic_residue_label):
    """
    Given an aromatic residue label return triplet of atoms on the ring

    Parameters
    ----------
    traj_frag_molid: int
        Identifier to simulation fragment in VMD
    frame_idx: int
        Frame number to query
    aromatic_residue_label: string
        ie "A:PHE:222"

    Returns
    -------
    aromatic_atom_triplet: list of strings
        ie ["A:PHE:329:CE1:55228", "A:PHE:329:CE2:55234", "A:PHE:329:CG:55225"]
    """


    residue_to_atom_names = {"PHE": "CG CE1 CE2", "TRP": "CD2 CZ2 CZ3", "TYR": "CG CE1 CE2"}
    chain, resname, resid = aromatic_residue_label.split(":")
    evaltcl("set aromatic_atoms [atomselect %s \" (chain %s) and (resname %s) and (resid %s) and (name %s)\" frame %s]" %(traj_frag_molid, chain, resname, resid, residue_to_atom_names[resname], frame_idx))
    aromatic_atom_triplet = get_atom_selection_labels("aromatic_atoms")
    evaltcl('$aromatic_atoms delete')
    return aromatic_atom_triplet


def compute_aromatics(traj_frag_molid, frame_idx, index_to_label, sele_id, itype, SOFT_DISTANCE_CUTOFF, DISTANCE_CUTOFF, ANGLE_CUTOFF, PSI_ANGLE_CUTOFF):
    """
    Compute aromatic interactions in a frame of simulation

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
    itype: string
        Specify which type of aromatics ("ps" or "ts") to compute
    SOFT_DISTANCE_CUTOFF: float
        Soft distance cutoff to find candidate aromatic pairs
    DISTANCE_CUTOFF: float
        Cutoff distance between aromatic centers
    ANGLE_CUTOFF: float
        Cutoff angle for the angle between normal vectors of aromatic planes
    PSI_ANGLE_CUTOFF: float
        Cutoff angle for how aligned two aromatic planes are

    Returns
    -------
    aromatics = list of tuples, [(frame_index, atom1_label, atom2_label, itype), ...]
        itype = "ps" or "ts"
    """

    aromatics = []

    if(sele_id == None):
        aromatic_atom_sel = "set aromatic_atoms [atomselect %s \" ((resname PHE) and (name CG CE1 CE2)) or ((resname TRP) and (name CD2 CZ2 CZ3)) or ((resname TYR) and (name CG CE1 CE2)) \" frame %s]" % (traj_frag_molid, frame_idx)
    else:
        aromatic_atom_sel = "set aromatic_atoms [atomselect %s \" ((resname PHE) and (name CG CE1 CE2) and (%s)) or ((resname TRP) and (name CD2 CZ2 CZ3) and (%s)) or ((resname TYR) and (name CG CE1 CE2) and (%s)) \" frame %s]" % (traj_frag_molid, sele_id, sele_id, sele_id, frame_idx)

    evaltcl(aromatic_atom_sel)
    contacts = evaltcl("measure contacts %s $aromatic_atoms" % (SOFT_DISTANCE_CUTOFF))
    evaltcl("$aromatic_atoms delete")

    ### Calculate set of distinct aromatic candidate pairs that may have pi-stacking
    contact_index_pairs = parse_contacts(contacts)
    res_pairs = set()
    residue_to_atom_labels = {}
    for aromatic1_index, aromatic2_index in contact_index_pairs:
        aromatic1_label = index_to_label[aromatic1_index]
        aromatic2_label = index_to_label[aromatic2_index]

        ### Check if the two atoms belong to same aromatic group
        aromatic1_res = ":".join(aromatic1_label.split(":")[0:3])
        aromatic2_res = ":".join(aromatic2_label.split(":")[0:3])
        if(aromatic1_res == aromatic2_res): continue

        if(aromatic1_res not in residue_to_atom_labels):
            residue_to_atom_labels[aromatic1_res] = get_aromatic_triplet(traj_frag_molid, frame_idx, aromatic1_res)
        if(aromatic2_res not in residue_to_atom_labels):
            residue_to_atom_labels[aromatic2_res] = get_aromatic_triplet(traj_frag_molid, frame_idx, aromatic2_res)

        k1 = (aromatic1_res, aromatic2_res)
        k2 = (aromatic2_res, aromatic1_res)
        if(k1 not in res_pairs and k2 not in res_pairs):
            res_pairs.add((aromatic1_res, aromatic2_res))

    ### Perform strict geometric criterion on candidate aromatic pairs
    for aromatic1_res, aromatic2_res in res_pairs:
        aromatic1_atom_labels = residue_to_atom_labels[aromatic1_res]
        aromatic2_atom_labels = residue_to_atom_labels[aromatic2_res]

        ### Distance between two aromatic centers must be below DISTANCE_CUTOFF
        arom1_atom1_coord = get_coord(traj_frag_molid, frame_idx, aromatic1_atom_labels[0])
        arom1_atom2_coord = get_coord(traj_frag_molid, frame_idx, aromatic1_atom_labels[1])
        arom1_atom3_coord = get_coord(traj_frag_molid, frame_idx, aromatic1_atom_labels[2])

        arom2_atom1_coord = get_coord(traj_frag_molid, frame_idx, aromatic2_atom_labels[0])
        arom2_atom2_coord = get_coord(traj_frag_molid, frame_idx, aromatic2_atom_labels[1])
        arom2_atom3_coord = get_coord(traj_frag_molid, frame_idx, aromatic2_atom_labels[2])

        aromatic1_centroid = calc_geom_centroid(arom1_atom1_coord, arom1_atom2_coord, arom1_atom3_coord)
        aromatic2_centroid = calc_geom_centroid(arom2_atom1_coord, arom2_atom2_coord, arom2_atom3_coord)
        aromatic_centers_distance = calc_geom_distance(aromatic1_centroid, aromatic2_centroid)
        if(aromatic_centers_distance > DISTANCE_CUTOFF): continue

        ### Angle between vectors normal to each aromatic plane must be below cutoff
        aromatic1_normal_vector = calc_geom_normal_vector(arom1_atom1_coord, arom1_atom2_coord, arom1_atom3_coord)
        aromatic2_normal_vector = calc_geom_normal_vector(arom2_atom1_coord, arom2_atom2_coord, arom2_atom3_coord)
        aromatic_plane_alignment_angle = calc_angle_between_vectors(aromatic1_normal_vector, aromatic2_normal_vector)
        if(itype == "ps"):
            aromatic_plane_alignment_angle = min(math.fabs(aromatic_plane_alignment_angle - 0), math.fabs(aromatic_plane_alignment_angle - 180))
            if(aromatic_plane_alignment_angle > ANGLE_CUTOFF): continue
        elif(itype == "ts"):
            aromatic_plane_perpendicular_angle = math.fabs(calc_angle_between_vectors(aromatic1_normal_vector, aromatic2_normal_vector) - 90)
            if(aromatic_plane_perpendicular_angle > ANGLE_CUTOFF): continue

        ### Psi Angle cutoff
        psi_angle1 = calc_geom_psi_angle(aromatic1_centroid, aromatic2_centroid, aromatic1_normal_vector)
        psi_angle2 = calc_geom_psi_angle(aromatic2_centroid, aromatic1_centroid, aromatic2_normal_vector)
        psi_angle = min(psi_angle1, psi_angle2)
        if(psi_angle > PSI_ANGLE_CUTOFF): continue

        ### Return aromatic interaction atom pairs
        for arom1_atom_label in aromatic1_atom_labels:
            for arom2_atom_label in aromatic2_atom_labels:
                aromatics.append([frame_idx, arom1_atom_label, arom2_atom_label, itype])

    return aromatics


def compute_pi_stacking(traj_frag_molid, frame_idx, index_to_label, sele_id, PI_STACK_CUTOFF_DISTANCE=7.0, PI_STACK_CUTOFF_ANGLE=30, PI_STACK_PSI_ANGLE=45):
    """
    Compute pi-stacking interactions in a frame of simulation

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
    PI_STACK_CUTOFF_DISTANCE: float, default = 7.0 angstroms
        cutoff for distance between centroids of two aromatic rings
    PI_STACK_CUTOFF_ANGLE: float, default = 30 degrees
        cutoff for angle between the normal vectors projecting
        from each aromatic plane.
    PI_STACK_PSI_ANGLE: float, default = 45 degrees
        cutoff for angle between normal vector projecting from
        aromatic plane 1 and vector between the two aromatic centroids

    Returns
    -------
    pi_stacking = list of tuples, [(frame_index, atom1_label, atom2_label, itype), ...]
        itype = "ps"
    """

    PI_STACK_SOFT_DISTANCE_CUTOFF = 10.0 # angstroms
    pi_stacking = compute_aromatics(traj_frag_molid, frame_idx, index_to_label, sele_id, "ps", PI_STACK_SOFT_DISTANCE_CUTOFF, PI_STACK_CUTOFF_DISTANCE, PI_STACK_CUTOFF_ANGLE, PI_STACK_PSI_ANGLE)
    return pi_stacking


def compute_t_stacking(traj_frag_molid, frame_idx, index_to_label, sele_id, T_STACK_CUTOFF_DISTANCE=5.0, T_STACK_CUTOFF_ANGLE=30, T_STACK_PSI_ANGLE=45):
    """
    Compute t-stacking interactions in a frame of simulation

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
    T_STACK_CUTOFF_DISTANCE: float, default = 5.0 angstroms
        cutoff for distance between centroids of two aromatic rings
    T_STACK_CUTOFF_ANGLE: float, default = 30 degrees
        cutoff for angle between the normal vectors projecting
        from each aromatic plane minus 90 degrees
    T_STACK_PSI_ANGLE: float, default = 45 degrees
        cutoff for angle between normal vector projecting from
        aromatic plane 1 and vector between the two aromatic
        centroids

    Returns
    -------
    t_stacking = list of tuples, [(frame_index, atom1_label, atom2_label, itype), ...]
        itype = "ts"
    """

    T_STACK_SOFT_DISTANCE_CUTOFF = 6.0 # angstroms
    t_stacking = compute_aromatics(traj_frag_molid, frame_idx, index_to_label, sele_id, "ts", T_STACK_SOFT_DISTANCE_CUTOFF, T_STACK_CUTOFF_DISTANCE, T_STACK_CUTOFF_ANGLE, T_STACK_PSI_ANGLE)
    return t_stacking


