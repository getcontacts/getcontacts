
from vmd import *
from .contact_utils import *

__all__ = ['compute_hydrophobics']


def compute_hydrophobics(traj_frag_molid, frame_idx, index_to_atom, sele1, sele2, geom_criteria):
    """
    Compute hydrophobic interactions in a frame of simulation

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
    geom_criteria: dict
        Container for geometric criteria

    Returns
    -------
    list of tuples, [(frame_index, atom1_label, atom2_label, itype), ...]
        itype = "hp"
    """
    epsilon = geom_criteria['VDW_EPSILON']
    res_diff = geom_criteria['VDW_RES_DIFF']

    aa_sel = "(resname ALA CYS PHE GLY ILE LEU MET PRO VAL TRP or ligand) and (element C S)"

    if sele1 == sele2:
        evaltcl("set hp_atoms [atomselect %s \"(%s) and (%s)\" frame %s]" % (traj_frag_molid, aa_sel, sele1, frame_idx))
        contacts = evaltcl("measure contacts %s $hp_atoms" % (epsilon + 2 * 1.8))
        evaltcl("$hp_atoms delete")
    else:
        evaltcl("set hp_atoms1 [atomselect %s \"(%s) and (%s)\" frame %s]" % (traj_frag_molid, aa_sel, sele1, frame_idx))
        evaltcl("set hp_atoms2 [atomselect %s \"(%s) and (%s)\" frame %s]" % (traj_frag_molid, aa_sel, sele2, frame_idx))
        contacts = evaltcl("measure contacts %s $hp_atoms1 $hp_atoms2" % (epsilon + 2 * 1.8))
        evaltcl("$hp_atoms1 delete")
        evaltcl("$hp_atoms2 delete")

    ret = []
    contact_index_pairs = parse_contacts(contacts)
    for atom1_index, atom2_index in contact_index_pairs:
        # Convert to atom label
        atom1, atom2 = index_to_atom[atom1_index], index_to_atom[atom2_index]

        if atom1.chain == atom2.chain and abs(atom1.resid - atom2.resid) < res_diff:
            continue

        # Perform distance cutoff with atom indices
        distance = compute_distance(traj_frag_molid, frame_idx, atom1_index, atom2_index)
        vanderwaal_cutoff = atom1.vdwradius + atom2.vdwradius + epsilon
        if distance < vanderwaal_cutoff:
            ret.append([frame_idx, "hp", atom1.get_label(), atom2.get_label()])

    return ret


__author__ = 'Anthony Ma <anthonyma27@gmail.com>, Rasmus Fonseca <fonseca.rasmus@gmail.com>'
__license__ = "Apache License 2.0"

