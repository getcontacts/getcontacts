
from vmd import *
from .contact_utils import *

__all__ = ['compute_hydrophobics']

sele1_hp_indices = None
sele2_hp_indices = None


def prepare_indices(molid, index_to_atom, sele1, sele2):
    global sele1_hp_indices, sele2_hp_indices
    aa_sel = "(resname ALA PHE GLY ILE LEU PRO VAL TRP or ligand) and (element C)"

    evaltcl("set allcarbons1 [atomselect %s \"%s and (%s)\" frame %s]" % (molid, aa_sel, sele1, 0))
    all_carbons1 = get_atom_selection_indices("allcarbons1")
    evaltcl("$allcarbons1 delete")
    # print(all_carbons1)

    sele1_hp_indices = []
    for c_idx in all_carbons1:
        evaltcl("set cneighbors [atomselect %s \"within 1.7 of (index %d)\" frame %s]" % (molid, c_idx, 0))
        c_neighbor_indices = get_atom_selection_indices("cneighbors")
        evaltcl("$cneighbors delete")
        # print("neighbors of", c_idx, ":", c_neighbor_indices)

        c_neighbor_elements = set([index_to_atom[n].element for n in c_neighbor_indices])
        if c_neighbor_elements | set(["C", "H"]) == set(["C", "H"]):
            sele1_hp_indices.append(str(c_idx))
    sele1_hp_indices = " ".join(sele1_hp_indices)

    evaltcl("set allcarbons2 [atomselect %s \"%s and (%s)\" frame %s]" % (molid, aa_sel, sele2, 0))
    all_carbons2 = get_atom_selection_indices("allcarbons2")
    evaltcl("$allcarbons2 delete")

    sele2_hp_indices = []
    for c_idx in all_carbons2:
        evaltcl("set cneighbors [atomselect %s \"within 1.7 of (index %d)\" frame %s]" % (molid, c_idx, 0))
        c_neighbor_indices = get_atom_selection_indices("cneighbors")
        evaltcl("$cneighbors delete")

        c_neighbor_elements = set([index_to_atom[n].element for n in c_neighbor_indices])
        if c_neighbor_elements == set(["C", "H"]):
            sele2_hp_indices.append(str(c_idx))
    sele2_hp_indices = " ".join(sele2_hp_indices)

    # print("preparing hp")
    # print("sele1:",sele1_hp_indices)
    # print("sele2:",sele2_hp_indices)


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

    if sele1_hp_indices is None:
        prepare_indices(traj_frag_molid, index_to_atom, sele1, sele2)

    if not sele1_hp_indices:
        contacts = ""
    elif sele1 == sele2:
        evaltcl("set hp_atoms [atomselect %s \"index %s\" frame %s]" % (traj_frag_molid, sele1_hp_indices, frame_idx))
        contacts = evaltcl("measure contacts %s $hp_atoms" % (epsilon + 2 * 1.7))
        evaltcl("$hp_atoms delete")
    else:
        evaltcl("set hp_atoms1 [atomselect %s \"index %s\" frame %s]" % (traj_frag_molid, sele1_hp_indices, frame_idx))
        evaltcl("set hp_atoms2 [atomselect %s \"index %s\" frame %s]" % (traj_frag_molid, sele2_hp_indices, frame_idx))
        # evaltcl("set hp_atoms1 [atomselect %s \"(%s) and (%s)\" frame %s]" % (traj_frag_molid, aa_sel, sele1, frame_idx))
        # evaltcl("set hp_atoms2 [atomselect %s \"(%s) and (%s)\" frame %s]" % (traj_frag_molid, aa_sel, sele2, frame_idx))
        contacts = evaltcl("measure contacts %s $hp_atoms1 $hp_atoms2" % (epsilon + 2 * 1.7))
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

