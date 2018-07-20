
from vmd import *
from .contact_utils import *

__all__ = ['compute_hydrophobics']

ATOM_RADIUS = {'C': 1.70,
               'S': 1.80}


def compute_hydrophobics(traj_frag_molid, frame_idx, index_to_label, sele_id1, sele_id2, ligands,
                         VDW_EPSILON, VDW_RES_DIFF):
    """
    Compute hydrophobic interactions in a frame of simulation

    Parameters
    ----------
    traj_frag_molid: int
        Identifier to simulation fragment in VMD
    frame_idx: int
        Frame number to query
    index_to_label: dict
        Maps VMD atom index to label "chain:resname:resid:name:index"
        {11205: "A:ASP:114:CA:11205, ...}
    sele_id1: string, default = None
        Compute contacts on subset of atom selection based on VMD query
    sele_id2: string, default = None
        If second VMD query is specified, then compute contacts between atom selection 1 and 2 
    ligands: list of string
        Residue names of ligands
    VDW_EPSILON: float, default = 0.5 angstroms
        amount of padding for calculating vanderwaals contacts
    VDW_RES_DIFF: int, default = 2
        minimum residue distance for which to consider computing interactions

    Returns
    -------
    list of tuples, [(frame_index, atom1_label, atom2_label, itype), ...]
        itype = "hp"
    """

    lig_sel = "" if ligands is None else " ".join(ligands)
    aa_sel = "((resname ALA CYS PHE GLY ILE LEU MET PRO VAL TRP " + lig_sel + ") and (element C S))"
    sel1 = "" if sele_id1 is None else " and (" + sele_id1 + ")"
    sel2 = "" if sele_id2 is None else " and (" + sele_id2 + ")"

    evaltcl("set hp_atoms1 [atomselect %s \"%s%s\" frame %s]" % (traj_frag_molid, aa_sel, sel1, frame_idx))
    evaltcl("set hp_atoms2 [atomselect %s \"%s%s\" frame %s]" % (traj_frag_molid, aa_sel, sel2, frame_idx))
    contacts = evaltcl("measure contacts %s $hp_atoms1 $hp_atoms2" % (VDW_EPSILON + 2 * ATOM_RADIUS['S']))
    evaltcl("$hp_atoms1 delete")
    evaltcl("$hp_atoms2 delete")

    ret = []
    contact_index_pairs = parse_contacts(contacts)
    for atom1_index, atom2_index in contact_index_pairs:
        # Convert to atom label
        atom1_label, atom2_label = index_to_label[atom1_index], index_to_label[atom2_index]

        # Another optimization is only convert to labels after passing initial checks.
        atom1_label_split = atom1_label.split(":")
        atom2_label_split = atom2_label.split(":")

        chain1 = atom1_label_split[0]
        chain2 = atom2_label_split[0]
        resi1 = int(atom1_label_split[2])
        resi2 = int(atom2_label_split[2])
        if chain1 == chain2 and abs(resi1-resi2) < VDW_RES_DIFF:
            continue

        element1 = atom1_label_split[3][0]
        element2 = atom2_label_split[3][0]

        if element1 not in ATOM_RADIUS:
            element1 = 'C'
        if element2 not in ATOM_RADIUS:
            element2 = 'C'

        # Perform distance cutoff with atom indices
        distance = compute_distance(traj_frag_molid, frame_idx, atom1_index, atom2_index)
        vanderwaal_cutoff = ATOM_RADIUS[element1] + ATOM_RADIUS[element2] + VDW_EPSILON
        if distance < vanderwaal_cutoff:
            ret.append([frame_idx, "hp", atom1_label, atom2_label])

    return ret


__author__ = 'Anthony Ma <anthonyma27@gmail.com>, Rasmus Fonseca <fonseca.rasmus@gmail.com>'
__license__ = "Apache License 2.0"

