##############################################################################
# MDContactNetworks: A Python Library for computing non-covalent contacts
#                    throughout Molecular Dynamics Trajectories. 
#
# Contact: Anthony Kai Kwang Ma, anthonyma27@gmail.com
##############################################################################

############################################################################
# Imports
############################################################################

from .contact_utils import *

__all__ = ["stratify_ligand_hbond_subtypes"]

############################################################################
# Functions
############################################################################


def ligand_residue_vs_water_hbonds(hbonds, solvent_resn, ligand):
    """
    Split hbonds into those involving residue and ligand directly and those
    mediated by water molecules.
    """
    ligand_residue_hbonds, water_hbonds = [], []
    for hbond in hbonds:
        frame_idx, atom1_label, atom2_label, itype = hbond
        if solvent_resn in atom1_label or solvent_resn in atom2_label:
            water_hbonds.append(hbond)
        else:
            ligand_residue_hbonds.append(hbond)

    return ligand_residue_hbonds, water_hbonds


def stratify_ligand_residue_hbonds(ligand_residue_hbonds, ligand):
    """
    Stratify ligand to residue hbonds into those involving sidechain
    or backbone atoms.
    """
    backbone_atoms = ['N', 'O']
    hls, hlb = [], []

    # Iterate through each ligand residue hbond and bin into appropriate subtype
    for frame_idx, atom1_label, atom2_label, itype in ligand_residue_hbonds:
        if ligand in atom1_label and ligand in atom2_label:
            continue
        elif ligand not in atom1_label and ligand not in atom2_label:
            continue
        elif ligand in atom1_label and ligand not in atom2_label:
            lig_atom = atom1_label
            res_atom = atom2_label
        elif ligand not in atom1_label and ligand in atom2_label:
            lig_atom = atom2_label
            res_atom = atom1_label

        res_atom_name = res_atom.split(":")[3]
        if res_atom_name in backbone_atoms:
            hlb.append([frame_idx, lig_atom, res_atom, "hlb"])
        else:
            hls.append([frame_idx, lig_atom, res_atom, "hls"])

    return hls, hlb


def stratify_ligand_vs_protein(atom_list, ligand):
    ligand_atoms, protein_atoms = [], []
    for atom in atom_list:
        if ligand in atom:
            ligand_atoms.append(atom)
        else:
            protein_atoms.append(atom)
    return ligand_atoms, protein_atoms


def stratify_ligand_water_bridge(water_hbonds, solvent_resn, ligand):
    """
    Compute water bridges between ligand and binding pocket residue
    """
    frame_idx, water_to_ligand_residues, _ = calc_water_to_residues_map(water_hbonds, solvent_resn)
    ligand_water_bridges = set()

    for water in water_to_ligand_residues:
        ligand_atoms, protein_atoms = stratify_ligand_vs_protein(water_to_ligand_residues[water], ligand)

        # Form ligand -- water -- protein pairs
        for lig_atom in ligand_atoms:
            for res_atom in protein_atoms:
                ligand_water_bridges.add((frame_idx, lig_atom, water, res_atom, "lwb"))

    lwb = sorted([list(entry) for entry in ligand_water_bridges])
    return lwb


def stratify_extended_ligand_water_bridge(water_hbonds, solvent_resn, ligand):
    """
    Compute extended water bridges between ligand and binding pocket residues
    """
    frame_idx, water_to_ligand_residues, solvent_bridges = calc_water_to_residues_map(water_hbonds, solvent_resn)

    extended_ligand_water_bridges = set()
    for water1, water2 in solvent_bridges:
        if water1 not in water_to_ligand_residues or water2 not in water_to_ligand_residues:
            continue
        lig_res_atom1_list, lig_res_atom2_list = water_to_ligand_residues[water1], water_to_ligand_residues[water2]
        ligand_atoms1, protein_atoms1 = stratify_ligand_vs_protein(lig_res_atom1_list, ligand)
        ligand_atoms2, protein_atoms2 = stratify_ligand_vs_protein(lig_res_atom2_list, ligand)

        for lig_atom1 in ligand_atoms1:
            for res_atom2 in protein_atoms2:
                extended_ligand_water_bridges.add((frame_idx, lig_atom1, water1, water2, res_atom2, "lwb2"))

        for lig_atom2 in ligand_atoms2:
            for res_atom1 in protein_atoms1:
                extended_ligand_water_bridges.add((frame_idx, lig_atom2, water2, water1, res_atom1, "lwb2"))

    lwb2 = sorted([list(entry) for entry in extended_ligand_water_bridges])
    return lwb2


def stratify_ligand_hbond_subtypes(hbonds, solvent_resn, ligand):
    """
    Stratify the full ligand hbonds list into the following subtypes: ligand-sidechain,
    ligand-backbone, ligand water-bridge, and extended ligand water-bridge

    Parameters
    ----------
    hbonds: list, [[frame_idx, atom1_label, atom2_label, itype], ...]
        List of all hydrogen bond contacts in a single frame. itype = "lhb"

    solvent_resn: string, default = TIP3
        Denotes the resname of solvent in simulation
    ligand: string
        Denotes the resname of ligand molecule

    Returns
    -------
    hbond_subtypes: list, [[frame_idx, atom1_label, atom2_label, itype], ...]
        List of all hydrogen contacts with itype = "hls", "hlb", "lwb", or "lwb2"
        corresponding to ligand-sidechain, ligand-backbone, ligand water bridge,
        and extended ligand water bridge respectively.
    """
    ligand_residue_hbonds, water_hbonds = ligand_residue_vs_water_hbonds(hbonds, solvent_resn, ligand)
    hls, hlb = stratify_ligand_residue_hbonds(ligand_residue_hbonds, ligand)
    lwb = stratify_ligand_water_bridge(water_hbonds, solvent_resn, ligand)
    lwb2 = stratify_extended_ligand_water_bridge(water_hbonds, solvent_resn, ligand)
    hbonds = hls + hlb + lwb + lwb2

    return hbonds
