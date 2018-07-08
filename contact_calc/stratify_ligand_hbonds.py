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


def stratify_ligand_residue_hbonds(ligand_residue_hbonds, ligands):
    """
    Stratify ligand to residue hbonds into those involving sidechain
    or backbone atoms.
    """
    backbone_atoms = ['N', 'O']
    hbls, hblb = [], []

    # Iterate through each ligand residue hbond and bin into appropriate subtype
    for frame_idx, atom1_label, atom2_label, itype in ligand_residue_hbonds:
        a1_is_ligand = any([l in atom1_label for l in ligands])
        a2_is_ligand = any([l in atom2_label for l in ligands])

        if a1_is_ligand and a2_is_ligand:
            continue
        elif not a1_is_ligand and not a2_is_ligand:
            continue
        elif a1_is_ligand:
            lig_atom = atom1_label
            res_atom = atom2_label
        elif a2_is_ligand:
            lig_atom = atom2_label
            res_atom = atom1_label

        res_atom_name = res_atom.split(":")[3]
        if res_atom_name in backbone_atoms:
            hblb.append([frame_idx, "hblb", lig_atom, res_atom])
        else:
            hbls.append([frame_idx, "hbls", lig_atom, res_atom])

    return hbls, hblb


def stratify_ligand_vs_protein(atom_list, ligands):
    ligand_atoms, protein_atoms = [], []
    for atom in atom_list:
        atom_is_ligand = any([l in atom for l in ligands])
        if atom_is_ligand:
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
                ligand_water_bridges.add((frame_idx, "lwb", lig_atom, res_atom, water))

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
                extended_ligand_water_bridges.add((frame_idx, "lwb2", lig_atom1, res_atom2, water1, water2))

        for lig_atom2 in ligand_atoms2:
            for res_atom1 in protein_atoms1:
                extended_ligand_water_bridges.add((frame_idx, "lwb2", lig_atom2, res_atom1, water2, water1))

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
    ligand: list of string
        Denotes the resname of ligand molecule

    Returns
    -------
    hbond_subtypes: list, [[frame_idx, atom1_label, atom2_label, itype], ...]
        List of all hydrogen contacts with itype = "hls", "hlb", "lwb", or "lwb2"
        corresponding to ligand-sidechain, ligand-backbone, ligand water bridge,
        and extended ligand water bridge respectively.
    """
    ligand_residue_hbonds, water_hbonds = ligand_residue_vs_water_hbonds(hbonds, solvent_resn, ligand)
    hbls, hblb = stratify_ligand_residue_hbonds(ligand_residue_hbonds, ligand)
    lwb = stratify_ligand_water_bridge(water_hbonds, solvent_resn, ligand)
    lwb2 = stratify_extended_ligand_water_bridge(water_hbonds, solvent_resn, ligand)
    hbonds = hbls + hblb + lwb + lwb2

    return hbonds
