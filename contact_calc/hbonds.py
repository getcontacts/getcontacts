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

from collections import defaultdict
from itertools import product
from vmd import *
from .contact_utils import *

__all__ = ['compute_hydrogen_bonds']

WATER_SHELL_RAD = 3.5

int_pattern = re.compile(r'\d+')


def extract_donor_acceptor(hbond_string):
    """
    Extract donors and acceptors from a vmd hbond output

    Parameters
    ----------
    hbond_string: str
        The output of vmds `measure hbonds` command; three lists indicating indices of donors, acceptors and hydrogens
        respectively. For example: "{106 91} {91 55} {107 92}"

    Returns
    -------
    set of (int, int)
        Set of donor and acceptor indices. For example: `set([(106, 91), (91, 55)])`
    """
    atom_indices = [int(index) for index in int_pattern.findall(hbond_string)]
    third = len(atom_indices) // 3
    return set(zip(atom_indices[0:third], atom_indices[third:third*2]))


def filter_dual_selection(sele1_atoms, sele2_atoms, idx1, idx2):
    """
    Filter interactions between selection 1 and selection 2

    Parameters
    ----------
    sele1_atoms: list
        List of atom label strings for all atoms in selection 1
    sele2_atoms: list
        List of atom label strings for all atoms in selection 2
    idx1: int
        Atom index for cation
    idx2: int
        Atom index for aromatic atom
    Returns
    -------
    bool
        True if interaction should be included
    """
    return ((idx1 in sele1_atoms) and (idx2 in sele2_atoms)) or ((idx1 in sele2_atoms) and (idx2 in sele1_atoms))


def compute_hydrogen_bonds(molid, frame, index_to_atom, solvent_resn, ligand_resn, sele1, sele2,
                           sele1_atoms, sele2_atoms, geom_criteria):
    """
    Compute hydrogen bonds involving protein for a single frame of simulation

    Parameters
    ----------
    molid: int
        Specifies which trajectory fragment in VMD to perform computations upon
    frame: int
        Specify frame index with respect to the smaller trajectory fragment
    index_to_atom: dict
        Maps VMD atom index to Atom
    solvent_resn: string
        Denotes the resname of solvent in simulation
    ligand_resn: string
        Denotes the resname(s) of ligand(s) in simulation
    sele1: string, default = None
        Compute contacts on subset of atom selection based on VMD query
    sele2: string, default = None
        If second VMD query is specified, then compute contacts between atom selection 1 and 2
    sele1_atoms: set
        Indices of atoms in selection 1
    sele2_atoms: set
        Indices of atoms in selection 2
    geom_criteria: dict
        Container for geometric criteria

    Return
    ------
    list of tuples, [(frame_idx, atom1_label, atom2_label, itype), ...]
    """
    cutoff_distance = geom_criteria['HBOND_CUTOFF_DISTANCE']
    cutoff_angle = geom_criteria['HBOND_CUTOFF_ANGLE']
    res_diff = geom_criteria['HBOND_RES_DIFF']

    if sele1 == sele2:
        evaltcl("set selunion [atomselect %s \"%s and not (carbon or sulfur)\" frame %s]" % (molid, sele1, frame))
        evaltcl("set shell [atomselect %s \"(solv) and within %s of (%s)\" frame %s]" %
                (molid, WATER_SHELL_RAD, sele1, frame))
    else:
        evaltcl("set selunion [atomselect %s \"((%s) or (%s)) and not (carbon or sulfur)\" frame %s]" % (molid, sele1, sele2, frame))
        evaltcl("set shell [atomselect %s \"(solv) and within %s of ((%s) or (%s))\" frame %s]" %
                (molid, WATER_SHELL_RAD, sele1, sele2, frame))

    sel_sel = extract_donor_acceptor(evaltcl("measure hbonds %s %s $selunion" % (cutoff_distance, cutoff_angle)))
    sel_sel = [(d, a) for (d, a) in sel_sel if filter_dual_selection(sele1_atoms, sele2_atoms, d, a)]

    sel_solv = set([])
    sel_solv |= extract_donor_acceptor(evaltcl("measure hbonds %s %s $selunion $shell" % (cutoff_distance, cutoff_angle)))
    sel_solv |= extract_donor_acceptor(evaltcl("measure hbonds %s %s $shell $selunion" % (cutoff_distance, cutoff_angle)))

    solv_solv = extract_donor_acceptor(evaltcl("measure hbonds %s %s $shell" % (cutoff_distance, cutoff_angle)))
    evaltcl("$selunion delete")
    evaltcl("$shell delete")

    hbonds = []
    # Stratify hbonds within sele1 and sele2
    for d_idx, a_idx in sel_sel:
        d_atom = index_to_atom[d_idx]
        a_atom = index_to_atom[a_idx]

        # Filter away local interactions
        if d_atom.chain == a_atom.chain and abs(d_atom.resid - a_atom.resid) < res_diff:
            continue

        d_bb = d_atom.is_bb()
        a_bb = a_atom.is_bb()
        d_lig = d_atom.resname in ligand_resn
        a_lig = a_atom.resname in ligand_resn

        if d_lig and a_lig:
            hb_type = "hbll"
        elif d_lig:
            if a_bb:
                hb_type = "hblb"
            else:
                hb_type = "hbls"
        elif a_lig:
            if d_bb:
                hb_type = "hblb"
            else:
                hb_type = "hbls"
        elif d_bb and a_bb:
            hb_type = "hbbb"
        elif not d_bb and not a_bb:
            hb_type = "hbss"
        else:
            hb_type = "hbsb"

        hbonds.append([frame, hb_type, d_atom.get_label(), a_atom.get_label()])

    # Build dictionary for water bridges
    water_dict = defaultdict(set)  # Associates water ids to lists of neighbors
    for d_idx, a_idx in sel_solv | solv_solv:
        d_atom = index_to_atom[d_idx]
        a_atom = index_to_atom[a_idx]

        # Filter away local interactions
        if d_atom.chain == a_atom.chain and d_atom.resid == a_atom.resid:
            continue

        d_solv = d_atom.resname in solvent_resn
        a_solv = a_atom.resname in solvent_resn

        if d_solv:
            water_dict[d_idx].add(a_idx)

        if a_solv:
            water_dict[a_idx].add(d_idx)

    for w_idx in water_dict:
        w_atom = index_to_atom[w_idx]
        for n1, n2 in product(water_dict[w_idx], repeat=2):
            n1_atom = index_to_atom[n1]
            n2_atom = index_to_atom[n2]

            n1_solv = n1_atom.resname in solvent_resn
            n2_solv = n2_atom.resname in solvent_resn

            if n1_solv and n2_solv:
                continue

            if n1_solv:
                n1, n2 = n2, n1
                n1_solv, n2_solv = n2_solv, n1_solv
                n1_atom, n2_atom = n2_atom, n1_atom

            n1_lig = n1_atom.resname in ligand_resn

            if n2_solv:
                # Check for water bridges between n1 and any neighbor of n2
                for n2_n in water_dict[n2]:
                    n2_n_atom = index_to_atom[n2_n]
                    if n2_n_atom.resname not in solvent_resn:
                        # The neighbor of n2 is a water, so check for potential extended water bridges
                        if not filter_dual_selection(sele1_atoms, sele2_atoms, n1, n2_n):
                            continue
                        if n1_atom.chain == n2_n_atom.chain and abs(n1_atom.resid - n2_n_atom.resid) < res_diff:
                            continue

                        n2_n_lig = n2_n_atom.resname in ligand_resn
                        if n1_lig or n2_n_lig:
                            hb_type = "lwb2"
                        else:
                            hb_type = "wb2"

                        hbonds.append([frame, hb_type, n1_atom.get_label(), n2_n_atom.get_label(), w_atom.get_label(), n2_atom.get_label()])
            else:
                if not filter_dual_selection(sele1_atoms, sele2_atoms, n1, n2):
                    continue
                if n1_atom.chain == n2_atom.chain and abs(n1_atom.resid - n2_atom.resid) < res_diff:
                    continue

                n2_lig = n2_atom.resname in ligand_resn

                if n1_lig or n2_lig:
                    hb_type = "lwb"
                else:
                    hb_type = "wb"

                hbonds.append([frame, hb_type, n1_atom.get_label(), n2_atom.get_label(), w_atom.get_label()])

    return hbonds
