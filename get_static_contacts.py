#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Computes non-covalent contact networks for a protein structure.

The input is a path to a structure and a specification of interaction-types. The
output is a tab-separated file where each line (except the first two) records frame
(always 0), interaction-type, and atoms involved in the interaction, e.g.:

    # total_frames:20000 interaction_types:sb,pc,ps,ts,hb
    # Columns: frame, interaction_type, atom_1, atom_2[, atom_3[, atom_4]]
    0   sb     C:GLU:21:OE2    C:ARG:86:NH2
    0   ps     C:TYR:36:CG     C:TRP:108:CG
    0   ts     A:TYR:36:CG     A:TRP:108:CG
    0   hbss   A:GLN:53:NE2    A:GLN:69:OE1
    0   wb2    C:ASN:110:O     C:SER:111:OG    W:TIP3:1524:OH2    W:TIP3:2626:OH2
    0   hbsb   A:LYS:28:N      A:HIS:27:ND1
    0   hbsb   A:ASP:52:OD2    A:PHE:48:O
    ...

Interactions that involve more than two atoms (water bridges and extended water bridges)
have extra columns to denote the identities of the water molecules. For simplicity, all
stacking and pi-cation interactions involving an aromatic ring will be denoted by the
CG atom. If there are multiple conformers in the structure, all interactions between all
conformers are computed, but the output will not reflect which conformers an interaction
relates to.

Interaction types are denoted by the following abbreviations:
  · sb - salt bridges
  · pc - pi-cation
  · ps - pi-stacking
  · ts - T-stacking
  · vdw - van der Waals
  Hydrogen bond subtypes:
    · hbbb - Backbone-backbone hydrogen bonds
    · hbsb - Backbone-sidechain hydrogen bonds
    · hbss - Sidechain-sidechain hydrogen bonds
    · wb - Water-mediated hydrogen bond
    · wb2 - Extended water-mediated hydrogen bond
  Ligand-hydrogen bond subtypes
    · hlb - Ligand-backbone hydrogen bonds
    · hls - Ligand-sidechain hydrogen bonds
    · lwb - Ligand water-mediated hydrogen bond
    · lwb2 - Ligand extended water-mediated hydrogen bond


Examples:

Compute salt bridges and hydrogen bonds for residues 100 to 160:
    get_static_contacts.py --structure struc.pdb \\
                           --output output.tsv \\
                           --solv IP3 \\
                           --sele "chain A and resid 100 to 160" \\
                           --ligand EJ4 \\
                           --itypes sb hb lhb

Pi-cation, pi-stacking, and vanderwaals contacts in the entire protein:
    get_static_contacts.py --structure.psf \\
                           --output output.tsv \\
                           --itypes pc ps vdw

Salt bridges and hydrogen bonds in the entire protein with modified distance cutoffs:
    get_static_contacts.py --structure struc.pdb \\
                           --output output.tsv \\
                           --sb_cutoff_dist 5.0 \\
                           --hbond_cutoff_dist 4.5 \\
                           --itypes sb hb
"""

import contact_calc.argparsers as ap
from contact_calc.compute_contacts import *


def process_geometric_criterion_args(args):
    geom_criterion_values = {
        "SALT_BRIDGE_CUTOFF_DISTANCE": args.sb_cutoff_dist,
        "PI_CATION_CUTOFF_DISTANCE": args.pc_cutoff_dist,
        "PI_CATION_CUTOFF_ANGLE": args.pc_cutoff_ang,
        "PI_STACK_CUTOFF_DISTANCE": args.ps_cutoff_dist,
        "PI_STACK_CUTOFF_ANGLE": args.ps_cutoff_ang,
        "PI_STACK_PSI_ANGLE": args.ps_psi_ang,
        "T_STACK_CUTOFF_DISTANCE": args.ts_cutoff_dist,
        "T_STACK_CUTOFF_ANGLE": args.ts_cutoff_ang,
        "T_STACK_PSI_ANGLE": args.ts_psi_ang,
        "HBOND_CUTOFF_DISTANCE": args.hbond_cutoff_dist,
        "HBOND_CUTOFF_ANGLE": args.hbond_cutoff_ang,
        "HBOND_RES_DIFF": args.hbond_res_diff,
        "VDW_EPSILON": args.vdw_epsilon,
        "VDW_RES_DIFF": args.vdw_res_diff
    }
    return geom_criterion_values


def main(argv=None):
    # Parse arguments
    parser = ap.build_getcontact_parser(__doc__, False)
    args, unknown = parser.parse_known_args(argv)

    top = args.structure
    traj = args.structure
    output = args.output
    cores = 1
    ligand = args.ligand
    solv = args.solv
    sele = args.sele
    sele2 = args.sele2
    beg = 0
    end = 0
    stride = 1
    geom_criterion_values = process_geometric_criterion_args(args)

    # If user mistakenly includes only sele2
    if(sele == None and sele2 != None):
        sele, sele2 = sele2, sele

    # Check interaction types
    all_itypes = ["sb", "pc", "ps", "ts", "vdw", "hb", "lhb"]
    if "all" in args.itypes:
        itypes = all_itypes
    else:
        for itype in args.itypes:
            if itype not in all_itypes:
                parser.error(itype + " is not a valid interaction type")

        itypes = args.itypes

    # Begin computation
    tic = datetime.datetime.now()
    compute_contacts(top, traj, output, itypes, geom_criterion_values, cores, beg, end, stride, solv, sele, sele2, ligand)
    toc = datetime.datetime.now()
    print("\nTotal computation time: " + str((toc-tic).total_seconds()) + " seconds")

    print("topology=%s" % top)
    print("trajectory=%s" % traj)
    print("output=%s" % output)
    print("cores=%s" % cores)
    print("ligand=%s" % ",".join(ligand))
    print("solv=%s" % solv)
    print("sele=%s" % sele)
    print("sele2=%s" % sele2)
    print("stride=%s" % stride)


if __name__ == "__main__":
    main()

    # Suppress stdout from vmd as program terminates
    devnull = open("/dev/null", "w")
    os.dup2(devnull.fileno(), 1)

__author__ = 'Anthony Ma <anthonyma27@gmail.com>, Rasmus Fonseca <fonseca.rasmus@gmail.com>'
__license__ = "Apache License 2.0"
