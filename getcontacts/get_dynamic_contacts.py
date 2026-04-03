#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Computes non-covalent contact networks for molecular simulations. Detailed
description at https://github.com/getcontacts/getcontacts and further command-
line help by providing the --help argument.

Example usage (hbonds and salt-bridges between chain A loop and ligand):
    get_dynamic_contacts.py --topology TOP.pdb \\
                            --trajectory TRAJ.nc \\
                            --output output.tsv \\
                            --cores 12 \\
                            --solv IP3 \\
                            --sele "(chain A and resid 100 to 160)" \\
                            --sele2 "resname EJ4" \\
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
    parser = ap.build_getcontact_parser(__doc__, True)
    args = parser.parse_args(argv)

    top = args.topology
    traj = args.trajectory
    output = args.output
    cores = args.cores
    ligand = args.ligand
    solv = args.solv
    lipid = args.lipid
    sele1 = args.sele
    sele2 = args.sele2
    beg = args.beg
    end = args.end if args.end else sys.maxsize
    stride = args.stride
    distout = args.distout
    geom_criteria = process_geometric_criterion_args(args)

    # If sele2 is None set it to sele1
    if sele2 is None:
        sele2 = sele1

    # Check interaction types
    all_itypes = ["hp", "sb", "pc", "ps", "ts", "vdw", "hb"]
    if "all" in args.itypes:
        itypes = all_itypes
    else:
        for itype in args.itypes:
            if itype not in all_itypes:
                parser.error(itype + " is not a valid interaction type")

        itypes = args.itypes

    # Begin computation
    tic = datetime.datetime.now()
    compute_contacts(top, traj, output, itypes, geom_criteria, cores, beg, end, stride, distout,
                     ligand, solv, lipid, sele1, sele2)
    toc = datetime.datetime.now()
    total_time = (toc-tic).total_seconds()
    print("\nTotal computation time:", total_time, "seconds")

    print("topology=%s" % top)
    print("trajectory=%s" % traj)
    print("output=%s" % output)
    print("cores=%s" % cores)
    print("solv=%s" % solv)
    print("sele=%s" % sele1)
    print("sele2=%s" % sele2)
    print("stride=%s" % stride)


if __name__ == "__main__":
    main()

    # Suppress stdout from vmd as program terminates
    devnull = open("/dev/null", "w")
    os.dup2(devnull.fileno(), 1)

__author__ = 'Anthony Ma <anthonyma27@gmail.com>, Rasmus Fonseca <fonseca.rasmus@gmail.com>'
__license__ = "Apache License 2.0"
