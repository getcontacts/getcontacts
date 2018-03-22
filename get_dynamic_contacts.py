#!/usr/bin/env python3

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

import argparse
from contact_calc.compute_contacts import *

HELP_STR = """
 ============================================================================
| get_dynamic_contacts: A Python Library for computing non-covalent contacts |
|                      throughout Molecular Dynamics Trajectories.           |
|                      Version 1.1.0                                         |
|                                                                            |
| Contact: Anthony Kai Kwang Ma, anthonyma27@gmail.com                       |
 ============================================================================

Command was:
python3 get_dynamic_contacts.py --help

usage: python3 get_dynamic_contacts.py [--help] [--topology TOPOLOGY] 
                      [--trajectory TRAJECTORY]
                      [--output OUTPUT_PATH] 
                      [--cores NUM_CORES]
                      [--solv SOLVENT]
                      [--sele SELECTION]
                      [--ligand LIGAND]
                      [--itype INTERACTION_TYPES]
                      [--sb_cutoff_dist SALT_BRIDGE_CUTOFF_DISTANCE]
                      [--pc_cutoff_dist PI_CATION_CUTOFF_DISTANCE]
                      [--pc_cutoff_ang PI_CATION_CUTOFF_ANGLE]
                      [--ps_cutoff_dist PI_STACK_CUTOFF_DISTANCE]
                      [--ps_cutoff_ang PI_STACK_CUTOFF_ANGLE]
                      [--ps_psi_ang PI_STACK_PSI_ANGLE]
                      [--ts_cutoff_dist T_STACK_CUTOFF_DISTANCE]
                      [--ts_cutoff_ang T_STACK_CUTOFF_ANGLE]
                      [--ts_psi_ang T_STACK_PSI_ANGLE]
                      [--hbond_cutoff_dist HBOND_CUTOFF_DISTANCE]
                      [--hbond_cutoff_ang HBOND_CUTOFF_ANGLE]
                      [--vdw_epsilon VDW_EPSILON]


required arguments:
    --topology TOPOLOGY             path to topology file 
    --trajectory TRAJECTORY         path to trajectory file
    --output OUTPUT_PATH            path to output file
    --itype INTERACTION_TYPES       list of interaction type flags

optional arguments:
    --help                  show this help message and exit
    --cores NUM_CORES       number of cpu cores to parallelize upon [default = 6]
    --solv SOLVENT          resname of solvent molecule [default = "TIP3"]
    --sele SELECTION        atom selection query in VMD [default = None]
    --ligand LIGAND         resname of ligand molecule [default = None]
    --stride STRIDE         skip frames with specified frequency [default = 1]
    --skip SKIP             skip specified number of frames at beggining of trajectory [default = 0]

geometric criteria options:
    --sb_cutoff_dist SALT_BRIDGE_CUTOFF_DISTANCE
                    cutoff for distance between anion and cation 
                    atoms [default = 4.0 angstroms]
    --pc_cutoff_dist PI_CATION_CUTOFF_DISTANCE
                    cutoff for distance between cation and centroid
                    of aromatic ring [default = 6.0 angstroms]
    --pc_cutoff_ang PI_CATION_CUTOFF_ANGLE
                    cutoff for angle between normal vector projecting
                    from aromatic plane and vector from aromatic center
                    to cation atom [default = 60 degrees]
    --ps_cutoff_dist PI_STACK_CUTOFF_DISTANCE
                    cutoff for distance between centroids of two aromatic
                    rings [default = 7.0 angstroms]
    --ps_cutoff_ang PI_STACK_CUTOFF_ANGLE
                    cutoff for angle between the normal vectors projecting
                    from each aromatic plane [default = 30 degrees]
    --ps_psi_ang PI_STACK_PSI_ANGLE
                    cutoff for angle between normal vector projecting from
                    aromatic plane 1 and vector between the two aromatic
                    centroids [default = 45 degrees]
    --ts_cutoff_dist T_STACK_CUTOFF_DISTANCE
                    cutoff for distance between centroids of two aromatic
                    rings [default = 5.0 angstroms]
    --ts_cutoff_ang T_STACK_CUTOFF_ANGLE
                    cutoff for angle between the normal vectors projecting
                    from each aromatic plane minus 90 degrees [default = 30 degrees]
    --ts_psi_ang T_STACK_PSI_ANGLE
                    cutoff for angle between normal vector projecting from
                    aromatic plane 1 and vector between the two aromatic
                    centroids [default = 45 degrees]
    --hbond_cutoff_dist HBOND_CUTOFF_DISTANCE
                    cutoff for distance between donor and acceptor atoms 
                    [default = 3.5 angstroms]
    --hbond_cutoff_ang HBOND_CUTOFF_ANGLE
                    cutoff for angle between donor hydrogen acceptor 
                    [default = 70 degrees]
    --vdw_epsilon VDW_EPSILON
                    amount of padding for calculating vanderwaals contacts 
                    [default = 0.5 angstroms]
    --vdw_res_diff VDW_RES_DIFF
                    minimum residue distance for which to consider computing 
                    vdw interactions [default = 2]


interaction type flags:
    sb             salt bridges 
    pc             pi-cation 
    ps             pi-stacking
    ts             t-stacking
    vdw            vanderwaals
    hb             hydrogen bonds
    lhb            ligand hydrogen bonds

output interaction subtypes:
    hbbb           backbone-backbone hydrogen bonds
    hbsb           backbone-sidechain hydrogen bonds
    hbss           sidechain-sidechain hydrogen bonds
    wb             water bridges
    wb2            extended water bridges 
    lwb            ligand water bridges
    lwb2           extended ligand water bridges 
    hls            ligand-sidechain residue hydrogen bonds 
    hlb            ligand-backbone residue hydrogen bonds 


Examples
========

Compute salt bridges and hydrogen bonds for residues 100 to 160:
    python get_dynamic_contacts.py --topology TOP.pdb \
    --trajectory TRAJ.nc \
    --output output.tsv \
    --cores 12 \
    --solv IP3 \
    --sele "chain A and resid 100 to 160" \
    --ligand EJ4 \
    --itype sb hb lhb

Pi-cation, pi-stacking, and vanderwaals contacts in the entire protein:
python get_dynamic_contacts.py --topology TOP.psf --trajectory TRAJ.dcd --output output.tsv --cores 6 --itype pc ps vdw

Salt bridges and hydrogen bonds in the entire protein with modified distance cutoffs:
python get_dynamic_contacts.py --topology TOP.mae --trajectory TRAJ.dcd --output output.tsv --cores 6 --sb_cutoff_dist 5.0 --hbond_cutoff_dist 4.5 --itype sb hb
"""

DESCRIPTION = "Computes non-covalent contact networks in MD simulations."


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
        "VDW_EPSILON": args.vdw_epsilon,
        "VDW_RES_DIFF": args.vdw_res_diff
    }
    return geom_criterion_values


def main(traj_required=True):
    if "--help" in sys.argv or "-h" in sys.argv:
        print(HELP_STR)
        exit(1)

    # Parse required and optional arguments
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--topology', type=str, required=True, help='path to topology file ')
    parser.add_argument('--trajectory', type=str, required=traj_required, default=None, help='path to trajectory file')
    parser.add_argument('--output', type=str, required=True, help='path to output file')
    parser.add_argument('--cores', type=int, default=6, help='number of cpu cores to parallelize upon')
    parser.add_argument('--solv', type=str, default="TIP3", help='resname of solvent molecule')
    parser.add_argument('--sele', type=str, default=None, help='atom selection query in VMD')
    parser.add_argument('--stride', type=int, default=1, help='skip frames with specified frequency')
    parser.add_argument('--skip', type=int, default=0, help='skip specified number of frames at beggining of trajectory')
    parser.add_argument('--ligand', type=str, nargs="+", default=[], help='resname of ligand molecule')

    # Parse geometric criterion arguments
    parser.add_argument('--sb_cutoff_dist', type=float, default=4.0, help='cutoff for distance between anion and cation atoms [default = 4.0 angstroms]')
    parser.add_argument('--pc_cutoff_dist', type=float, default=6.0, help='cutoff for distance between cation and centroid of aromatic ring [default = 6.0 angstroms]')
    parser.add_argument('--pc_cutoff_ang', type=float, default=60, help='cutoff for angle between normal vector projecting from aromatic plane and vector from aromatic center to cation atom [default = 60 degrees]')
    parser.add_argument('--ps_cutoff_dist', type=float, default=7.0, help='cutoff for distance between centroids of two aromatic rings [default = 7.0 angstroms]')
    parser.add_argument('--ps_cutoff_ang', type=float, default=30, help='cutoff for angle between the normal vectors projecting from each aromatic plane [default = 30 degrees]')
    parser.add_argument('--ps_psi_ang', type=float, default=45, help='cutoff for angle between normal vector projecting from aromatic plane 1 and vector between the two aromatic centroids [default = 45 degrees]')
    parser.add_argument('--ts_cutoff_dist', type=float, default=5.0, help='cutoff for distance between centroids of two aromatic rings [default = 5.0 angstroms]')
    parser.add_argument('--ts_cutoff_ang', type=float, default=30, help='cutoff for angle between the normal vectors projecting from each aromatic plane minus 90 degrees [default = 30 degrees]')
    parser.add_argument('--ts_psi_ang', type=float, default=45, help='cutoff for angle between normal vector projecting from aromatic plane 1 and vector between the two aromatic centroids [default = 45 degrees]')
    parser.add_argument('--hbond_cutoff_dist', type=float, default=3.5, help='cutoff for distance between donor and acceptor atoms [default = 3.5 angstroms]')
    parser.add_argument('--hbond_cutoff_ang', type=float, default=70, help='cutoff for angle between donor hydrogen acceptor [default = 70 degrees]')
    parser.add_argument('--vdw_epsilon', type=float, default=0.5, help='amount of padding for calculating vanderwaals contacts [default = 0.5 angstroms]')
    parser.add_argument('--vdw_res_diff', type=int, default=2, help='minimum residue distance for which to consider computing vdw interactions')


    parser.add_argument('--itypes',
                        required=True,
                        type=str,
                        nargs="+",
                        metavar="ITYPE",
                        help='Compute only these interaction types. Valid choices are: \n'
                             '* all (default), \n'
                             '* sb (salt-bridges), \n'
                             '* pc (pi-cation), \n'
                             '* ps (pi-stacking), \n'
                             '* ts (t-stacking), \n'
                             '* vdw (van der Waals), \n'
                             '* hb (hydrogen bonds)\n'
                             '* lhb (ligand hydrogen bonds)')

    args, unknown = parser.parse_known_args()

    top = args.topology
    traj = args.trajectory
    output = args.output
    cores = args.cores
    ligand = args.ligand
    solv = args.solv
    sele = args.sele
    stride = args.stride
    skip = args.skip
    geom_criterion_values = process_geometric_criterion_args(args)

    # Check interaction types
    all_itypes = ["sb", "pc", "ps", "ts", "vdw", "hb", "lhb"]
    if "all" in args.itypes:
        itypes = all_itypes
    else:
        for itype in args.itypes:
            if itype not in all_itypes:
                print("Error: " + itype + " is not a valid interaction type")
                exit(1)

        itypes = args.itypes

    # Begin computation
    tic = datetime.datetime.now()
    compute_contacts(top, traj, output, itypes, geom_criterion_values, cores, stride, skip, solv, sele, ligand)
    toc = datetime.datetime.now()
    print("Computation time: " + str((toc-tic).total_seconds()) + " seconds")

    print("topology=%s" % top)
    print("trajectory=%s" % traj)
    print("output=%s" % output)
    print("cores=%s" % cores)
    print("ligand=%s" % ",".join(ligand))
    print("solv=%s" % solv)
    print("sele=%s" % sele)
    print("stride=%s" % stride)


if __name__ == "__main__":
    main()

    # Suppress stdout from vmd as program terminates
    devnull = open('/dev/null', "w")
    os.dup2(devnull.fileno(), 1)
