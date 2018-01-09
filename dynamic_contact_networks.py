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

import os
import errno
import sys
import datetime
import argparse
from contact_calc.compute_contacts import *

HELP_STR = """
 ==========================================================================
| MDContactNetworks: A Python Library for computing non-covalent contacts  |
|                    throughout Molecular Dynamics Trajectories.           |
|                    Version 1.1.0                                         |
|                                                                          |
| Contact: Anthony Kai Kwang Ma, anthonyma27@gmail.com                     |
 ==========================================================================

Command was:
python3 contact_networks.py --help

usage: python3 contact_networks.py [--help] [--topology TOPOLOGY] 
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
    hls            ligand-sidechain residue hydrogen bonds 
    hlb            ligand-backbone residue hydrogen bonds 
    lwb            ligand water bridges
    lwb2           extended ligand water bridges 


examples:
Salt bridges and hydrogen bonds for residues 100 to 160:
python contact_networks.py --topology TOP.pdb --trajectory TRAJ.nc --output output.tsv --cores 12 --solv IP3 --sele "chain A and resid 100 to 160" --ligand EJ4 --itype sb hb lhb

Pi-cation, pi-stacking, and vanderwaals contacts in the entire protein:
python contact_networks.py --topology TOP.psf --trajectory TRAJ.dcd --output output.tsv --cores 6 --itype pc ps vdw

Salt bridges and hydrogen bonds in the entire protein with modified distance cutoffs:
python contact_networks.py --topology TOP.mae --trajectory TRAJ.dcd --output output.tsv --cores 6 --sb_cutoff_dist 5.0 --hbond_cutoff_dist 4.5 --itype sb hb
"""

DESCRIPTION = "Computes non-covalent contact networks in MD simulations."


def clean_path(path):
    if path[-1] != '/':
        path += '/'
    return path


def open_dir(dirname):
    try:
        os.makedirs(dirname)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def validate_itypes(ITYPES):
    """
    Prints error and halts program if specified interaction type is invalid
    """
    valid_itypes = ["sb", "pc", "ps", "ts", "vdw", "hb", "lhb"]

    for itype in ITYPES:
        if itype not in valid_itypes:
            print("%s not a valid interaction type ..." % (itype))
            exit(1)
    if not ITYPES:
        print("Need to specify at least one interaction type ...")
        exit(1)


def process_geometric_criterion_args(args):
    SALT_BRIDGE_CUTOFF_DISTANCE = args.sb_cutoff_dist
    PI_CATION_CUTOFF_DISTANCE = args.pc_cutoff_dist
    PI_CATION_CUTOFF_ANGLE = args.pc_cutoff_ang
    PI_STACK_CUTOFF_DISTANCE = args.ps_cutoff_dist
    PI_STACK_CUTOFF_ANGLE = args.ps_cutoff_ang
    PI_STACK_PSI_ANGLE = args.ps_psi_ang
    T_STACK_CUTOFF_DISTANCE = args.ts_cutoff_dist
    T_STACK_CUTOFF_ANGLE = args.ts_cutoff_ang
    T_STACK_PSI_ANGLE = args.ts_psi_ang
    HBOND_CUTOFF_DISTANCE = args.hbond_cutoff_dist
    HBOND_CUTOFF_ANGLE = args.hbond_cutoff_ang
    VDW_EPSILON = args.vdw_epsilon

    geom_criterion_values = {}
    geom_criterion_values['SALT_BRIDGE_CUTOFF_DISTANCE'] = SALT_BRIDGE_CUTOFF_DISTANCE
    geom_criterion_values['PI_CATION_CUTOFF_DISTANCE'] = PI_CATION_CUTOFF_DISTANCE
    geom_criterion_values['PI_CATION_CUTOFF_ANGLE'] = PI_CATION_CUTOFF_ANGLE
    geom_criterion_values['PI_STACK_CUTOFF_DISTANCE'] = PI_STACK_CUTOFF_DISTANCE
    geom_criterion_values['PI_STACK_CUTOFF_ANGLE'] = PI_STACK_CUTOFF_ANGLE
    geom_criterion_values['PI_STACK_PSI_ANGLE'] = PI_STACK_PSI_ANGLE
    geom_criterion_values['T_STACK_CUTOFF_DISTANCE'] = T_STACK_CUTOFF_DISTANCE
    geom_criterion_values['T_STACK_CUTOFF_ANGLE'] = T_STACK_CUTOFF_ANGLE
    geom_criterion_values['T_STACK_PSI_ANGLE'] = T_STACK_PSI_ANGLE
    geom_criterion_values['HBOND_CUTOFF_DISTANCE'] = HBOND_CUTOFF_DISTANCE
    geom_criterion_values['HBOND_CUTOFF_ANGLE'] =  HBOND_CUTOFF_ANGLE
    geom_criterion_values['VDW_EPSILON'] = VDW_EPSILON
    return geom_criterion_values


def main():
    if "--help" in sys.argv or "-h" in sys.argv:
        print (HELP_STR)
        exit(1)

    # Parse required and optional arguments
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--topology', type=str, required=True, help='path to topology file ')
    parser.add_argument('--trajectory', type=str, required=False, default=None, help='path to trajectory file')
    parser.add_argument('--output', type=str, required=True, help='path to output file')
    parser.add_argument('--cores', type=int, default=6, help='number of cpu cores to parallelize upon')
    parser.add_argument('--solv', type=str, default="TIP3", help='resname of solvent molecule')
    parser.add_argument('--sele', type=str, default=None, help='atom selection query in VMD')
    parser.add_argument('--stride', type=int, default=1, help='skip frames with specified frequency')
    parser.add_argument('--ligand', type=str, default=None, help='resname of ligand molecule')

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

    # Parse interaction types
    class ITypeAction(argparse.Action):
        def __call__(self, parser, ns, values, option):
            if self.dest == "all":
                ns.itypes = set([])
            if not hasattr(ns, "itype"):
                ns.itypes = set()
            ns.itypes.add(self.dest)
    parser.add_argument('--salt-bridge', '-sb', dest='sb', action=ITypeAction, nargs=0, help="Compute salt bridge interactions")
    parser.add_argument('--pi-cation', '-pc', dest='pc', action=ITypeAction, nargs=0, help="Compute pi-cation interactions")
    parser.add_argument('--pi-stacking', '-ps', dest='ps', action=ITypeAction, nargs=0, help="Compute pi-stacking interactions")
    parser.add_argument('--t-stacking', '-ts', dest='ts', action=ITypeAction, nargs=0, help="Compute t-stacking interactions")
    parser.add_argument('--vanderwaals', '-vdw', dest='vdw', action=ITypeAction, nargs=0, help="Compute van der Waals interactions (warning: there will be many)")
    parser.add_argument('--hbond', '-hb', dest='hb', action=ITypeAction, nargs=0, help="Compute hydrogen bond interactions")
    parser.add_argument('--ligand-hbond', '-lhb', dest='lhb', action=ITypeAction, nargs=0, help="Compute ligand hydrogen bond interactions")
    parser.add_argument('--all-interactions', '-all', dest='all', action='store_true', help="Compute all types of interactions")

    args, unknown = parser.parse_known_args()

    top = args.topology
    traj = args.trajectory
    output = args.output
    cores = args.cores
    ligand = args.ligand
    solv = args.solv
    sele = args.sele
    stride = args.stride
    geom_criterion_values = process_geometric_criterion_args(args)

    if args.all:
        itypes = ["sb", "pc", "ps", "ts", "vdw", "hb", "lhb"]
    elif not hasattr(args, "itypes"):
        print("Error: at least one interaction type or '--all-interactions' is required")
        exit(1)
    else:
        itypes = args.itypes
    validate_itypes(itypes)

    # Begin computation
    tic = datetime.datetime.now()
    compute_contacts(top, traj, output, itypes, geom_criterion_values, cores, stride, solv, sele, ligand)
    toc = datetime.datetime.now()
    print("Computation Time: " + str((toc-tic).total_seconds()))

    print("topology=%s" % top)
    print("trajectory=%s" % traj)
    print("output=%s" % output)
    print("cores=%s" % cores)
    print("ligand=%s" % ligand)
    print("solv=%s" % solv)
    print("sele=%s" % sele)
    print("stride=%s" % stride)


if __name__ == "__main__":
    main()
