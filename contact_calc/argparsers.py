# -*- coding: utf-8 -*-

import sys
import argparse


class PrintUsageParser(argparse.ArgumentParser):
    """
    A custom argparser that prints usage and error (not only error) when encountering a problem.
    Also sets the formatter_class to `RawTextHelpFormatter` so `__doc__` can be used as description.
    """
    def __init__(self, *args, **kwargs):
        kwargs['formatter_class'] = argparse.RawTextHelpFormatter
        super(PrintUsageParser, self).__init__(*args, **kwargs)

    def error(self, message):
        """
        Prints full program usage and exits with errorcode 2 when error occurs
        """
        self.print_help(sys.stderr)
        sys.stderr.write('\nError: %s\n' % message)
        sys.exit(2)


def build_getcontact_parser(program_description, trajectory):
    parser = PrintUsageParser(description=program_description, add_help=False)

    required_group = parser.add_argument_group("Required arguments")

    if trajectory:
        required_group.add_argument("--topology", type=str, metavar="PATH", required=True,
                                    help="path to topology file ")
        required_group.add_argument("--trajectory", type=str, metavar="PATH", required=True,
                                    help="path to trajectory file")
    else:
        required_group.add_argument("--structure", type=str, metavar="PATH", required=True,
                                    help="path to structure file ")

    required_group.add_argument("--output", type=str, metavar="PATH", required=True,
                                help="path to output file")
    required_group.add_argument("--itypes", type=str, nargs="+", required=True, metavar="ITYPE",
                                help="Compute only these interaction types. Valid choices are: \n"
                                     "· all, \n"
                                     "· sb (salt-bridges), \n"
                                     "· pc (pi-cation), \n"
                                     "· ps (pi-stacking), \n"
                                     "· ts (t-stacking), \n"
                                     "· vdw (van der Waals), \n"
                                     "· hb (hydrogen bonds)\n"
                                     "· lhb (ligand hydrogen bonds)")

    optional_group = parser.add_argument_group("Optional arguments")
    optional_group.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS,
                                help="show this help message and exit.")
    optional_group.add_argument("--solv", type=str, metavar="SOLVNAME", default="TIP3",
                                help="resname of solvent molecule")
    optional_group.add_argument("--sele", type=str, metavar="VMDSEL", default=None,
                                help="atom selection 1 query in VMD")
    optional_group.add_argument("--sele2", type=str, metavar="VMDSEL", default=None,
                                help="atom selection 2 query in VMD")
    optional_group.add_argument("--ligand", type=str, metavar="LIGNAME", nargs="*", default=[],
                                help="resname of ligand molecules")

    if trajectory:
        optional_group.add_argument("--cores", type=int, metavar="INT", default=6,
                                    help="number of cpu cores to parallelize on [default = 6]")
        optional_group.add_argument("--beg", type=int, metavar="INT", default=0,
                                    help="first frame to read [default = 0]")
        optional_group.add_argument("--end", type=int, metavar="INT", default=None,
                                    help="last frame to read (unless end-of-trajectory is met first) [default = infty]")
        optional_group.add_argument("--stride", type=int, metavar="INT", default=1,
                                    help="steps between computed frames [default = 1]")

    # Parse geometric criterion arguments
    geometric_group = parser.add_argument_group("Geometric criteria")
    geometric_group.add_argument("--sb_cutoff_dist", type=float, metavar="FLOAT", default=4.0,
                                 help="cutoff for distance between anion and cation atoms [default = 4.0 angstroms]")
    geometric_group.add_argument("--pc_cutoff_dist", type=float, metavar="FLOAT", default=6.0,
                                 help="cutoff for distance between cation and centroid of aromatic ring "
                                      "[default = 6.0 angstroms]")
    geometric_group.add_argument("--pc_cutoff_ang", type=float, metavar="FLOAT", default=60,
                                 help="cutoff for angle between normal vector projecting from aromatic plane and "
                                      "vector from aromatic center to cation atom [default = 60 degrees]")
    geometric_group.add_argument("--ps_cutoff_dist", type=float, metavar="FLOAT", default=7.0,
                                 help="cutoff for distance between centroids of two aromatic rings "
                                      "[default = 7.0 angstroms]")
    geometric_group.add_argument("--ps_cutoff_ang", type=float, metavar="FLOAT", default=30,
                                 help="cutoff for angle between the normal vectors projecting from each aromatic plane "
                                      "[default = 30 degrees]")
    geometric_group.add_argument("--ps_psi_ang", type=float, metavar="FLOAT", default=45,
                                 help="cutoff for angle between normal vector projecting from aromatic plane 1 and "
                                      "vector between the two aromatic centroids [default = 45 degrees]")
    geometric_group.add_argument("--ts_cutoff_dist", type=float, metavar="FLOAT", default=5.0,
                                 help="cutoff for distance between centroids of two aromatic rings "
                                      "[default = 5.0 angstroms]")
    geometric_group.add_argument("--ts_cutoff_ang", type=float, metavar="FLOAT", default=30,
                                 help="cutoff for angle between the normal vectors projecting from each aromatic "
                                      "plane minus 90 degrees [default = 30 degrees]")
    geometric_group.add_argument("--ts_psi_ang", type=float, metavar="FLOAT", default=45,
                                 help="cutoff for angle between normal vector projecting from aromatic plane 1 and "
                                      "vector between the two aromatic centroids [default = 45 degrees]")
    geometric_group.add_argument("--hbond_cutoff_dist", type=float, metavar="FLOAT", default=3.5,
                                 help="cutoff for distance between donor and acceptor atoms [default = 3.5 angstroms]")
    if(trajectory):
        geometric_group.add_argument("--hbond_cutoff_ang", type=float, metavar="FLOAT", default=70,
                                 help="cutoff for angle between donor hydrogen acceptor [default = 70 degrees]")
    else:
        geometric_group.add_argument("--hbond_cutoff_ang", type=float, metavar="FLOAT", default=180,
                                 help="cutoff for angle between donor hydrogen acceptor [default = 180 degrees]")
    geometric_group.add_argument("--vdw_epsilon", type=float, metavar="FLOAT", default=0.5,
                                 help="amount of padding for calculating vanderwaals contacts "
                                      "[default = 0.5 angstroms]")
    geometric_group.add_argument("--vdw_res_diff", type=int, metavar="INT", default=2,
                                 help="minimum residue distance for which to consider computing vdw interactions "
                                      "[default = 2]")

    return parser

