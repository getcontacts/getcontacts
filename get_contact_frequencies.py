#!/usr/bin/env python
"""
Determines the frequencies of residue pair contacts in molecular dynamics simulations.
Given one or more get_dynamic_contacts outputs, this script determines the frequency of
each unique interaction of the form (itype, residue 1, residue2), weighted by number of
frames, across all inputs.

The inputs are one or more get_dynamic_countacts output file paths as well as an
output path. The user may also specify a subset of interaction types to compute
frequencies for. The user may additionally provide a label file to convert residue
labellings (typically for the use of aligning sequences for performing frequency
comparisons with other trajectories).

The output is a single tsv file with each row indicating residue id 1, residue id 2,
and contact frequency.
"""

from __future__ import division
import sys
import argparse
from contact_calc.transformations import *


# def parse_labelfile(label_file):
#     """
#     Parses a label-file and returns a dictionary with the residue label mappings. Unless prepended with a comment-
#     indicator (#), each line is assumed to have a valid residue identifier (e.g. "A:ALA:1") and a label which the
#     residue should be mapped to (e.g. "A1").
#
#     Example:
#         parse_labelfile(["A:ALA:1\tA1")
#         # Returns {"A:ALA:1": "A1"}
#
#     Parameters
#     ----------
#     label_file: Iterable[str]
#         Lines with tab-separated residue identifier and label
#
#     Returns
#     -------
#     dict of str: str
#         Mapping from residue-id in contact-file to label of any format
#
#     """
#     ret = {}
#     for line in label_file:
#         line = line.strip()
#         # Ignore line if empty or comment
#         if line[0] == "#" or len(line) == 0:
#             continue
#
#         tokens = line.split("\t")
#         ret[tokens[0]] = tokens[1]
#
#     return ret


def main(argv=None):
    # Parse command line arguments
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            # Prints full program help when error occurs
            self.print_help(sys.stderr)
            sys.stderr.write('\nError: %s\n' % message)
            sys.exit(2)

    parser = MyParser(description=__doc__,
                      formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--input_files',
                        type=argparse.FileType('r'),
                        required=True,
                        nargs='+',
                        metavar='FILE.tsv',
                        help="Path to one or more contact-file outputs")
    parser.add_argument('--label_file',
                        type=argparse.FileType('r'),
                        required=False,
                        metavar='FILE.tsv',
                        help="A label file for standardizing residue names between different proteins")
    parser.add_argument('--output_file',
                        type=argparse.FileType('w'),
                        required=True,
                        metavar='FILE.tsv',
                        help="Path to output file")
    parser.add_argument('--itypes',
                        required=False,
                        default="all",
                        type=str,
                        nargs="+",
                        metavar="ITYPE",
                        help='Include only these interaction types in frequency computation. Valid choices are: \n'
                             '* all (default), \n'
                             '* sb (salt-bridges), \n'
                             '* pc (pi-cation), \n'
                             '* ps (pi-stacking), \n'
                             '* ts (t-stacking), \n'
                             '* vdw (van der Waals), \n'
                             '* hbbb, hbsb, hbss, (hydrogen bonds with specific backbone/side-chain profile)\n'
                             '* wb, wb2 (water-bridges and extended water-bridges) \n'
                             '* hbls, hblb (ligand-sidechain and ligand-backbone hydrogen bonds), \n'
                             '* lwb, lwb2 (ligand water-bridges and extended water-bridges)')

    # results, unknown = parser.parse_known_args()
    args = parser.parse_args(argv)

    # Update itypes if "all" is specified
    if "all" in args.itypes:
        args.itypes = ["sb", "pc", "ps", "ts", "vdw", "hb", "lhb", "hbbb", "hbsb",
                       "hbss", "wb", "wb2", "hbls", "hblb", "lwb", "lwb2", "br"]

    output_file = args.output_file
    input_files = args.input_files
    itypes = args.itypes
    # labels = parse_labelfile(args.label_file) if args.label_file else None
    labels, colors = (None, None)
    if args.label_file is not None:
        labels, colors = parse_residuelabels(args.label_file)
        args.label_file.close()

    counts = []
    for input_file in input_files:
        contacts, num_frames = parse_contacts(input_file, itypes)
        input_file.close()
        residue_contacts = res_contacts(contacts)
        residue_contacts = relabel(residue_contacts, labels)
        counts.append((num_frames, gen_counts(residue_contacts)))

    total_frames, frequencies = gen_frequencies(counts)

    output_file.write('#\ttotal_frames:%d\tinteraction_types:%s\n' % (total_frames, ','.join(itypes)))
    output_file.write('#\tColumns:\tresidue_1,\tresidue_2\tframe_count\tcontact_frequency\n')
    for (res1, res2), (count, frequency) in frequencies.items():
        output_file.write('\t'.join([res1, res2, "%.3f" % frequency]) + "\n")

    output_file.close()


if __name__ == '__main__':
    main()
