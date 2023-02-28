#!/usr/bin/env python
"""
Compute a hierarchical clustering dendrogram from a set of residue contact frequencies.
"""

from __future__ import division
import sys
import argparse
import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt


def parse_frequencyfiles(freq_files, freq_cutoff):
    ret = []
    for freq_file in freq_files:
        freq_file_interactions = set()
        for line in freq_file:
            line = line.strip()
            if len(line) == 0 or line[0] == "#":
                continue

            tokens = line.split("\t")
            res1 = tokens[0]
            res2 = tokens[1]
            freq = float(tokens[2])

            if freq > freq_cutoff:
                if res1 > res2:
                    res1, res2 = res2, res1
                freq_file_interactions.add((res1, res2))

        ret.append(freq_file_interactions)

    return ret


def build_dist_matrix(contact_sets):
    n = len(contact_sets)
    ret = np.zeros((n, n))

    for i in range(n):
        for j in range(i):
            jac = len(contact_sets[i] & contact_sets[j]) / float(len(contact_sets[i] | contact_sets[j]))
            dist = 1-jac
            ret[i][j] = dist
            ret[j][i] = dist

    return ret


def main():
    # Parse command line arguments
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            # Prints full program help when error occurs
            self.print_help(sys.stderr)
            sys.stderr.write('\nError: %s\n' % message)
            sys.exit(2)

    parser = MyParser(description=__doc__,
                      formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--input_frequencies',
                        type=argparse.FileType('r'),
                        required=True,
                        metavar="FILE.tsv",
                        nargs='+',
                        help="Paths to one or more residue frequency files")
    parser.add_argument('--frequency_cutoff',
                        type=float,
                        required=False,
                        metavar="FLOAT",
                        default=0.6,
                        help="Only interactions occurring at least this frequently will be plotted (default: 0.6)")
    parser.add_argument('--column_headers',
                        type=str,
                        required=False,
                        metavar="STR",
                        nargs='+',
                        help="Header column labels. If nothing is specified, the input_frequencies filenames are used")
    parser.add_argument('--dendrogram_output',
                        type=str,
                        required=True,
                        metavar="FILE.svg",
                        default=None,
                        help="The dendrogram will be written to this svg-file")

    # results, unknown = parser.parse_known_args()
    args = parser.parse_args()

    if args.column_headers is None:
        args.column_headers = [f.name for f in args.input_frequencies]

    interaction_sets = parse_frequencyfiles(args.input_frequencies, args.frequency_cutoff)
    distance_matrix = squareform(build_dist_matrix(interaction_sets))
    linkage_matrix = linkage(distance_matrix, "single")

    plt.rcParams['lines.linewidth'] = 0.5
    plt.figure(figsize=(20,4))
    dendrogram(linkage_matrix, labels=args.column_headers)
    plt.tight_layout()
    plt.savefig(args.dendrogram_output)
    print("Done, dendrogram SVG saved to "+args.dendrogram_output)


if __name__ == '__main__':
    main()
