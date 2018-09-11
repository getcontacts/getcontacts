#!/usr/bin/env python
"""
Take a set of residue-frequency files generated by `get_contact_frequencies.py`,
group them into a single table file by matching residue pair ids, and plot them
as a clustered heat-map, a tsv table, or a multi-flare.
"""

from __future__ import division
import contact_calc.argparsers as ap
import argparse
import numpy as np
from contact_calc.flare import compose_frequencytable, write_json
from contact_calc.transformations import parse_frequencyfiles


def write_frequencytable(freq_table, col_labels, fname, cluster_columns=True, cluster_rows=True):
    row_labels = [(r1, r2) for r1, r2 in freq_table]
    freq_matrix = np.array([freq_table[(r1, r2)] for (r1, r2) in freq_table])
    m, n = freq_matrix.shape
    if cluster_rows:
        from scipy.cluster.hierarchy import linkage, leaves_list
        l = linkage(freq_matrix, method='single')
        row_ordering = leaves_list(l)
    else:
        row_ordering = range(m)

    if cluster_columns:
        from scipy.cluster.hierarchy import linkage, leaves_list
        l = linkage(freq_matrix.T, method='single')
        col_ordering = leaves_list(l)
    else:
        col_ordering = range(n)

    freq_matrix = freq_matrix[row_ordering]
    freq_matrix = freq_matrix[:, col_ordering]
    row_labels = [row_labels[i] for i in row_ordering]
    col_labels = [col_labels[i] for i in col_ordering]

    with open(fname, "w") as out_file:
        out_file.write("\t".join(["", ""] + col_labels) + "\n")
        for i in range(m):
            res1, res2 = row_labels[i]
            freq_strings = [str(freq) for freq in freq_matrix[i]]
            out_file.write("\t".join([res1, res2] + freq_strings) + "\n")


def write_pymol_distances(multiflare, fname):
    """
    TODO: Document
    """
    from collections import defaultdict
    # num_frames = max(multiflare["edges"], lambda e: e["frames"][-1]) + 1
    num_frames = 0
    for e in multiflare["edges"]:
        num_frames = max(num_frames, e["frames"][-1] + 1)
    iprofiles = map(lambda e: (e["name1"], e["name2"], str(e["frames"])), multiflare["edges"])
    iprofile_hist = defaultdict(list)
    for (n1, n2, iprofile) in iprofiles:
        iprofile_hist[iprofile].append((n1, n2))

    with open(fname, "w") as f:
        prioritized_iprofiles = sorted(iprofile_hist.keys(), key=lambda k: len(iprofile_hist[k]), reverse=True)
        for iprofile in prioritized_iprofiles:
            iprofile_dec = ["-" for _ in range(num_frames)]
            for cond in map(int, iprofile.strip('[]').split(",")):
                print(cond)
                iprofile_dec[cond] = '+'
            iprofile_dec = "row_" + "".join(iprofile_dec)

            for (n1, n2) in iprofile_hist[iprofile]:
                if n1 == n2:
                    continue
                chain1 = n1.split(":")[0]
                resi1 = n1.split(":")[2]
                chain2 = n2.split(":")[0]
                resi2 = n2.split(":")[2]
                f.write("distance %s, ///%s/%s/CA, ///%s/%s/CA\n" % (iprofile_dec, chain1, resi1, chain2, resi2))


def plot_frequencies(freq_table, col_labels, out_file, cluster_columns):
    import pandas as pd
    import matplotlib
    import os
    # if "DISPLAY" not in os.environ:
        # matplotlib.use('agg')
    matplotlib.use('Agg')

    import seaborn as sns; 
    sns.set(color_codes=True)
    sns.set(font_scale=1.5)

    freq_matrix = np.array([freq_table[(r1, r2)] for (r1, r2) in freq_table])
    row_labels = [r1 + " - " + r2 for (r1, r2) in freq_table]
    pdframe = pd.DataFrame(freq_matrix, index=row_labels, columns=col_labels)

    # Scale down figsize if too large
    figsize = [pdframe.shape[1], pdframe.shape[0]]
    if figsize[1] > 320:
        figsize[0] *= 320 / figsize[1]
        figsize[1] *= 320 / figsize[1]

    # Create clustermap
    fingerprints = sns.clustermap(pdframe,
                                  figsize=figsize,
                                  annot=False,
                                  col_cluster=cluster_columns,
                                  linewidths=0.5,
                                  linecolor='black',
                                  cmap='Greens')

    # Remove color bar
    # fingerprints.cax.set_visible(False)

    import matplotlib.pyplot as plt
    plt.setp(fingerprints.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.setp(fingerprints.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)

    fingerprints.savefig(out_file)


def main(argv=None):
    # Parse command line arguments
    parser = ap.PrintUsageParser(description=__doc__)
    parser.add_argument('--input_frequencies',
                        type=argparse.FileType('r'),
                        required=True,
                        nargs='+',
                        help="Paths to one or more residue frequency files")
    parser.add_argument('--frequency_cutoff',
                        type=float,
                        required=False,
                        default=0.6,
                        help="Only interactions occurring at least this frequently will be plotted (default: 0.6)")
    parser.add_argument('--column_headers',
                        type=str,
                        required=False,
                        nargs='+',
                        help="Header column labels. If nothing is specified, the input_frequencies filenames are used")
    parser.add_argument('--cluster_columns',
                        type=bool,
                        required=False,
                        default=False,
                        help="Perform hierarchical clustering on the columns (default: False)")
    parser.add_argument('--table_output',
                        type=str,
                        required=False,
                        default=None,
                        help="If specified, the tab-separated frequency table will be written to this file")
    parser.add_argument('--plot_output',
                        type=str,
                        required=False,
                        default=None,
                        help="If specified, the heatmap will be written to this file (supports svg and png formats)")
    parser.add_argument('--flare_output',
                        type=str,
                        required=False,
                        default=None,
                        help="If specified, a compare-flare will be written to this json-file")
    parser.add_argument('--pymol_output',
                        type=str,
                        required=False,
                        default=None,
                        help="If specified, a distance-selection will be written to this pml-file")

    args = parser.parse_args(argv)

    freq_table = parse_frequencyfiles(args.input_frequencies, args.frequency_cutoff)

    # Determine column headers and exit on error
    column_headers = [f.name for f in args.input_frequencies] if args.column_headers is None else args.column_headers
    if len(column_headers) != len(args.input_frequencies):
        parser.error("--column_header arguments must match length of --input_frequencies")

    # Check output format and call corresponding function(s)
    if all(a is None for a in [args.table_output, args.flare_output, args.plot_output, args.pymol_output]):
        parser.error("--table_output, --flare_output, or --plot_output must be specified")

    if args.table_output is not None:
        write_frequencytable(freq_table, column_headers, args.table_output)
        print("Wrote frequency table to "+args.table_output)

    if args.flare_output is not None:
        compare_flare = compose_frequencytable(freq_table, column_headers, args.frequency_cutoff)
        write_json(compare_flare, args.flare_output)
        print("Wrote multi flare to "+args.flare_output)

    if args.plot_output is not None:
        plot_frequencies(freq_table, column_headers, args.plot_output, args.cluster_columns)
        print("Wrote fingerprint heatmap to "+args.plot_output)

    if args.pymol_output is not None:
        compare_flare = compose_frequencytable(freq_table, column_headers, args.frequency_cutoff)
        write_pymol_distances(compare_flare, args.pymol_output)
        print("Wrote pymol file to "+args.pymol_output)

    for f in args.input_frequencies:
        f.close()


if __name__ == '__main__':
    main()


__author__ = 'Rasmus Fonseca <fonseca.rasmus@gmail.com>'
__license__ = "Apache License 2.0"
