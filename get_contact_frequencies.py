#!/usr/bin/env python
"""
Determines the frequencies of residue pair contacts in molecular
dynamics simulations. Given one or more MDContact outputs, this
script determines the frequency of each unique interaction of the
form (itype, residue 1, residue2), weighted by number of frames,
across all inputs.

The inputs are one or more MDContact output file paths as well as an
output path. The user may also specify a subset of interaction types
to compute frequencies for. The user may additionally provide a label
file to convert residue labellings (typically for the use of aligning
sequences for performing frequency comparisons with other
trajectories).

The output is a single tsv file with each row indicating residue
id 1, residue id 2, and contact frequency.
"""

from __future__ import division
from collections import defaultdict
import sys
import argparse


def atomid_to_resid(atom):
    return atom[0:atom.rfind(":")]
    # return ':'.join(atom.split(':')[1:3])


def gen_counts(input_lines, interaction_types, residuelabels=None):
    """
    Parse each line in `input_lines` as a line from MDContacts and return interaction-counts for each residue pair. If
    `residuelabels` is defined it is used to modify residue identifiers and to filter out residues not indicated.

    For example:
        inputs = [
            "# total_frames: 3",
            "\t".join(["0", "hbbb", "A:ALA:1:N", "A:ARG:4:O"]),
            "\t".join(["0", "vdw", "A:ALA:1:CB", "A:ARG:4:CA"]),
            "\t".join(["1", "vdw", "A:ALA:1:N", "A:CYS:5:CA"]),
            "\t".join(["2", "hbbb", "A:THR:2:N", "A:CYS:5:O"]),
            "\t".join(["2", "hbss", "A:ALA:1:N", "A:CYS:5:O"])
        ]
        labels = {"A:ALA:1": "A1", "A:ARG:4": "R4", "A:CYS:5": "C5"}

        # Only consider hbbb and vdw, filter away THR, and map to single-letter labels
        gen_counts(inputs, ["hbbb", "vdw"], labels)
        #  Returns: { ("A1", "R4"): 1, ("A1", "C5"): 1 }

    Parameters
    ----------
    input_lines: Iterable[str]
        Interactions formatted as MDContacts output, e.g. ["0\thbbb\tA:ALA:1:N\tA:ARG:4:H", ...]
    interaction_types: list of str
        Which interaction types to consider
    residuelabels: dict of (str: str)
        Remaps and filters residuelabels, e.g. {"A:ARG:4": "R4"}

    Returns
    -------
    (int, dict of (str, str): int)
        Total frame-count and mapping of residue-residue interactions to frame-count
    """
    # Maps residue pairs to set of frames in which they're present
    rescontact_frames = defaultdict(set)
    total_frames = 0

    for line in input_lines:
        line = line.strip()
        if "total_frames" in line:
            tokens = line.split(" ")
            total_frames = int(tokens[1][tokens[1].find(":")+1:])

        if len(line) == 0 or line[0] == "#":
            continue

        tokens = line.split("\t")

        # Check that the interaction type is specified
        itype = tokens[1]
        if itype not in interaction_types:
            continue

        frame = int(tokens[0])
        res1 = atomid_to_resid(tokens[2])
        res2 = atomid_to_resid(tokens[3])

        # Change residue id according to `residuelabels` or skip if any of the residues are not present
        if residuelabels is not None:
            if res1 not in residuelabels or res2 not in residuelabels:
                continue
            res1 = residuelabels[res1]
            res2 = residuelabels[res2]

        # Ensure lexicographical order of residue names
        if res2 < res1:
            res1, res2 = res2, res1

        rescontact_frames[(res1, res2)].add(frame)

    # Insted of returning list of frames for each interaction, only return number of frames
    rescontact_counts = {(res1, res2): len(frames) for (res1, res2), frames in rescontact_frames.items()}
    return total_frames, rescontact_counts


def parse_labelfile(label_file):
    """
    Parses a label-file and returns a dictionary with the residue label mappings. Unless prepended with a comment-
    indicator (#), each line is assumed to have a valid residue identifier (e.g. "A:ALA:1") and a label which the
    residue should be mapped to (e.g. "A1").

    Example:
        parse_labelfile(["A:ALA:1\tA1")
        # Returns {"A:ALA:1": "A1"}

    Parameters
    ----------
    label_file: Iterable[str]
        Lines with tab-separated residue identifier and label

    Returns
    -------
    dict of str: str
        Mapping from residue-id in contact-file to label of any format

    """
    ret = {}
    for line in label_file:
        line = line.strip()
        # Ignore line if empty or comment
        if line[0] == "#" or len(line) == 0:
            continue

        tokens = line.split("\t")
        ret[tokens[0]] = tokens[1]

    return ret


def gen_frequencies(count_list):
    """
    Take a list of residue contact counts (see output of `gen_counts`) and compute total counts and frequencies.

    Example:
        clist = [
            (4, {("A1", "R4"): 4, ("A1", "C5"): 3}),  # First simulation has 4 frames and two contacts
            (3, {("A1", "R4"): 2})                    # Second simulation has 3 frames and one contact
        ]
        gen_frequencies(clist)
        # Returns: (7, {("A1", "R4"): (6, 0.857), ("A1", "C5"): (3, 0.429)})

    Parameters
    ----------
    count_list: list of (int, dict of (str, str): int)
        List with individual frame counts and dictionaries mapping residue pairs to frame-counts

    Return
    ------
    (int, dict of (str, str): (int, float))
        Total framecount and mapping of residue ID pairs to the number of frames in which they contact and the frequency
    """
    rescontact_count = defaultdict(int)
    total_frames = 0
    for frames, rescount_dict in count_list:
        total_frames += frames

        for (res1, res2), count in rescount_dict.items():
            rescontact_count[(res1, res2)] += count

    respair_freqs = {respair: (count, float(count) / total_frames) for respair, count in rescontact_count.items()}
    return total_frames, respair_freqs


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
                             '* hls, hlb (ligand-sidechain and ligand-backbone hydrogen bonds), \n'
                             '* lwb, lwb2 (ligand water-bridges and extended water-bridges)')

    # results, unknown = parser.parse_known_args()
    args = parser.parse_args()

    # Update itypes if "all" is specified
    if "all" in args.itypes:
        args.itypes = ["sb", "pc", "ps", "ts", "vdw", "hb", "lhb", "hbbb", "hbsb",
                       "hbss", "wb", "wb2", "hls", "hlb", "lwb", "lwb2"]

    output_file = args.output_file
    input_files = args.input_files
    itypes = args.itypes
    labels = parse_labelfile(args.label_file) if args.label_file else None

    counts = [gen_counts(input_file, itypes, labels) for input_file in input_files]
    total_frames, frequencies = gen_frequencies(counts)

    output_file.write('#\ttotal_frames:%d\tinteraction_types:%s\n' % (total_frames, ','.join(itypes)))
    output_file.write('#\tColumns:\tresidue_1,\tresidue_2\tframe_count\tcontact_frequency\n')
    for (res1, res2), (count, frequency) in frequencies.items():
        output_file.write('\t'.join([res1, res2, str(count), "%.3f" % frequency]) + "\n")


if __name__ == '__main__':
    main()
