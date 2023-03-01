#!/usr/bin/env python3

from contact_calc.transformations import parse_contacts, multi_to_single_contact
import argparse as ap

"""
Converts a list of atomic contacts across multiple frames into a single-frame contact
file where atomic interactions are filtered out if their corresponding residues do not
frequently form a interactions of that particular kind.
"""

if __name__ == "__main__":
    # Parse command line arguments
    parser = ap.ArgumentParser(description=__doc__, formatter_class=ap.RawTextHelpFormatter)
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    parser._action_groups.append(optional)  # added this line

    required.add_argument('--input',
                          required=True,
                          type=ap.FileType('r'),
                          metavar='FILE',
                          help='A contact-file or frequency file')
    optional.add_argument('--output',
                          required=False,
                          metavar='FILE',
                          type=ap.FileType('w'),
                          help='The name of the output contact-file')
    optional.add_argument('--min_frequency',
                          required=False,
                          metavar='FLOAT',
                          default=0.6,
                          type=float,
                          help='Minimum residue frequency')

    args = parser.parse_args()

    contacts, num_frames = parse_contacts(args.input)
    contacts = multi_to_single_contact(contacts, args.min_frequency * num_frames)
    contacts = [str(c[0]) + "\t" + "\t".join(c[1:]) for c in contacts]

    # Write to output
    if args.output:
        args.output.write("\n".join(contacts))
        args.output.close()
        print("Wrote residue contact file to " + args.output.name)
    else:
        print("\n".join(contacts))


__license__ = "Apache License 2.0"
__maintainer__ = "Rasmus Fonseca"
__email__ = "fonseca.rasmus@gmail.com"
