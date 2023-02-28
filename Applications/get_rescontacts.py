#!/usr/bin/env python3

"""
Takes a list of atomic contacts as input and generates a list of residue contacts
as output, represented by the CA atom (in case of amino acids) or the first atom
(for any other type of residue).

"""

from collections import defaultdict


if __name__ == "__main__":
    import sys
    from os import path
    sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
    from contact_calc.transformations import parse_contacts, parse_frequencyfiles, ParseError

    # Parse command line arguments
    import argparse as ap
    parser = ap.ArgumentParser(description=__doc__, formatter_class=ap.RawTextHelpFormatter)
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    parser._action_groups.append(optional)  # added this line

    required.add_argument('--input',
                          required=True,
                          type=ap.FileType('r'),
                          metavar='FILE',
                          help='A contact-file or frequency file')
    required.add_argument('--output',
                          required=False,
                          metavar='FILE',
                          type=ap.FileType('w'),
                          help='The name of the output contact-file')

    args = parser.parse_args()
    contacts, total_frames = parse_contacts(args.input)

    # Build a dict mapping residue ids to a set of all its atom names
    representative_atoms = defaultdict(list)
    for c in contacts:
        resi1 = c[2][0:c[2].rindex(":")]
        name1 = c[2][c[2].rindex(":")+1:]
        representative_atoms[resi1].append(name1)
        resi2 = c[3][0:c[3].rindex(":")]
        name2 = c[3][c[3].rindex(":")+1:]
        representative_atoms[resi2].append(name2)

    # Reduce the sets of atoms to a single representative atoms
    representative_atoms = {resi: "CA" if "CA" in names else names[0] for (resi, names) in representative_atoms.items()}

    # Build set of residue contacts (convert atom to representative and remove duplicates
    rescontacts = set()
    for c in contacts:
        frame = str(c[0])
        itype = c[1]
        resi1 = ":".join(c[2].split(":")[0:3])
        resi2 = ":".join(c[3].split(":")[0:3])
        atom1 = resi1 + ":" + representative_atoms[resi1]
        atom2 = resi2 + ":" + representative_atoms[resi2]
        rescontacts.add((frame, itype, atom1, atom2))

    rescontacts = sorted(rescontacts, key=lambda contact: (int(contact[0]), contact[1]))
    rescontacts = ["\t".join(contact) for contact in rescontacts]

    # Write to output
    if args.output:
        args.output.write("\n".join(rescontacts))
        args.output.close()
        print("Wrote residue contact file to " + args.output.name)
    else:
        print("\n".join(rescontacts))


__license__ = "Apache License 2.0"
__maintainer__ = "Rasmus Fonseca"
__email__ = "fonseca.rasmus@gmail.com"
