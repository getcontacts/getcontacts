#!/bin/sh
# Shebang-hack for launching pymol
''':'
exec pymol -q "$0" -- "$@"
'''

__doc__ = """
Draws a structure and it's atomic interactions in pymol. Requires at least a structure file 
and a contact-file as input. Additionally, a pymol selection can be provided as well.

Example usages:
    pymol_frequencies.py ../example/5xnd_topology.pdb ../example/5xnd_all-contacts.tsv
    pymol_frequencies.py ../example/5xnd_topology.pdb ../example/5xnd_all-contacts.tsv "resi 22-30"

    # Alternatively pymol can also be called directly on this script using the "--" argument
    pymol pymol_frequencies.py -- ../example/5xnd_topology.pdb ../example/5xnd_all-contacts.tsv
"""

import re
from collections import defaultdict
import sys

# Check cmd-line arguments
if len(sys.argv) not in [3,4] or "pymol" not in sys.modules:
    print("Usage: "+sys.argv[0]+" <structurefile> <contactfile> [selection]")
    print("or:    pymol "+sys.argv[0]+" -- <structurefile> <contactfile> [selection]")
    print(__doc__)
    sys.exit(1)

# Ensure that optional selection is set to "all" by default
if len(sys.argv) == 3:
    sys.argv.append("all")

# Read structure
cmd.load(sys.argv[1])

# Parse contact-file
interaction_frames = defaultdict(set)
total_frames = 0
with open(sys.argv[2]) as cfile:
    for line in cfile.readlines():
        line = line.strip()
        if not line:  # Ignore empty lines
            continue

        # Parse header
        if line[0] == "#":
            total_frame_match = re.search(r'total_frames:(\d+)', line)
            if total_frame_match:
                total_frames = int(total_frame_match.group(1))
            continue

        # Regex that matches a contact line and groups the frame number and the two first residues,
        frame_atom_match = re.match(r'^(\d+)\t.*?\t([^:]+:[^:]+:\d+:[^\s]+)\t([^:]+:[^:]+:\d+:[^\s]+)', line)
        frame = int(frame_atom_match.group(1))
        atom1 = frame_atom_match.group(2)
        atom2 = frame_atom_match.group(3)

        if atom2 < atom1:
            atom1, atom2 = atom2, atom1

        interaction_frames[(atom1, atom2)].add(frame)



# Write contacts
selection = sys.argv[3]
cgos = [[],[],[]]
for (atom1, atom2), frames in interaction_frames.items():
    frequency = len(frames) / float(total_frames)
    a1chain, _, a1resi, a1name = atom1.split(":")
    a2chain, _, a2resi, a2name = atom2.split(":")
    pmatom1 = "/".join(["","","",a1chain, a1resi, a1name])
    pmatom2 = "/".join(["","","",a2chain, a2resi, a2name])
    c1 = cmd.get_model(pmatom1).atom[0].coord
    c2 = cmd.get_model(pmatom2).atom[0].coord
  
    # Check that either pmatom1 or pmatom2 are in the selection
    if len(cmd.get_model(pmatom1 + " & "+selection).atom)==0 and \
       len(cmd.get_model(pmatom2 + " & "+selection).atom)==0:
        continue

    rad = frequency * 0.10 + 0.05
    if frequency > 0.75:
        col = [0.19, 0.56, 0.18]
        cgo_idx = 0
    elif frequency > 0.25:
        col = [1.0, 0.89, 0.35]
        cgo_idx = 1
    else:
        col = [0.83, 0.33, 0.32]
        cgo_idx = 2

    cgos[cgo_idx] += [CYLINDER] + c1 + c2 + [rad] + col + col

cmd.load_cgo(cgos[0], "high_freq")
cmd.load_cgo(cgos[1], "medium_freq")
cmd.load_cgo(cgos[2], "low_freq")
cmd.color("gray30", "elem C")
cmd.bg_color("white")
cmd.orient()


