"""
Draws a structure and it's atomic interactions in pymol.

Requires at least a structure file and a contact-file as input. Additionally, a pymol selection can be provided as well.
Note that "--" must be added for pymol to pass command line arguments correctly (see examples)

Example usages:
    pymol pymol_frequencies.py -- ../example/5xnd_topology.pdb ../example/5xnd_all-contacts.tsv

    pymol pymol_frequencies.py -- ../example/5xnd_topology.pdb ../example/5xnd_all-contacts.tsv "resi 22-30"
"""

import pymol
pymol.finish_launching()

import sys
import re
from collections import defaultdict

print sys.argv

if len(sys.argv) not in [3,4]:
    print "Usage: pymol "+sys.argv[0]+" -- <structurefile> <contactfile> [selection]"
    cmd.quit()

# Ensure that optional selection is set to "all" by default
if len(sys.argv) == 3:
    sys.argv.append("all")

cmd.load(sys.argv[1])


interaction_frames = defaultdict(set)
total_frames = 0
max_frame = 0
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

        if max_frame < frame:
            max_frame = frame
        interaction_frames[(atom1, atom2)].add(frame)

total_frames = max(total_frames, max_frame + 1)


# Write contacts
cgos = []
for (atom1, atom2), frames in interaction_frames.items():
    frequency = len(frames) / float(total_frames)
    a1chain, _, a1resi, a1name = atom1.split(":")
    a2chain, _, a2resi, a2name = atom2.split(":")
    pmatom1 = "/".join(["","","",a1chain, a1resi, a1name])
    pmatom2 = "/".join(["","","",a2chain, a2resi, a2name])
    c1 = cmd.get_model(pmatom1).atom[0].coord
    c2 = cmd.get_model(pmatom2).atom[0].coord
    if len(cmd.get_model(pmatom1 + " & "+sys.argv[3]).atom)==0 and len(cmd.get_model(pmatom2 + " & "+sys.argv[3]).atom)==0:
        continue
    # c2 = cmd.get_model(pmatom2).atom[0].coord
    rad = frequency * 0.18
    col = [0.2+0.6*frequency, 0.2+0.4*frequency, 0.2+0.7*frequency]

    cgos += [CYLINDER] + c1 + c2 + [rad] + col + col

cmd.load_cgo(cgos, "contacts")


