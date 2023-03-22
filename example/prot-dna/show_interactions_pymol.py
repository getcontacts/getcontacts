from collections import defaultdict

sele = "all"
struc = "6cvo.pdb"
conts = "contacts.tsv"

# Read structure
cmd.load(struc)

# Parse contact-file
interaction_frames = defaultdict(set)
total_frames = 0
with open(conts) as cfile:
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
selection = sele
cgos = []
for (atom1, atom2), frames in interaction_frames.items():
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

    rad = 0.10 + 0.05
    col = [0.83, 0.33, 0.32]
    cgos += [CYLINDER] + c1 + c2 + [rad] + col + col

cmd.load_cgo(cgos, "contacts")
cmd.color("gray80", "elem C")
cmd.orient()


