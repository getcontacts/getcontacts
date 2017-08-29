from collections import defaultdict

def get_residue_from_atom(atom_name):
	split_atom_name = atom_name.split(":")
	return ":".join(split_atom_name[0:3])

def parse_stitched_file(stitched_filename):
	with open(stitched_filename) as f:
		file_lines = f.readlines()
	file_lines_split = [x.strip().split() for x in file_lines] 
	return file_lines_split

def get_interaction_set(stitched):
	interaction_to_frames = defaultdict(set)
	for instance in stitched:
		frame = instance[0]
		atom1 = instance[1]
		atom2 = instance[2]
		interaction_pair = get_residue_from_atom(atom1) + "\t" + get_residue_from_atom(atom2)
		interaction_type = instance[3]
		interaction_to_frames[interaction_pair].add(frame)
	return interaction_to_frames

def create_frequencies_file(itype, OUTPUT_DIR, stitched_filename, sim_length, interaction_to_frames):
	frequencies_filename = OUTPUT_DIR + "/" + itype + "_frequencies.txt" 
	fq = open(frequencies_filename, 'w')

	fq.write('NumFrames:\t%d\n' % (sim_length))

	for interaction in interaction_to_frames:
		interaction_frequency = float(len(interaction_to_frames[interaction])) / sim_length
		output_string = "%s\t%.4f\t%d\n" % (interaction, interaction_frequency, len(interaction_to_frames[interaction]))
		fq.write(output_string)

	fq.close()

def make_frequencies_file(itype, OUTPUT_DIR, stitched_filename, sim_length):
	stitched = parse_stitched_file(stitched_filename)
	interaction_to_frames = get_interaction_set(stitched)

	create_frequencies_file(itype, OUTPUT_DIR, stitched_filename, sim_length, interaction_to_frames)

