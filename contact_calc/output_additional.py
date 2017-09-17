from __future__ import division
from collections import defaultdict
import json

def open_dir(dirname):
	try:
		os.makedirs(dirname)
	except OSError as e:
		if e.errno != errno.EEXIST:
			raise

def get_residue_from_atom(atom_name):
	split_atom_name = atom_name.split(':')
	return "%s:%s" % (split_atom_name[1], split_atom_name[2])

def create_dynamic_jsons(output_filename, interaction_to_frames):
	#assemble and write out .json
	edges = []
	for interaction_pair in interaction_to_frames: #loop over all respairs
		int1, int2 = interaction_pair.split(',')
		edges.append({"name1":int1, "name2":int2, "frames":list(interaction_to_frames[interaction_pair])})
	with open(output_filename, 'w+') as fq:
		fq.write(json.dumps({"edges":edges}))

def create_frequencies_file(frequency_filename, sim_length, respair_to_frames):
	#write out .csv file
	with open(frequency_filename, 'w+') as fq:
		fq.write('Res1,Res2,Freq,NumFrames,,TotalFrames:,%d\n' % sim_length)
		for interaction in respair_to_frames:
			interaction_frequency = len(respair_to_frames[interaction])/sim_length
			fq.write("%s,%.8f,%d\n" % (interaction, interaction_frequency, len(respair_to_frames[interaction])))

def get_atompair_set(itype, stitched_lines):
	atompair_to_frames = defaultdict(set)
	for stitched_line in stitched_lines:
		frame = int(stitched_line[0])
		atom1 = stitched_line[1]
		atom2 = stitched_line[-2]
		atom_pair = ','.join(sorted([atom1, atom2]))
		atompair_to_frames[atom_pair].add(frame)
	return atompair_to_frames

def get_respair_set(itype, stitched_lines):
	respair_to_frames = defaultdict(set)
	for stitched_line in stitched_lines:
		frame = int(stitched_line[0])
		atom1 = stitched_line[1]
		atom2 = stitched_line[-2]
		res_pair = ','.join(sorted([get_residue_from_atom(atom1), get_residue_from_atom(atom2)]))
		respair_to_frames[res_pair].add(frame)
	return respair_to_frames

def make_additional_files(itype, output_dir, stitched_filename, sim_length):
	with open(stitched_filename) as stitched_open:
		stitched_lines = [line.strip().split() for line in stitched_open.readlines()] 

	respair_to_frames = get_respair_set(itype, stitched_lines)
	atompair_to_frames = get_atompair_set(itype, stitched_lines)

	frequency_filename = output_dir + '/' + itype + "_frequencies.csv"
	byres_filename = output_dir + '/' + itype + "_byres.json"
	byatom_filename = output_dir + '/' + itype + "_byatom.json"

	create_frequencies_file(frequency_filename, sim_length, respair_to_frames)
	create_dynamic_jsons(byres_filename, respair_to_frames)
	create_dynamic_jsons(byatom_filename, atompair_to_frames)

