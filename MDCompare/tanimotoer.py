from __future__ import division
import sys
import json
import re
# import pandas as pd
# import seaborn as sns
import os
from utils import *
import subprocess

def produce_tanimoto(filename):
	if re.search("vdw", filename):
		return
	print "Producing tanimoto file for: %s" % filename
	with open(filename, 'r') as json_open:
		d = json.load(json_open)
	output_filename = os.path.splitext(filename)[0] + "_tanimoto.csv"
	print "Output filename: %s" % output_filename

	#Make a dict from a residue-pair to the set of frames that the residue-pair interaction appears in
	interaction_to_frames = {}
	for edge in d["edges"]:
		frames_set = set(edge["frames"])
		interaction = "%s-%s-%.2f" % (edge["name1"], edge["name2"], len(frames_set)/1500)
		if len(frames_set) > 1500*0.5 or len(frames_set) < 1500*0.2:
			continue
		interaction_to_frames[interaction] = frames_set

	if len(interaction_to_frames) == 0:
		return

	#Loop through all residue-pairs pairwise and generate tanimoto score
	sorted_interactions = sorted(interaction_to_frames.keys())
	with open(output_filename, 'w+') as wopen:
		wopen.write("%s\n" % ','.join([""] + sorted_interactions))
		for interaction1 in sorted_interactions:
			write_line = [interaction1]
			frames1 = interaction_to_frames[interaction1]
			for interaction2 in sorted_interactions:
				frames2 = interaction_to_frames[interaction2]
				tanimoto = len(frames1.intersection(frames2)) / len(frames1.union(frames2))
				write_line.append("%.8f" % tanimoto)
			wopen.write("%s\n" % ','.join(write_line))

def main(argv):
	paths = []

	with open("./sample_inputs/sample_input.csv", 'r') as ropen:
		matrix_lines = [line.strip() for line in ropen.readlines()]
	for line in matrix_lines:
		if len(line) == 0:
			continue
		descriptor, path, _ = line.split(',')
		if re.search(argv[1], descriptor):
			paths.append(clean_path(path))

	for path in paths:
		for filename in listdir_files(path):
			if re.search("_byatom\.json", filename):
				produce_tanimoto(filename)
			if re.search("_byres_generic\.json", filename):
				produce_tanimoto(filename)

	

if __name__ == "__main__":
	main(sys.argv)
