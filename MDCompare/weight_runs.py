from __future__ import division
import sys
import os
import re
import json
from collections import defaultdict
from utils import *

def compile_interaction_reps(basename, repitition_paths):
	totFrames = 0 #counter for total number of frames across repititions
	interaction_stats = defaultdict(int) #counter for instances of an interaction across all repititions

	#Open each repitition folder and find the .txt file for this interaction type
	for repitition_path in repitition_paths:
		freq_filename = "%s%s" % (repitition_path, basename)
		with open(freq_filename, 'r') as open_freq_file:
			interaction_lines = [line.strip().split(',') for line in open_freq_file.readlines()]
		header = interaction_lines[0]
		totFrames += int(header[-1].split(':')[1])
		for interaction_line in interaction_lines[1:]:
			respair = ','.join(interaction_line[0:2])
			interaction_stats[respair] += int(interaction_line[3])
	return totFrames, interaction_stats

def smash_raw_outputs(output_dir, repitition_paths):
	search_pattern = ".*_frequencies\.csv$"
	for filename in os.listdir(repitition_paths[0]):
		if re.match(".*_frequencies_generic\.csv$", filename):
			search_pattern = ".*_frequencies_generic\.csv$"
			break
	new_filenames = [output_dir + filename for filename in os.listdir(repitition_paths[0]) if re.match(search_pattern, filename) and not filename[0] == '.']
	for new_filename in new_filenames:
		totFrames, interaction_stats = compile_interaction_reps(os.path.basename(new_filename), repitition_paths)
		with open(new_filename, 'w+') as wopen:
			wopen.write('Res1,Res2,Freq,NumFrames,TotalFrames:%d\n' % totFrames)
			for respair in interaction_stats:
				wopen.write('%s,%.8f,%d,%d\n' % (respair, interaction_stats[respair]/totFrames, interaction_stats[respair], totFrames))

def main(argv):
	output_dir = clean_path(argv[1])
	repitition_paths = [clean_path(arg) for arg in argv[2:]]
	open_dir(output_dir)
	smash_raw_outputs(output_dir, repitition_paths)

if __name__ == "__main__":
	main(sys.argv)
