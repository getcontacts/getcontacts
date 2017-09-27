from collections import defaultdict
import sys
import os
import re
import json
from utils import *

def compile_freqs(metadata, output_directory):
	interaction_database = defaultdict(dict)
	for descriptor in metadata:
		# print "Descriptor: %s" % descriptor
		descriptor_dir = metadata[descriptor]
		# print "Descriptor directory: %s" % descriptor_dir
		for filename in listdir_files(descriptor_dir):
			# print "Filename: %s" % filename
			if re.match(".*_frequencies_generic\.csv$", filename) and os.path.basename(filename)[0] != '.':
				interaction_type = os.path.basename(filename).split('_')[0]
				# print "Interaction type: %s" % interaction_type
				with open(filename, 'r') as ropen:
					lines = [line.strip().split(',') for line in ropen.readlines()]
				header = lines[0]
				for interaction_line in lines[1:]:
					res1 = interaction_line[0]
					res2 = interaction_line[1]
					frequency = float(interaction_line[2])
					interaction = '_'.join([res1, res2])
					# print "Interaction: %s" % interaction
					if interaction not in interaction_database[interaction_type]:
						interaction_database[interaction_type][interaction] = {}
					interaction_database[interaction_type][interaction][descriptor] = frequency

	for interaction_type in interaction_database:
		output_filename = "%s%s.csv" % (output_directory, interaction_type)
		with open(output_filename, 'w+') as wopen:
			wopen.write("%s\n" % (','.join(["ResPair"] + sorted(metadata.keys()))))
			for interaction in interaction_database[interaction_type]:
				writeLine = False
				output_line = [interaction]
				for descriptor in sorted(metadata.keys()):
					if descriptor in interaction_database[interaction_type][interaction]:
						frequency = interaction_database[interaction_type][interaction][descriptor]
						if frequency >= 0.5:
							writeLine = True
						output_line.append("%.8f" % frequency)
					else:
						output_line.append("%.8f" % 0.0)
				if writeLine:
					wopen.write("%s\n" % ','.join(output_line))

def compile_simulations(argv):
	output_directory = argv[1]
	metadata_filename = argv[2]
	open_dir(output_directory)

	with open(metadata_filename, 'r') as ropen:
		metadata_lines = [line.strip() for line in ropen.readlines()]

	metadata = {}
	for line in metadata_lines:
		if len(line) == 0:
			continue
		descriptor, descriptor_dir = line.split(',')
		metadata[descriptor] = descriptor_dir

	compile_freqs(metadata, output_directory)

if __name__ == "__main__":
	compile_simulations(sys.argv)
