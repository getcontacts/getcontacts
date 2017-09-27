from collections import defaultdict
import sys
import os
import errno
import re
import json
from utils import *

def compile_interaction_files(outputs_filename, interaction_files):
	descriptor_database = defaultdict(dict)
	interaction_database = set()
	for filename in interaction_files:
		interaction_type = os.path.basename(filename).split('.')[0]
		if interaction_type in ["vdw", "hbbb", "wb"]:
			continue
		with open(filename, 'r') as ropen:
			file_lines = [line.strip() for line in ropen.readlines()]
		header = file_lines[0]
		descriptors = header.split(',')[1:]
		for line in file_lines[1:]:
			line = line.split(',')
			if len(line) == 0:
				continue
			interaction = '-'.join(line[0:2] + [interaction_type])
			# print interaction_type
			interaction_frequencies = line[1:]
			interaction_database.add(interaction)

			for i in xrange(len(descriptors)):
				descriptor_database[descriptors[i]][interaction] = interaction_frequencies[i]

	with open(outputs_filename, 'w+') as wopen:
		wopen.write("%s\n" % ','.join(["ResPair"] + sorted(descriptor_database.keys())))
		for interaction in interaction_database:
			output_line = [interaction]
			for descriptor in sorted(descriptor_database.keys()):
				if interaction in descriptor_database[descriptor]:
					output_line.append(descriptor_database[descriptor][interaction])
				else:
					output_line.append('%.8f' % 0)
			wopen.write("%s\n" % ','.join(output_line))

def compile_interactions(argv):
	outputs_filename = argv[1]
	interaction_files = argv[2:]
	compile_interaction_files(outputs_filename, interaction_files)

if __name__ == "__main__":
	compile_interactions(sys.argv)
