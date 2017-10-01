'''
> python mdcompare.py input_file output_directory -gd generic_dict
'''

import sys
import re
import subprocess
import argparse
from collections import defaultdict
from utils import *

def extract_databases(input_file):
	paths_database = defaultdict(list)
	proteins_database = {}
	with open(input_file, 'r') as ropen:
		matrix_lines = [line.strip() for line in ropen.readlines()]
	for line in matrix_lines:
		if len(line) == 0:
			continue
		descriptor, path, protein = line.split(',')
		# print "Descriptor: %s, Path %s, Protein: %s" % (descriptor, path, protein)
		paths_database[descriptor].append(path)
		proteins_database[descriptor] = protein
	return paths_database, proteins_database

def mdcompare(argv):
	#Get and parse command line arguments
	parser = argparse.ArgumentParser(description='MDCompare companion to MDContactNetworks')
	parser.add_argument('input_file', nargs = 1, help="Correctly formatted .csv file of input")
	parser.add_argument('output_directory', nargs = 1, help="Directory for MDCompare outputs")
	parser.add_argument('-generic_dict', '-gd', dest='generic_dict', nargs='?', help="A correctly formatted file for standardizing residue names between different. A genericization dict must be provided if 2 or more proteins appear in the input file")
	results = parser.parse_args()
	input_file = results.input_file[0]
	output_directory = clean_path(results.output_directory[0])

	open_dir(output_directory) #open output directory if it doesn't already exist
	paths_database, proteins_database = extract_databases(input_file)
	should_genericize = len(set(proteins_database.values())) > 1 
	if should_genericize and results.generic_dict is None:
		assert "Input file has more than one protein, so a genericization file must be provided"

	#Genericize all files
	if should_genericize:
		generic_dict = results.generic_dict
		for descriptor in paths_database:
			print descriptor
			for path in paths_database[descriptor]:
				command = ["python", "genericize.py", path, proteins_database[descriptor], generic_dict]

				#ad hoc lines to deal with incorrectly numbered M3 case --- remove for all other applications
				if re.search("ACM3_HUMAN-Inactive", descriptor):
					print descriptor, "Special Case"
					command.append("0")

				sp = subprocess.Popen(command)
				sp.wait()
		print "All files genericized!"

	metadata = {}
	for descriptor in paths_database:
		print "Now weighting %s" % descriptor
		repitition_paths = paths_database[descriptor]
		weighted_dirname = clean_path("%sweighted_results/%s/" % (output_directory, descriptor))
		open_dir(weighted_dirname)
		metadata[descriptor] = weighted_dirname
		command = ["python", "weight_runs.py", weighted_dirname] + repitition_paths
		sp = subprocess.Popen(command)
		sp.wait()
	print "All simulations weighted!"

	metadata_filename = "%smetadata.csv" % output_directory
	with open(metadata_filename, 'w+') as wopen:
		for descriptor in metadata:
			wopen.write("%s\n" % ','.join([descriptor, metadata[descriptor]]))

	sp = subprocess.Popen(["python", "compile_simulations.py", output_directory, metadata_filename])
	sp.wait()
	print "Compiled simulation outputs in outputs folder"

	sp = subprocess.Popen(["python", "compile_interactions.py", "%s/compiled_interactions.csv" % output_directory] + listdir_files(output_directory))
	sp.wait()
	print "Produced compiled interactions file"
	print "All Done!"


if __name__ == '__main__':
	mdcompare(sys.argv)
