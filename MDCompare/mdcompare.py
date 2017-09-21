import sys
import re
import subprocess
import argparse
from collections import defaultdict
from utils import *

def extract_databases(input_matrix):
	paths_database = defaultdict(list)
	proteins_database = {}
	with open(input_matrix, 'r') as ropen:
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
	parser = argparse.ArgumentParser(description='MDCompare companion to MDContactNetworks')
	parser.add_argument('input_matrix', nargs = 1, help="Correctly formatted .csv file of input")
	parser.add_argument('output_directory', nargs = 1, help="Directory for MDCompare outputs")
	parser.add_argument('-generic_dict', '-gd', dest='generic_dict', nargs='?', help="A correctly formatted file for standardizing residue names between different. A genericization dict must be provided if 2 or more proteins appear in the input matrix")
	results = parser.parse_args()
	input_matrix = results.input_matrix[0]
	output_directory = clean_path(results.output_directory[0])
	should_genericize = results.generic_dict is not None
	open_dir(output_directory)
	paths_database, proteins_database = extract_databases(input_matrix)

	if should_genericize:
		generic_dict = results.generic_dict
		# genericization_processes = []
		for descriptor in paths_database:
			for path in paths_database[descriptor]:
				command = ["python", "genericize.py", path, proteins_database[descriptor], generic_dict]
				if re.search("M3_Inactive", descriptor):
					command.append("0")
				sp = subprocess.Popen(command)
				sp.wait()
		print "All files genericized!"

	# weighting_processes = []
	# for descriptor in paths_database:
	# 	weighting_processes.append(weight_runs())
	# 	for weighting_process in weighting_processes:
	# 		weighting_process.wait()

	# compare_frequencies()
	# compare_interactions()
	# generate_heatmaps()




if __name__ == '__main__':
	mdcompare(sys.argv)
