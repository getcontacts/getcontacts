import sys
import os
import re
import json
from utils import *
from collections import defaultdict

'''Process correctly formatted genericization dictionary'''
def build_database(generic_dict):
	database = defaultdict(dict)
	with open(generic_dict, 'r') as generic_open:
		lines = [line.strip().split(',') for line in generic_open.readlines() if len(line.strip()) > 0]
	protein_names = lines[0][1:]
	residue_names = lines[1:]
	for line in residue_names:
		generic_name = line[0]
		specific_names = line[1:]
		for i in xrange(len(specific_names)):
			specific_name = specific_names[i]
			protein_name = protein_names[i]
			if specific_name != '-':
				database[protein_name][specific_name] = generic_name
	return database

def get_generic_res(old_res, protein):
	if (old_res, protein) in to_omit:
		return "No residue found"
	aa, pos = old_res.split(':')
	if special_case:
		pos = str(int(pos) + 1)
	new_res = seq1[aa] + pos

	if new_res not in database[protein]:
		print "Residue %s in protein %s is not in genericization database, and interactions involving this residue will be omitted in genericized files" % (new_res, protein)
		to_omit.append((old_res, protein))
		return "No residue found"
	return database[protein][new_res]

def convert_json(filename, protein):
	with open(filename, 'r') as to_convert_open:
		d = json.load(to_convert_open)
	new_edges = []
	for interaction in d["edges"]:
		generic_res1 = get_generic_res(interaction["name1"], protein)
		generic_res2 = get_generic_res(interaction["name2"], protein)
		if generic_res1 == "No residue found" or generic_res2 == "No residue found":
			continue
		new_edges.append({"name1":generic_res1, "name2":generic_res2, "frames":interaction["frames"]}) 
	with open(append_to_filename(filename, "generic"), 'w+') as write_file_open:
		write_file_open.write(json.dumps({"edges":new_edges}))

def convert_csv(filename, protein):
	with open(filename, 'r') as to_convert_open:
		file_lines = [line.strip().split(',') for line in to_convert_open.readlines() if len(line) > 0]
	header = file_lines[0]
	new_lines = []
	for line in file_lines[1:]:
		generic_res1 = get_generic_res(line[0], protein)
		generic_res2 = get_generic_res(line[1], protein)
		if generic_res1 == "No residue found" or generic_res2 == "No residue found":
			continue
		line[0] = generic_res1
		line[1] = generic_res2
		new_lines.append(line)
	with open(append_to_filename(filename, "generic"), 'w+') as write_file_open:
		write_file_open.write("%s\n" % ','.join(header))
		for line in new_lines:
			write_file_open.write("%s\n" % ','.join(line))

def genericize(argv):
	#declare global arguments
	global seq1
	global database
	global special_case
	global to_omit

	#grab command line arguments
	dirname = argv[1]
	protein = argv[2]
	generic_dict = argv[3]
	special_case = len(argv) > 4 #kind of hacky, to remove
	to_omit = []

	with open("./utils/seq1.json", 'r') as seq1_open:
		seq1 = json.load(seq1_open) #defines mapping from 3-letter amino acid codes to 1-letter codes
	database = build_database(generic_dict)

	for filename in listdir_files(dirname):
		if re.match(".*_byres\.json$", filename):
			convert_json(filename, protein)
			print "%s genericized successfully." % filename
		elif re.match(".*_frequencies\.csv$", filename):
			convert_csv(filename, protein)
			print "%s genericized successfully." % filename

if __name__ == "__main__":
	genericize(sys.argv)
