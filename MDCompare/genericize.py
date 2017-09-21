import sys
import os
import re
import json
from utils import *
from collections import defaultdict

def build_database(generic_dict):
	database = defaultdict(dict)
	with open(generic_dict, 'r') as generic_open:
		lines = [line.strip().split(',') for line in generic_open.readlines()]
	generic_names = lines[0][1:]
	protein_names = lines[1:]
	for line in protein_names:
		protein = line[0]
		res_names = line[1:]
		for i in xrange(len(res_names)):
			res_name = res_names[i]
			if res_name is not '-':
				database[protein][res_name] = generic_names[i]
	return database

def get_generic_res(old_res, protein):
	aa, pos = old_res.split(':')
	if special_case:
		pos = str(int(pos) + 1)
	new_res = seq1[aa] + pos

	if new_res not in database[protein]:
		print "Residue %s in protein %s is not in genericization database, and interactions involving this residue will be omitted in genericized files" % (new_res, protein)
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
	with open(prepend_to_filename(filename, "generic"), 'w+') as write_file_open:
		write_file_open.write(json.dumps({"edges":new_edges}))

def convert_csv(filename, protein):
	with open(filename, 'r') as to_convert_open:
		file_lines = [line.strip().split(',') for line in to_convert_open.readlines() if len(line) > 0]
	header = file_lines[0]
	new_lines = []
	for line in file_lines[1:]:
		# print line
		line[0] = get_generic_res(line[0], protein)
		line[1] = get_generic_res(line[1], protein)
		new_lines.append(line)
	with open(prepend_to_filename(filename, "generic"), 'w+') as write_file_open:
		write_file_open.write("%s\n" % ','.join(header))
		for line in new_lines:
			write_file_open.write("%s\n" % ','.join(line))

def genericize(argv):
	global seq1
	global database
	global special_case

	dirname = argv[1]
	protein = argv[2]
	generic_dict = argv[3]
	special_case = len(argv) > 4

	with open("./utils/seq1.json", 'r') as seq1_open:
		seq1 = json.load(seq1_open)
	database = build_database(generic_dict)

	for filename in listdir_files(dirname):
		if re.match(".*_byres\.json$", filename):
			convert_json(filename, protein)
		elif re.match(".*\.csv$", filename) and not re.match(".*_generic\.csv$", filename):
			convert_csv(filename, protein)
		print "%s genericized successfully." % filename

if __name__ == "__main__":
	genericize(sys.argv)
