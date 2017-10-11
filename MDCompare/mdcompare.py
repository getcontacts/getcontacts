'''
The MDCompare utility enables users to compare the frequencies of interactions across multiple MDContactNetworks
outputs. 
'''

from __future__ import division
import sys
import os
import argparse
import re
import json
from utils.utils import *

def get_write_lines(respair_to_simcond_to_frequency, simulation_conditions):
	write_lines = []
	for respair in respair_to_simcond_to_frequency:
		write_line = ['-'.join(respair)]
		should_write = False
		for simcond in simulation_conditions:
			if simcond in respair_to_simcond_to_frequency[respair]:
				frequency = respair_to_simcond_to_frequency[respair][simcond][0] / respair_to_simcond_to_frequency[respair][simcond][1]
				write_line.append("%.4f" % frequency)
				if frequency >= 0.5:
					should_write = True
			else:
				write_line.append("%.4f" % 0.0)
		if should_write:
			write_lines.append(write_line)
	return write_lines

def weight_multiple_IDs(simcond_data_list):
	inttype_to_frequencies = {}

	# Find the largest common set of interaction types present in all datasets for this simulation condition
	inttype_set = set(simcond_data_list[0].keys())
	for ID_data in simcond_data_list:
		inttype_set = inttype_set.intersection(set(ID_data.keys()))
	inttype_list = sorted(list(inttype_set))

	for inttype in inttype_list:
		respair_to_frequency = {}
		for i in xrange(len(simcond_data_list)):
			ID_inttype_to_frequencies = simcond_data_list[i]
			inttype_specific_ID_inttype_to_frequencies = ID_inttype_to_frequencies[inttype]
			for respair in inttype_specific_ID_inttype_to_frequencies:
				numFrames, totFrames = inttype_specific_ID_inttype_to_frequencies[respair]
				if respair in respair_to_frequency:
					old_numFrames, old_totFrames = respair_to_frequency[respair]
					new_numFrames = numFrames + old_numFrames
					new_totFrames = totFrames + old_totFrames
					respair_to_frequency[respair] = (new_numFrames, new_totFrames)
				else:
					respair_to_frequency[respair] = (numFrames, totFrames)
		inttype_to_frequencies[inttype] = respair_to_frequency
	return inttype_to_frequencies

def genericize_res(old_res, protein):
''' Given a residue in the format 'Ala:379', will convert the residue name to a generic name
using the user-provided generic_dict
'''
	global generic_database

	aa, pos = old_res.split(':')
	new_res = seq1[aa] + pos
	if new_res not in generic_database[protein]:
		print "Residue %s in protein %s is not in genericization database" % (new_res, protein)
		return None
	else:
		return generic_database[protein][new_res]

def genericize(inttype_to_frequencies, protein):
	new_inttype_to_frequencies = {}
	for inttype in inttype_to_frequencies:
		new_respair_to_frequency = {}
		for respair in inttype_to_frequencies[inttype]:
			res1, res2 = respair
			new_res1, new_res2 = genericize_res(res1, protein), genericize_res(res2, protein)
			if not new_res1 or not new_res2:
				continue
			new_respair_to_frequency[(new_res1, new_res2)] = inttype_to_frequencies[inttype][respair]
		new_inttype_to_frequencies[inttype] = new_respair_to_frequency
	return new_inttype_to_frequencies


def build_generic_database(generic_dict):
	generic_database = {}
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
			if protein_name not in generic_database:
				generic_database[protein_name] = {}
			if specific_name != '-':
				generic_database[protein_name][specific_name] = generic_name
	return generic_database

def extract_data(path):
	global first_update_to_interactions_to_study
	global interactions_to_study

	inttype_to_frequencies = {}
	for filename in listdir_files(path):
		if not re.match(".*_frequencies\.csv", os.path.basename(filename)):
			continue
		inttype, _ = os.path.basename(filename).split('_')
		assert(inttype not in inttype_to_frequencies)
		with open(filename, 'r') as ropen:
			file_lines = [line.strip().split(',') for line in ropen.readlines() if len(line.strip()) > 0]
		respair_to_frequency = {}
		for line in file_lines[1:]:
			res1 = line[0]
			res2 = line[1]
			numFrames = int(line[3])
			totFrames = int(line[4])
			respair_to_frequency[(res1, res2)] = (numFrames, totFrames)
		inttype_to_frequencies[inttype] = respair_to_frequency
	if first_update_to_interactions_to_study:
		interactions_to_study = set(inttype_to_frequencies.keys())
		first_update_to_interactions_to_study = False
	else:
		interactions_to_study = interactions_to_study.intersection(set(inttype_to_frequencies.keys()))

	return inttype_to_frequencies

def extract_input_file(input_file):
	simcond_to_id = {}
	id_to_path = {}
	simcond_to_protein = {}

	with open(input_file, 'r') as ropen:
		matrix_lines = [line.strip() for line in ropen.readlines() if len(line.strip()) > 0]
	for line in matrix_lines[1:]:
		ID, simcond, path, protein = line.split(',')
		ID = int(ID)
		if simcond in simcond_to_id:
			assert(ID not in simcond_to_id[simcond])
			simcond_to_id[simcond].append(ID)
		else:
			simcond_to_id[simcond] = [ID]

		id_to_path[ID] = clean_path(path)

		if simcond in simcond_to_protein:
			assert(simcond_to_protein[simcond] == protein)
		else:
			simcond_to_protein[simcond] = protein
	return simcond_to_id, id_to_path, simcond_to_protein

def mdcompare(argv):
	# Define globals
	global generic_database
	global seq1
	global interactions_to_study
	global first_update_to_interactions_to_study
	first_update_to_interactions_to_study = True

	# Parse command line arguments
	parser = argparse.ArgumentParser(description='MDCompare companion to MDContactNetworks')
	parser.add_argument('-i', dest='input_file', required=True, nargs = 1, help="Correctly formatted .csv file of input")
	parser.add_argument('-o', dest='output_directory', required=True, nargs = 1, help="Directory for MDCompare outputs")
#	parser.add_argument('-s', dest='interactions_to_study', required=True, nargs = 1, help="Quotation-enclosed list of interaction flags")
	parser.add_argument('-g', dest='generic_dict', nargs='?', help="A correctly formatted file for standardizing residue names between different proteins.")
	results = parser.parse_args()
	input_file = results.input_file[0]
	output_directory = clean_path(results.output_directory[0])
#	interactions_to_study = results.interactions_to_study
	generic_dict = results.generic_dict

	open_dir(output_directory)
	simcond_to_id, id_to_path, simcond_to_protein = extract_input_file(input_file)
	should_genericize = generic_dict is not None
	if len(set(simcond_to_protein.values())) > 1 and not should_genericize:
		assert("Generic dict must be provided because more than one protein appears in input file")
	if should_genericize:
		generic_database = build_generic_database(generic_dict)
		# seq1 is a dict from 3-letter amino acid codes to 1-letter codes
		with open("./utils/seq1.json", 'r') as seq1_open:
			seq1 = json.load(seq1_open)			

	# Construct a nested dict from simulation_condition -> interaction type -> residue pair -> tuple of
	# (# of frames with this interaction, total frames)
	simcond_to_data = {}
	for simcond in simcond_to_id:
		simcond_data_list = []
		for ID in simcond_to_id[simcond]:
			path = id_to_path[ID]
			ID_data = extract_data(path)
			if should_genericize:
				protein = simcond_to_protein[simcond]
				ID_data = genericize(ID_data, protein)
			simcond_data_list.append(ID_data)
		# If more than one ID correspond to this simulation condition, sum the numFrames and totFrames of
		# every interaction type->residue pair between its appearances in the data from different ID's
		if len(simcond_data_list) == 1:
			simcond_to_data[simcond] = simcond_data_list[0]
		else:
			simcond_to_data[simcond] = weight_multiple_IDs(simcond_data_list)

	simulation_conditions = sorted(simcond_to_data.keys())
	inttype_to_respair_to_simcond_to_frequency = {}
	for inttype in interactions_to_study:
		respair_to_simcond_to_frequency = {}
		for simcond in simcond_to_data:
			for respair in simcond_to_data[simcond][inttype]:
				if respair not in respair_to_simcond_to_frequency:
					respair_to_simcond_to_frequency[respair] = {}
				respair_to_simcond_to_frequency[respair][simcond] = simcond_to_data[simcond][inttype][respair]

		inttype_to_respair_to_simcond_to_frequency[inttype] = respair_to_simcond_to_frequency

		write_lines = get_write_lines(respair_to_simcond_to_frequency, simulation_conditions)
		with open("%s%s.csv" % (output_directory, inttype), 'w+') as wopen:
			wopen.write("%s\n" % ','.join(["Residue Pair"] + simulation_conditions))
			for write_line in write_lines:
				wopen.write("%s\n" % ','.join(write_line))

	with open("%scompiled_interactions.csv" % output_directory, 'w+') as wopen:
		wopen.write("%s\n" % ','.join(["Residue Pair"] + simulation_conditions))
		for inttype in interactions_to_study:
			if inttype in ['hbbb', 'wb', 'wb2', 'vdw']:
				continue
			else:
				write_lines = get_write_lines(inttype_to_respair_to_simcond_to_frequency[inttype], simulation_conditions)
				for write_line in write_lines:
					write_line[0] += "-%s" % inttype
					wopen.write("%s\n" % ','.join(write_line))



if __name__ == '__main__':
	mdcompare(sys.argv)
