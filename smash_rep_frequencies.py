from collections import defaultdict
import sys
import os
import errno
sys.path.insert(0, '/scratch/PI/rondror/augustine/')
from run_script import *
import re

def open_dir(dir_name):
	try:
		os.makedirs(dir_name)
	except OSError as e:
		if e.errno != errno.EEXIST:
			raise

def clean_path(path):
	if path[-1] != '/':
		path += '/'
	return path

def listdir_files(path):
	return [path + filename for filename in os.listdir(path) if os.path.isfile(path + filename)]

def compile_interaction_reps(interaction_type, rep_paths):
	interaction_stats = dict()
	numFrames = 0

	for rep_path in rep_paths:
		freq_filename = rep_path + interaction_type + '_' + 'frequencies.txt'

#		if freq_filename != '/scratch/PI/rondror/augustine/MuscarinicComparisons/outputs/m4_active_ach/1/sb_frequencies.txt':
#			continue

		with open(freq_filename, 'r') as open_freq_file:
			interaction_lines = [line.strip().split() for line in open_freq_file.readlines()]
			for interaction_line in interaction_lines:
#				print interaction_line
				if interaction_line[0] == 'NumFrames:':
					numFrames += float(interaction_line[1])
				else:
					interaction = '\t'.join(interaction_line[:2])
					freq = float(interaction_line[3])

					if interaction in interaction_stats:
						interaction_stats[interaction] += freq
					else:
						interaction_stats[interaction] = freq

	return numFrames, interaction_stats

#	print interaction_stats

def read_path_files(rep_paths, outputs_path, interaction_types):
	receptors = {'m1': 'ACM1_HUMAN', 'm2': 'ACM2_HUMAN', 'm3': 'ACM3_HUMAN', 'm4': 'ACM4_HUMAN'}
	states = ['active', 'inactive']

	for interaction_type in interaction_types:
		for state in states:
			open_dir('%s%s/%s/' % (outputs_path, interaction_type, state))

	compiled_filenames = []

	for filename in listdir_files(rep_paths):
		file_descriptor = os.path.splitext(os.path.basename(filename))[0].split('_')
		receptor = receptors[file_descriptor[0]]
		state = file_descriptor[1]
		print file_descriptor, receptor, state
		'''
		if os.path.basename(filename) != 'm1_active_ach.txt':
			continue
			'''
		with open(filename, 'r') as read_file:
			rep_paths = [clean_path(line.strip()) for line in read_file.readlines()]
		for interaction_type in interaction_types:
			numFrames, interaction_stats = compile_interaction_reps(interaction_type, rep_paths)
			compiled_filename = '%s%s/%s/%s' % (outputs_path, interaction_type, state, os.path.basename(filename))
			compiled_filenames.append([compiled_filename, receptor])

#			print compiled_filename

			with open(compiled_filename, 'w+') as int_file_open:
#				int_file_open.write('numFrames:\t%d\n' % (numFrames))
				for interaction in interaction_stats:
					int_file_open.write('%s\t%.8f\n' % (interaction, interaction_stats[interaction] / numFrames))

	with open('%sto_convert.txt' % (outputs_path), 'w+') as to_convert_open:
		for compiled_filename in compiled_filenames:
			to_convert_open.write('%s\n' % ('\t'.join(compiled_filename)))

def main(argv):
	rep_paths = clean_path(argv[1])
	outputs_path = clean_path(argv[2])
	interaction_types = [int_type[1:] for int_type in argv[3:]]
	read_path_files(rep_paths, outputs_path, interaction_types)

if __name__ == "__main__":
	main(sys.argv)
