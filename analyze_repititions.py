'''
Usage: python analyze_repititions.py paths_folder output_path interaction_type

Args:
	1) paths_folder - Each file in this directory is a text file corresponding to a permutation of
		receptor/state/condition. Each such text file contains a list of paths. Each such path leads
		to a directory of the simulation of the receptor/state/condition
	2) output_path - The directory where outputs of MDContactNetworks are routed. For a file in
		paths_folder, and for a repitition paths in the file, the outputs of MDContactNetworks will be
		saved under output_path/file/repitition_path/
	3) interaction_type - The type(s) of interaction you want to study
'''

import sys
import os
import errno
sys.path.insert(0, '/scratch/PI/rondror/augustine/')
from run_script import *
import re

MDContactNetworks_command_parts = ['TOP=', '\nTRAJ=', '\nOUTPUT_DIR=', '\nINTERACTION_TYPES=', '\ncd /scratch/PI/rondror/augustine/MDContactNetworks\n/share/PI/rondror/software/miniconda/bin/python dynamic_contact_networks.py --topology $TOP --trajectory $TRAJ --output_dir $OUTPUT_DIR --cores 12 --stride 1 --itype $INTERACTION_TYPES']

def get_file_descriptor(full_filename):
	return os.path.splitext(os.path.basename(full_filename))[0]

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

'''Given valid parameters, runs MDContactNetworks'''
def run_MDContactNetworks(name, TOP, TRAJ, OUTPUT_DIR, INTERACTION_TYPE, output_path):
	open_dir('%ssbatch_outputs/' % (output_path))
	command = MDContactNetworks_command_parts[0] + TOP + MDContactNetworks_command_parts[1] + TRAJ  \
			+ MDContactNetworks_command_parts[2] + OUTPUT_DIR + MDContactNetworks_command_parts[3] \
			+ INTERACTION_TYPE + MDContactNetworks_command_parts[4]
	output_filename = '%ssbatch_outputs/%s_OUTPUT.out' % (output_path, name)
	error_filename = '%ssbatch_outputs/%s_ERROR.out' % (output_path, name)
	create_script(name, name, output_filename, error_filename, command)

'''Given the path to a MD simulation, extract the filenames of the topology and trajectory files
of the simulation'''
def extract_rep_parameters(repitition_path):
	for md_filename in listdir_files(repitition_path):
		if re.search('\.psf$', md_filename):
			TOP = md_filename
		if re.search('reimaged\.nc$', md_filename) != None:
			TRAJ = md_filename
	assert len(TOP) > 0, "Topology file could not be found in directory %s" % (repitition_path)
	assert len(TRAJ) > 0, "Trajectory file could not be found in directory %s" % (repitition_path)
	return TOP, TRAJ

def read_path_files(paths_folder, output_path, interaction_types):
	#If output directories don't exist, create them
	open_dir(output_path)
	open_dir('%srepitition_paths/' % (output_path))

	#Loop through all files in the paths_folder
	for paths_filename in listdir_files(paths_folder):
		#Get every repitition path
		with open(paths_filename, 'r') as paths_filename_open:
			repitition_paths = [clean_path(line.strip()) for line in paths_filename_open.readlines()]

		#Run MDContactNetworks on each repitition
		output_dirname = clean_path('%s%s/' % (output_path, get_file_descriptor(paths_filename)))
		for repitition_path in repitition_paths:
			TOP, TRAJ = extract_rep_parameters(repitition_path)
			OUTPUT_DIR = clean_path(output_dirname + os.path.basename(repitition_path[:-1]))
			run_name = '%s-%s' % (get_file_descriptor(paths_filename), os.path.basename(repitition_path[:-1]))
			print "Now queueing %s:" % (run_name)
			run_MDContactNetworks(run_name, TOP, TRAJ, OUTPUT_DIR, interaction_types, output_path)

		#Write path to repitition output to corresponding file in the repitition_paths directory
		with open('%srepitition_paths/%s' % (output_path, os.path.basename(paths_filename)), 'w+') as write_file:
			for repitition_path in repitition_paths:
				write_file.write(clean_path(output_dirname + os.path.basename(repitition_path[:-1])) + '\n')

	print "All repititions successfully queued. Check squeue for progress."

def main(argv):
	paths_folder = clean_path(argv[1])
	output_path = clean_path(argv[2])
	interaction_types = '\'' + ' '.join(argv[3:]) + '\''

	read_path_files(paths_folder, output_path, interaction_types)

if __name__ == "__main__":
	main(sys.argv)
