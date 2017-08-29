from collections import defaultdict
import sys
import os
import errno
import re

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


def compile_interactions(compiled_path, output_filename):
	interactions = defaultdict(dict)
	respairs = set()
	for inter_file in listdir_files(compiled_path):
		file_descriptor = get_file_descriptor(inter_file)
		interaction_type, _ = file_descriptor.split('_')

		with open(inter_file, 'r') as open_inter_file:
			inter_file_lines = [line for line in open_inter_file.readlines()]

		receptor_states = inter_file_lines[0].strip().split()[1:]

		for inter_file_line in inter_file_lines[1:]:
			split_inter_file_line = inter_file_line.strip().split()
			respair = split_inter_file_line[0] + ":" + interaction_type
			respairs.add(respair)
			respair_frequencies = split_inter_file_line[1:]
			
			for i in xrange(len(receptor_states)):
				interactions[receptor_states[i]][respair] = respair_frequencies[i]

	with open(output_filename, 'w+') as output_file_open:
		output_file_open.write("ResPair")
		for receptor_state in sorted(interactions.keys()):
			output_file_open.write('\t%s' % receptor_state)
		output_file_open.write('\n')
		for respair in respairs:
			output_file_open.write(respair)
			for receptor_state in sorted(interactions.keys()):
				if respair in interactions[receptor_state]:
					output_file_open.write('\t%s' % interactions[receptor_state][respair])
				else:
					output_file_open.write('\t%.8f' % 0)
			output_file_open.write('\n')
		
def main(argv):
	compiled_path = clean_path(argv[1])
	output_filename = '%scompiled_ints.txt' % (compiled_path)
	compile_interactions(compiled_path, output_filename)

if __name__ == "__main__":
	main(sys.argv)
