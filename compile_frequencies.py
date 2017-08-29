from collections import defaultdict
import sys
import os
import errno
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


def compile_freqs(dirname, output_filename):
	int_to_file_to_freq = defaultdict(defaultdict)

	for freq_file in listdir_files(dirname):
		short_filename = os.path.splitext(os.path.basename(freq_file))[0]

		print short_filename

		with open(freq_file, 'r') as open_freq_file:
			freq_file_lines = [line.strip().split() for line in open_freq_file.readlines()]
		
		for interaction in freq_file_lines:
			res1 = interaction[0]
			res2 = interaction[1]
			frequency = interaction[2]
			interaction_pair = res1 + "-" + res2

			int_to_file_to_freq[interaction_pair][short_filename] = frequency

	with open(output_filename, 'w+') as output_file_open:
		output_file_open.write("ResPair")
		for freq_file in sorted(listdir_files(dirname)):
			short_filename = os.path.splitext(os.path.basename(freq_file))[0]
			output_file_open.write("\t" + short_filename)
		output_file_open.write("\n")

		for interaction_pair in int_to_file_to_freq:
			writeLine = False
			output_line = interaction_pair
#			output_file_open.write(interaction_pair)

			for freq_file in sorted(listdir_files(dirname)):
				short_filename = os.path.splitext(os.path.basename(freq_file))[0]
				if short_filename in int_to_file_to_freq[interaction_pair]:
					frequency = float(int_to_file_to_freq[interaction_pair][short_filename])
					if frequency >= 0.5:
						writeLine = True
					output_line += "\t%s" % (frequency)
				else:
					output_line += "\t%.8f" % (0.0)
#				output_file_open.write(output_string)
			
			output_line += "\n"
#			output_file_open.write("\n")
			if writeLine:
				output_file_open.write(output_line)

def main(argv):
	smashed_path = clean_path(argv[1])
	outputs_path = clean_path(argv[2])
	interaction_types = [int_type[1:] for int_type in argv[3:]]
	open_dir(outputs_path)

	states = ['active', 'inactive']

	for interaction in interaction_types:
		for state in states:
			dirname = '%s%s/%s/' % (smashed_path, interaction, state)
			output_filename = '%s%s_%s.txt' % (outputs_path, interaction, state)
			print dirname
			compile_freqs(dirname, output_filename)
			print ""

if __name__ == "__main__":
	main(sys.argv)
