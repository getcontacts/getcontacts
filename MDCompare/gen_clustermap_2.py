from __future__ import division
import sys
import json
import re
import pandas as pd
import seaborn as sns
import os
from utils import *
import matplotlib.pyplot as plt

def produce_tanimoto(filename):
	if re.search("vdw", filename) or re.search("hbbb", filename) or re.search("wb2", filename):
		return
	sns.set(color_codes = True)
	print "Producing clustermap for %s:" % (filename)
	open_table = pd.read_csv(filename, header = 0, index_col = 0)
	if open_table.empty:
		print "Warning: %s is empty, so a clustermap will not be produced for this file" % filename
		return
	try:
		cluster = sns.clustermap(open_table, row_cluster = True, col_cluster = True, figsize=(open_table.shape[1], open_table.shape[0]), annot=True)
	except ValueError:
		return
	plt.setp(cluster.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
	plt.setp(cluster.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
	plt.title(os.path.basename(filename))

	cluster.savefig('%s_heatmap.png' % os.path.splitext(filename)[0])


def main(argv):
	paths = []

	with open("./sample_inputs/sample_input.csv", 'r') as ropen:
		matrix_lines = [line.strip() for line in ropen.readlines()]
	for line in matrix_lines:
		if len(line) == 0:
			continue
		descriptor, path, _ = line.split(',')
		if re.search(argv[1], descriptor):
			paths.append(clean_path(path))

	for path in paths:
		path = path.replace("/scratch/PI/rondror/augustine", "../..")
		for filename in listdir_files(path):
			if re.search("_byatom_tanimoto\.csv", filename):
				produce_tanimoto(filename)
			if re.search("_byres_generic_tanimoto\.csv", filename):
				produce_tanimoto(filename)

if __name__ == "__main__":
	main(sys.argv)
