##############################################################################
# MDContactNetworks: A Python Library for computing non-covalent contacts
#                    throughout Molecular Dynamics Trajectories. 
# Copyright 2016-2017 Stanford University and the Authors
#
# Authors: Anthony Kai Kwang Ma
# Email: anthony.ma@yale.edu, anthonyma27@gmail.com, akma327@stanford.edu
##############################################################################


##############################################################################
# Imports
##############################################################################

import os
import sys
import datetime
import argparse
from contact_calc.compute_contacts import *

HELP_STR = """
 ===============================================
|                 Anthony Ma, 2017              |
|               Stanford University             |
|                  Version 1.0.0                |
 ===============================================

Command was:
python dynamic_contact_networks.py --help

usage: python dynamic_contact_networks.py [--help] [--topology TOPOLOGY] 
					  [--trajectory TRAJECTORY]
					  [--output_dir OUTPUT_DIRECTORY] 
					  [--itype INTERACTION_TYPES]
					  [--cores NUM_CORES]
					  [--solv SOLVENT]
					  [--sele SELECTION]
					  [--ligand LIGAND]

optional arguments:
	--help				show this help message and exit
	--topology TOPOLOGY		path to topology file 
	--trajectory TRAJECTORY		path to trajectory file
	--output_dir OUTPUT_DIRECTORY	path to output directory
	--itype INTERACTION_TYPES	list of interaction type flags
	--cores NUM_CORES		number of cpu cores to parallelize upon
	--solv SOLVENT 			resname of solvent molecule
	--sele SELECTION 		atom selection query in VMD
	--ligand LIGAND 		resname of ligand molecule 

interaction type flags:
	-sb				salt bridges 
	-pc				pi-cation 
	-ps				pi-stacking
	-ts				t-stacking
	-vdw				vanderwaals
	-hb				hydrogen bonds
	-lhb				ligand hydrogen bonds
	-hbbb				backbone-backbone hydrogen bonds
	-hbsb				backbone-sidechain hydrogen bonds
	-hbss				sidechain-sidechain hydrogen bonds
	-wb				water bridges
	-wb2				extended water bridges 
	-hls				ligand-sidechain residue hydrogen bonds 
	-hlb				ligand-backbone residue hydrogen bonds 
	-lwb				ligand water bridges
	-lwb2 				extended ligand water bridges 

"""
DESCRIPTION="Computes non-covalent contact networks in MD simulations."


def validate_itypes(ITYPES):
	"""
	Throws error if user specified interaction type is mispelled
	"""
	valid_itypes = ["-sb", "-pc", "-ps", "-ts", "-vdw", "-hb", "-hlb", "-hbbb", 
			"-hbsb", "-hbss", "-wb", "-wb2", "-hls", "-hlb", "-lwb", "-lwb2"]

	for itype in ITYPES:
		if(itype not in valid_itypes):
			print("%s not a valid interaction type ..." % (itype))
			exit(1)

def process_args(args):
	topology = args.topology
	trajectory = args.trajectory
	output_dir = args.output_dir 
	cores = args.cores 
	ligand = args.ligand 
	solv = args.solv 
	sele = args.sele 
	stride = args.stride
	

	return topology, trajectory, output_dir, cores, ligand, solv, sele, stride

def main():
	if("--help" in sys.argv):
		print (HELP_STR)
		exit(1)
	if("--itype" not in sys.argv):
		print("Need to specifiy interaction types ...")
		exit(1)

	### Process arguments
	parser = argparse.ArgumentParser(prog='PROG', add_help=False)
	parser.add_argument('--topology', type=str, default=None, help='path to topology file ')
	parser.add_argument('--trajectory', type=str, default=None, help='path to trajectory file')
	parser.add_argument('--output_dir', type=str, default=None, help='path to output directory')
	parser.add_argument('--cores', type=int, default=6, help='number of cpu cores to parallelize upon')
	parser.add_argument('--solv', type=str, default="TIP3", help='resname of solvent molecule')
	parser.add_argument('--sele', type=str, default=None, help='atom selection query in VMD')
	parser.add_argument('--stride', type=int, default=1, help='skip frames with specified frequency')
	parser.add_argument('--ligand', type=str, default=None, help='resname of ligand molecule ')
	namespace = argparse.Namespace()
	args, unknown = parser.parse_known_args()

	TOP, TRAJ, OUTPUT_DIR, cores, ligand, solv, sele, stride = process_args(args)
	ITYPES = sys.argv[sys.argv.index('--itype') + 1:]
	if("-all" in ITYPES):
		ITYPES = ["-sb", "-pc", "-ps", "-ts", "-vdw", "-hb", "-hlb"]

	validate_itypes(ITYPES)

	### Begin computation
	tic = datetime.datetime.now()
	compute_dynamic_contacts(TOP, TRAJ, OUTPUT_DIR, ITYPES, cores, stride, solv, sele, ligand)
	toc = datetime.datetime.now()
	print("Computation Time: " + str((toc-tic).total_seconds()))

if __name__ == "__main__":
	main()

