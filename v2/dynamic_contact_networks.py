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
from contact_calc.compute_contacts import *

USAGE_STR = """

# Purpose 
# Computes non-covalent contacts in MD simulation for any protein of study

# Usage
# python DynamicContactNetworks.py <TOP> <TRAJ> <OUTPUT_DIR> <stride> <solv> <chain> <ligand> <INTERACTION_TYPES>

# Arguments
# <TOP> Absolute path to topology
# <TRAJ> Absolute path to reimaged MD trajectory
# <OUTPUT_DIR> Absolute path to output directory to store contacts data
# <INTERACTION_TYPES> -itypes flag followed by a list of non-covalent interaction types to compute (ie -sb -pc -ps -ts -vdw, etc)
# <optional -stride flag> To denote a stride value other than default of 1
# <optional -solv flag> To denote a solvent id other than default of TIP3
# <optional -chain flag> To denote the specific chain ID to query when using VMD's hydrogen bond calculator 
# <optional -ligand flag> To denote the resname of ligand in the simulation.

# Example
TOP=topology.pdb"
TRAJ="trajectory.dcd"
OUTPUT_DIR="output"
python dynamic_contact_networks.py $TOP $TRAJ $OUTPUT_DIR -itype -hb -hlb

"""

K_MIN_ARG = 4

if __name__ == "__main__":
	if(len(sys.argv) < K_MIN_ARG):
		print(USAGE_STR)
		exit(1)

	### Required Arguments
	(TOP, TRAJ, OUTPUT_DIR) = (sys.argv[1], sys.argv[2], sys.argv[3])
	ITYPES = sys.argv[sys.argv.index('-itype') + 1:]

	### Optional Arguments
	stride = 1 # default
	solvent_resn = "TIP3" # default
	chain_id = None
	ligand = None

	if("-stride" in sys.argv):
		stride_index = sys.argv.index("-stride")
		stride = int(sys.argv[stride_index + 1])

	if("-solv" in sys.argv):
		solv_index = sys.argv.index("-solv")
		solvent_resn = sys.argv[solvIndex + 1]

	if("-chain" in sys.argv):
		chain_index = sys.argv.index("-chain")
		chain_id = str(sys.argv[chain_index + 1])

	if("-ligand" in sys.argv):
		ligand_index = sys.argv.index("-ligand")
		ligand = str(sys.argv[ligand_index + 1])

	tic = datetime.datetime.now()
	compute_dynamic_contacts(TOP, TRAJ, OUTPUT_DIR, ITYPES, stride, solvent_resn, chain_id, ligand)
	toc = datetime.datetime.now()
	print("Computation Time: " + str((toc-tic).total_seconds()))


