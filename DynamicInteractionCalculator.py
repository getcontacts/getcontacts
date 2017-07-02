# Author: Anthony Kai Kwang Ma
# Email: akma327@stanford.edu
# Date: October 1, 2015

# DynamicInteractionCalculator.py

from __future__ import print_function
import os
import sys
import time
import itertools
import mdtraj as md
from contact_calc.saltbridge import *
from contact_calc.pication import *
from contact_calc.aromatic import *
from contact_calc.vanderwaal import *
from contact_calc.vmdMeasureHBonds import * 
from contact_calc.stratifyHBonds import *
from contact_calc.stratifyLigandHBonds import *

USAGE_STR= """
# Usage:
# python DynamicInteractionCalculator.py <TRAJ_PATH> <TOP_PATH> <OUTPUT_DIRECTORY_PATH> <OUTPUT_FILE_IDENTIFIER> <Interaction Type> <optional -process > <optional file path for hbond to post process> <optional -stride flag > <stride value if not 1> <optional -solv if not water> <solvent> <optional -chain> <chainID> <optional -ligand> <Ligand>

# Arguments:
# <TRAJ_PATH> Absolute Path to trajectory 
# <TOP_PATH> Absolute path to topology 
# <OUTPUT_DIRECTIORY_PATH> Absolute path to folder to store the dynamic interaction output 
# <OUTPUT_FILE_IDENTIFIER> Name of the file inserted into the OUTPUT_DIRECTORY_PATH 
# <Interaction Type> Flag to denote the type of interaction being calculated (ie -sb, -pc, -ps, -ts, etc)
# <optional -process flag> To denote whether extra file path for the raw water mediated hydrogen bond file needs to be processed
# <optional -stride flag> To denote a stride value other than default of 1
# <optional -solv flag> To denote a solvent id other than default of TIP3
# <optional -chain flag> To denote the specific chain ID to query when using VMD's hydrogen bond calculator 
# <optional -ligand flag> To denote the resname of ligand in the simulation.

# Example 
TRAJ_PATH="Prod_0_reimaged.nc"
TOP_PATH="step5_assembly.pdb"
OUTPUT_DIRECTORY_PATH="output_contacts"
OUTPUT_FILE_IDENTIFIER="hbsb.txt"
python DynamicInteractionCalculator.py $TRAJ_PATH $TOP_PATH $OUTPUT_DIRECTORY_PATH $OUTPUT_FILE_IDENTIFIER -hbsb

"""

K_MIN_ARG = 5

def createFileWriter(OUTPUT):
	"""
		Create file descripter to write out non-covalent interactions for a simulation fragment
	"""
	if not os.path.exists(os.path.dirname(OUTPUT)):
		os.makedirs(os.path.dirname(OUTPUT))
	f = open(OUTPUT, 'w')
	return f


def calcDynamicInteractions(TOP, TRAJ, OUTPUT, INTERACTION_TYPE, post_process_file=None, stride=1, solventId=None, solvent=None, chainId=None, ligand=None):
	"""
		Calculate non-covalent interaction of specified type for a particular simulation trajectory fragment. 
	"""
	print("Begin calculating non-covalent interaction ...")
	f = createFileWriter(OUTPUT)
	f.write("Stride:" + str(stride) + "\n")
	f.write("TrajectoryPath:" + TRAJ + "\n")
	f.write("TopologyPath:" + TOP + "\n")


	### MD Traj Interactions
	if(INTERACTION_TYPE in ['-sb', '-pc', '-ps', '-ts', '-vdw']):
		traj = md.load(TRAJ_PATH, top = TOP_PATH)
		traj = traj[::stride] #take stride value into account 
	if(INTERACTION_TYPE == '-sb'):
		calcSaltBridgeResults(traj, file_writer, TOP_PATH)
	if(INTERACTION_TYPE == '-pc'):
		calcPiCationResults(traj, file_writer, TOP_PATH)
	if(INTERACTION_TYPE == '-ps'):
		calcPiStackingResults(traj, file_writer, TOP_PATH)
	if(INTERACTION_TYPE == '-ts'):
		calcTStackingResults(traj, file_writer, TOP_PATH)
	if(INTERACTION_TYPE == '-vdw'):
		calcVanderwaalsResults(traj, file_writer, TOP_PATH)


	### VMD HBond Interactions
	if(INTERACTION_TYPE == '-hbw'):
		calcHydrogenBondWaterResults(file_writer, stride, TOP_PATH, TRAJ_PATH, solventId, chainId)
	if(INTERACTION_TYPE == '-hbbb'):
		calcBackBoneHBondResults(file_writer, post_process_file)
	if(INTERACTION_TYPE == '-hbsb'):
		calcSideChainBackBoneHBondResults(file_writer, post_process_file)
	if(INTERACTION_TYPE == '-hbss'):
		calcSideChainHBondResults(file_writer, post_process_file)
	if(INTERACTION_TYPE == '-rw'):
		calcResidueWaterBondResults(file_writer, post_process_file)
	if(INTERACTION_TYPE == '-wb'):
		calcWaterBondResults(file_writer, post_process_file)
	if(INTERACTION_TYPE == '-wb2'):
		calcWaterBond2Results(file_writer, post_process_file)

	### MVD Ligand HBond Interactions
	if(INTERACTION_TYPE == "-lhbw"):
		calcLigandHydrogenBondWaterResults(file_writer, stride, TOP_PATH, TRAJ_PATH, solventId, chainId, ligand)
	if(INTERACTION_TYPE == "-hlb"):
		calcLigandBackBoneHBondResults(file_writer, post_process_file, ligand)
	if(INTERACTION_TYPE == "-hls"):
		calcLigandSideChainHBondResults(file_writer, post_process_file, ligand)
	if(INTERACTION_TYPE == '-lw'):
		calcLigandToWaterBondResults(file_writer, post_process_file, ligand)
	if(INTERACTION_TYPE == "-lwb"):
		calcLigandWaterBondResults(file_writer, post_process_file, ligand)
	if(INTERACTION_TYPE == "-lwb2"):
		calcLigandWaterBond2Results(file_writer, post_process_file, ligand)


if __name__ == "__main__":
	if(len(sys.argv) < K_MIN_ARG):
		print(USAGE_STR)
		exit(1)

	### Required arguments
	(TOP, TRAJ, OUTPUT, INTERACTION_TYPE) = (sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

	### Optional arguments
	post_process_file = None
	stride = 1
	solvent = None
	solventId = None 
	chainId = None
	ligand = None

	if("-process" in sys.argv):
		process_index = sys.argv.index("-process")
		post_process_file = sys.argv[process_index + 1]

	if("-stride" in sys.argv):
		stride_index = sys.argv.index("-stride")
		stride = int(sys.argv[stride_index + 1])

	if("-solv" in sys.argv):
		solv_index = sys.argv.index("-solv")
		solventId = sys.argv[solvIndex + 1]
		solvent = "resname" + solventId

	if("-chain" in sys.argv):
		chain_index = sys.argv.index("-chain")
		chainId = str(sys.argv[chain_index + 1])

	if("-ligand" in sys.argv):
		ligand_index = sys.argv.index("-ligand")
		ligand = str(sys.argv[ligand_index + 1])

	calcDynamicInteractions(TOP, TRAJ, OUTPUT, INTERACTION_TYPE, post_process_file, stride, solventId, solvent, chainId, ligand)
