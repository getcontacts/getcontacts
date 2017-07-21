# Author: Anthony Kai Kwang Ma
# Email: akma327@stanford.edu

import sys
import os
import glob
import re

USAGE_STR= """
# This script stitches individual pieces of dynamic interactions computed from small trajectories
# and then returns a file containing all dynamic interactions across the frames for the larger 
# combined trajectory. 

# Usage: 
# python stitchFrames.py <EXP_REP_PATH> <FRAGMENT_ID> <TITRATION_ID> <FILE_NAME_ID> <OUTPUT_PATH> <OUTPUT_FILE_NAME>

# Arguments:
# <EXP_REP_PATH> Absolute path to the folder containing the dynamic interaction output for small pieces of trajectories
# <FRAGMENT_ID> Identifier for prefix in the name of folders storing fragments of trajectory and dynamic interaction output
# <TITRATION_ID> Name of folder that stores output for specific geometric criterion of certain interaction types 
# <FILE_NAME_ID> Name of files that we want to combine across all the fragments 
# <OUTPUT_PATH> Absolute path to output file containing the stitched dynamic interactions 
# <OUTPUT_FILE_NAME> Name of the file containing stitched dynamic interactions 


# Example:
EXP_REP_PATH="/scratch/PI/rondror/akma327/noncovalent_Interaction_Scripts/mOR-InteractionOutput/MOR_active_waters/mor_active_refine150_agonist_noNanobody/rep_1"
FRAGMENT_ID="Prod"
TITRATION_ID="HB_3.5A_70D"
FILE_NAME_ID="water_bond_result.txt"
OUTPUT_PATH="/scratch/PI/rondror/akma327/noncovalent_Interaction_Scripts/mOR-InteractionOutput/MOR_active_waters/mor_active_refine150_agonist_noNanobody/rep_1"
OUTPUT_FILE_NAME="water_bond_result_stitch.txt"
python stitchFrames.py $EXP_REP_PATH $FRAGMENT_ID $TITRATION_ID $FILE_NAME_ID $OUTPUT_PATH $OUTPUT_FILE_NAME

"""

K_MIN_ARG = 7

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split('(\d+)', text) ]


if __name__ == "__main__":
	if(len(sys.argv) < K_MIN_ARG):
		print(USAGE_STR)
		exit(1)
	exp_rep_path = sys.argv[1]
	folder_id = sys.argv[2]
	titration_id = sys.argv[3]
	file_name_id = sys.argv[4]
	output_path = sys.argv[5]
	output_file_name = sys.argv[6]
	if not os.path.exists(output_path):
		os.makedirs(output_path)
	f = open(output_path + "/" + output_file_name, 'w')
	if(titration_id == "None"):
		prod_files = glob.glob(exp_rep_path + "/" + folder_id + "*/" + file_name_id)
	else:
		prod_files = glob.glob(exp_rep_path + "/" + folder_id + "*/" + titration_id + "/" + file_name_id)
	prod_files.sort(key=natural_keys)
	frameIndex = -1
	totalComputeTime = 0
	for index, piece_file in enumerate(prod_files):
		p = open(piece_file, 'r')
		for line in p:
			if("Stride:" in line):
				if(index == 0): f.write(line)
				else: continue
			elif("TrajectoryPath:" in line):
				f.write(line)
			elif("TopologyPath:" in line):
				if(index == 0): f.write(line)
				else: continue
			elif("nFrames:" in line):
				continue
			elif("Salt Bridges:" in line):
				if(index == 0): f.write(line)
				else: continue
			elif("Pi-Cation:" in line):
				if(index == 0): f.write(line)
				else: continue
			elif("Pi-Stacking:" in line):
				if(index == 0): f.write(line)
				else: continue
			elif("T-Stacking:" in line):
				if(index == 0): f.write(line)
				else: continue
			elif("Van Der Waals:" in line):
				if(index == 0): f.write(line)
				else: continue
			elif("Hydrogen Bonds:" in line):
				if(index == 0): f.write(line)
				else: continue	
			elif("Hydrogen Bond-Water Mediated:" in line):
				if(index == 0): f.write(line)
				else:continue
			elif("Backbone Backbone Hydrogen Bonds:" in line):
				if(index == 0): f.write(line)
				else: continue	
			elif("Sidechain Backbond Hydrogen Bonds:" in line):
				if(index == 0): f.write(line)
				else: continue	
			elif("Sidechain Sidechain Hydrogen Bonds:" in line):
				if(index == 0): f.write(line)
				else: continue	
			elif("Residue Water Hydrogen Bonds:" in line):
				if(index == 0): f.write(line)
				else:continue
			elif("Water Bonds:" in line):
				if(index == 0): f.write(line)
				else: continue
			elif("Extended Water Bonds:" in line):
				if(index == 0): f.write(line)
				else:continue
			elif("Ligand Backbone Hydrogen Bonds:" in line):
				if(index == 0): f.write(line)
				else:continue
			elif("Ligand Sidechain Hydrogen Bonds:" in line):
				if(index == 0): f.write(line)
				else:continue
			elif("Ligand Water Hydrogen Bonds:" in line):
				if(index == 0): f.write(line)
				else:continue
			elif("Ligand Water Bonds:" in line):
				if(index == 0): f.write(line)
				else:continue
			elif("Ligand Extended Water Bonds:" in line):
				if(index == 0): f.write(line)
				else:continue
			elif("Salt Bridge Frame: " in line): 
				frameIndex += 1
				f.write("Salt Bridge Frame: " + str(frameIndex) + "\n")
			elif("Pi_Cation Frame: " in line):
				frameIndex += 1
				f.write("Pi_Cation Frame: " + str(frameIndex) + "\n")
			elif("Pi_Stacking Frame: " in line):
				frameIndex += 1
				f.write("Pi_Stacking Frame: " + str(frameIndex) + "\n")
			elif("T_Stacking Frame: " in line):
				frameIndex += 1
				f.write("T_Stacking Frame: " + str(frameIndex) + "\n")
			elif("Vanderwaal Frame: " in line):
				frameIndex += 1
				f.write("Vanderwaal Frame: " + str(frameIndex) + "\n")
			elif("Hydrogen Bond Water Frame: " in line):
				frameIndex += 1
				f.write("Hydrogen Bond Water Frame: " + str(frameIndex) + "\n")
			elif("Backbone Backbone Hydrogen Bond Frame:" in line):
				frameIndex += 1
				f.write("Backbone Backbone Hydrogen Bond Frame:" + str(frameIndex) + "\n")
			elif("Sidechain Backbond Hydrogen Bond Frame:" in line):
				frameIndex += 1
				f.write("Sidechain Backbond Hydrogen Bond Frame:" + str(frameIndex) + "\n")
			elif("Sidechain Sidechain Hydrogen Bond Frame:" in line):
				frameIndex += 1
				f.write("Sidechain Sidechain Hydrogen Bond Frame:" + str(frameIndex) + "\n")
			elif("Residue Water Hydrogen Bond Frame:" in line):
				frameIndex += 1
				f.write("Residue Water Hydrogen Bond Frame:" + str(frameIndex) + "\n")
			elif("Water Bond Frame:" in line):
				frameIndex += 1
				f.write("Water Bond Frame:" + str(frameIndex) + "\n")
			elif("Extended Water Bond Frame:" in line):
				frameIndex += 1
				f.write("Extended Water Bond Frame:" + str(frameIndex) + "\n")
			elif("Ligand Backbone Hydrogen Bond Frame:" in line):
				frameIndex += 1
				f.write("Ligand Backbone Hydrogen Bond Frame:" + str(frameIndex) + "\n")
			elif("Ligand Backbone Hydrogen Bond Frame:" in line):
				frameIndex += 1
				f.write("Ligand Backbone Hydrogen Bond Frame:" + str(frameIndex) + "\n")
			elif("Ligand Sidechain Hydrogen Bond Frame:" in line):
				frameIndex += 1
				f.write("Ligand Sidechain Hydrogen Bond Frame:" + str(frameIndex) + "\n")
			elif("Ligand Water Hydrogen Bond Frame:" in line):
				frameIndex += 1
				f.write("Ligand Water Hydrogen Bond Frame:" + str(frameIndex) + "\n")
			elif("Ligand Water Bond Frame:" in line):
				frameIndex += 1
				f.write("Ligand Water Bond Frame:" + str(frameIndex) + "\n")
			elif("Ligand Extended Water Bond Frame:" in line):
				frameIndex += 1
				f.write("Ligand Extended Water Bond Frame:" + str(frameIndex) + "\n")
			elif("Computing Time:" in line):
				totalComputeTime += float(line.split("Computing Time:")[1])
			else:
				f.write(line)
	f.write("nFrames:" + str(frameIndex) + "\n")
	f.write("Computing Time:" + str(totalComputeTime) + "\n")

