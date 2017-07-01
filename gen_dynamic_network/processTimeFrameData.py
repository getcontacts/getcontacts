# Author: Anthony Kai Kwang Ma
# Email: akma327@stanford.edu
# Date: February 18, 2016

# Process Noncovalent Interaction by Time Frame for MOR Trajectories- Anthony Ma 02/18/16

from __future__ import print_function
import os
import sys
import datetime
import time
import itertools
import copy
import mdtraj as md
from hbond import baker_hubbard
from saltbridge import *
from vanderwaal import *
from hydrophobe import *
from pication import *
from aromatic import *
from waterMediatedHBonds import * 

# Usage:
# python processTimeFrameData.py <Stitch Interaction Input Path> <Input File Name> <Binary Dictionary Output Path> <Output File Name> <interaction type><binarize flag>
# Example:
# StitchedInteractionPath="/scratch/PI/rondror/akma327/noncovalent_Interaction_Scripts/mOR-InteractionOutput/MOR_active_waters/mor_active_refine150_agonist_noNanobody_noNTerm/rep_2"
# INPUT_FILE_PATH="salt_bridge_result_stitch.txt"
# OUTPUT_PATH="/scratch/PI/rondror/akma327/noncovalent_Interaction_Scripts/mOR-InteractionOutputDictionary/MOR_active_waters/mor_active_refine150_agonist_noNanobody_noNTerm/rep_2"
# OUTPUT_FILE_NAME="salt_bridge_result_dict.txt"
# python processTimeFrameData.py $StitchedInteractionPath $INPUT_FILE_PATH $OUTPUT_PATH $OUTPUT_FILE_NAME -sb -b


FLAG_INTERACTION = {"-sb": "Salt Bridges", "-pc": "Pi-Cation", 
"-ps": "Pi-Stacking", "-ts": "T_Stacking", "-vdw": "Van Der Waals", 
"-hb": "Hydrogen Bonds", "-hbw": "Hydrogen Bond-Water Mediated",
"-hbbb": "Backbone Backbone Hydrogen Bonds",
"-hbsb": "Sidechain Backbond Hydrogen Bonds",
"-hbss": "Sidechain Sidechain Hydrogen Bonds",
"-rw": "Residue Water Hydrogen Bonds",
"-wb": "Water Bonds",
"-wb2": "Extended Water Bonds",
"-hlb": "Ligand Backbone Hydrogen Bonds",
"-hls": "Ligand Sidechain Hydrogen Bonds",
"-lw": "Ligand Water Hydrogen Bonds",
"-lwb": "Ligand Water Bonds",
"-lwb2": "Ligand Extended Water Bonds"}

def retrieveDataFileReader():
	print("retrieveDataFileReader()")
	input_path = sys.argv[1]
	input_filename = sys.argv[2]
	filename = input_path + "/" + input_filename
	print("Creating Binary Dictionary for: " + filename)
	if os.path.isfile(filename):
		f = open(filename, 'r')
		return f
	else:
		print("Interaction Results Not Found")
		return None

def createTimeFrameDict():
	numLinesRead = 0
	print("createTimeFrameDict()")
	f = retrieveDataFileReader()
	stride = 1
	TrajectoryPaths = []
	TopologyPath = ""
	maxDelta = 0
	appendcount = 0
	#dictionary to store pair of interaction key with list of time points
	timelapse_interaction_by_chain_dict = {}
	if(f != None):
		currFrame = 0
		currProtein = ""
		for line in f:
			# if(numLinesRead % 1000 == 0): 
			# 	print("Line: " + str(numLinesRead) + " ** " + str(line))
			# 	if(0 in timelapse_interaction_by_chain_dict):
			# 		print("len dict " + str(len(timelapse_interaction_by_chain_dict[0])))
			numLinesRead += 1
			if("Stride:" in line):
				stride = int(line.split("Stride:")[1])
			elif("TrajectoryPath:" in line):
				TrajectoryPaths.append(line.split("TrajectoryPath:")[1])
			elif("TopologyPath:" in line):
				TopologyPath = line.split("TopologyPath:")[1]
			elif("nFrames:" in line):
				continue
			elif("Computing Time:" in line):
				continue
			elif("Frame: " in line):
				currFrame = int(line.split("Frame: ")[1])
			elif("Salt Bridges:" in line):
				currProtein = (line.split("Salt Bridges:")[1])
			elif("Pi-Cation:" in line):
				currProtein = (line.split("Pi-Cation:")[1])
			elif("Pi-Stacking:" in line):
				currProtein = (line.split("Pi-Stacking:")[1])
			elif("T-Stacking:" in line):
				currProtein = (line.split("T-Stacking:")[1])
			elif("Van Der Waals:" in line):
				currProtein = (line.split("Van Der Waals:")[1])
			elif("Hydrogen Bonds:" in line):
				currProtein = (line.split("Hydrogen Bonds:")[1])
			elif("Hydrogen Bond-Water Mediated:" in line):
				currProtein = (line.split("Hydrogen Bond-Water Mediated:")[1])
			elif("Backbone Backbone Hydrogen Bonds:" in line):
				currProtein = (line.split("Backbone Backbone Hydrogen Bonds:")[1])
			elif("Sidechain Backbond Hydrogen Bonds:" in line):
				currProtein = (line.split("Sidechain Backbond Hydrogen Bonds:")[1])
			elif("Sidechain Sidechain Hydrogen Bonds:" in line):
				currProtein = (line.split("Sidechain Sidechain Hydrogen Bonds:")[1])
			elif("Residue Water Hydrogen Bonds:" in line):
				currProtein = (line.split("Residue Water Hydrogen Bonds:")[1])
			elif("Water Bonds:" in line):
				currProtein = (line.split("Water Bonds:")[1])
			elif("Extended Water Bonds:" in line):
				currProtein = (line.split("Extended Water Bonds:")[1])
			elif("Ligand Backbone Hydrogen Bonds:" in line):
				currProtein = (line.split("Ligand Backbone Hydrogen Bonds:")[1])
			elif("Ligand Sidechain Hydrogen Bonds:" in line):
				currProtein = (line.split("Ligand Sidechain Hydrogen Bonds:")[1])
			elif("Ligand Water Hydrogen Bonds:" in line):
				currProtein = (line.split("Ligand Water Hydrogen Bonds:")[1])
			elif("Ligand Water Bonds:" in line):
				currProtein = (line.split("Ligand Water Bonds:")[1])
			elif("Ligand Extended Water Bonds:" in line):
				currProtein = (line.split("Ligand Extended Water Bonds:")[1])
			
			

			elif(line != "\n"):
				# if(numLinesRead >100000):
				# 	break
				key_pair = tuple(line.strip().split(" -- "))
				split_key = key_pair[0].split("_")
				if(len(split_key) > 1):
					chain_id = int(split_key[1])
				else:
					chain_id = 0
				if(chain_id not in timelapse_interaction_by_chain_dict.keys()):
					timelapse_interaction_by_chain_dict[chain_id] = {} #create empty dictionary
					timelapse_interaction_by_chain_dict[chain_id][key_pair] = [currFrame]
				else:
					tic1 = datetime.datetime.now()
					if key_pair not in timelapse_interaction_by_chain_dict[chain_id].keys():
						toc1 = datetime.datetime.now()
						deltaLookup = toc1 - tic1 
						deltaLookup =deltaLookup.microseconds
						appendcount += 1
						tic = datetime.datetime.now()
						timelapse_interaction_by_chain_dict[chain_id][key_pair] = [currFrame]
						toc = datetime.datetime.now()
						delta = toc - tic
						delta = delta.microseconds
						# if (delta > maxDelta): maxDelta = delta
						# if (maxDelta == delta or appendcount %50 == 0): 
						if(numLinesRead %1000 ==0):
							print("Line Number: " + str(numLinesRead) + " Lookup Time1: " + str(deltaLookup))
							print("time for append: " + str(delta) + " Dict Size: " + str(len(timelapse_interaction_by_chain_dict[chain_id])))

					else:
						toc2 = datetime.datetime.now()
						timelapse_interaction_by_chain_dict[chain_id][key_pair].append(currFrame)
						deltaLookup = toc2 -tic1
						deltaLookup =deltaLookup.microseconds
						if(numLinesRead %1000 == 0):
							print("Line Number: " + str(numLinesRead) + " Lookup Time2: " + str(deltaLookup) + " Dict Size: " + str(len(timelapse_interaction_by_chain_dict[chain_id])))
	totalFrames = currFrame + 1
	return stride, TrajectoryPaths, TopologyPath, timelapse_interaction_by_chain_dict, totalFrames

def displayFrequencyMap(stride, TrajectoryPath, TopologyPath, timelapse_interaction_by_chain_dict):
	print("Stride:" + str(stride))
	print("TrajectoryPath:" + TrajectoryPath)
	print("TopologyPath:" + TopologyPath)
	for chain_id in timelapse_interaction_by_chain_dict.keys():
		print("Dictionary for Chain: " + str(chain_id))
		for interaction in timelapse_interaction_by_chain_dict[chain_id].keys():
			key_identifier = interaction[0] + " -- " + interaction[1]
			print(key_identifier + ":\n" + str(timelapse_interaction_by_chain_dict[chain_id][interaction]))

def writeFrequencyMap(stride, TrajectoryPaths, TopologyPath, timelapse_interaction_by_chain_dict, totalFrames):
	print("writeFrequencyMap()")
	output_path = sys.argv[3]
	output_filename = sys.argv[4]
	interaction_selection = sys.argv[5]
	filename = output_path + "/" + output_filename
	if not os.path.exists(os.path.dirname(filename)):
		os.makedirs(os.path.dirname(filename))
	f = open(filename, 'w')
	f.write("Stride:" + str(stride) + "\n")
	for TrajectoryPath in TrajectoryPaths:
		f.write("TrajectoryPath:" + TrajectoryPath + "\n")
	f.write("TopologyPath:" + TopologyPath + "\n")
	f.write(FLAG_INTERACTION[interaction_selection] + " Dictionary Heat Map " + "\n")
	f.write("TotalFrames: " + str(totalFrames) + "\n")
	for chain_id in timelapse_interaction_by_chain_dict.keys():
		f.write("Dictionary for Chain: " + str(chain_id) + "\n")
		for interaction in timelapse_interaction_by_chain_dict[chain_id].keys():
			key_identifier = interaction[0] + " -- " + interaction[1]
			f.write(key_identifier + "~" + str(timelapse_interaction_by_chain_dict[chain_id][interaction]) + "\n")
	print("Wrote results to: " + filename)

# converts a list of time points in which interaction is present to a binary
# representation over all time points. eg. [1,2,5] --> [0,1,1,0,0,1]
def binarizeTimePoints(timelapse_interaction_by_chain_dict, totalFrames):
	print("binarizeTimePoints()")
	def binarizeTime(timePoints, totalFrames):
		binaryVec = [0]*(totalFrames)
		for tp in timePoints:
			binaryVec[tp] = 1
		return binaryVec
	binary_timelapse_interaction_by_chain_dict = copy.deepcopy(timelapse_interaction_by_chain_dict)
	for chain_id in binary_timelapse_interaction_by_chain_dict.keys():
		for interaction in binary_timelapse_interaction_by_chain_dict[chain_id].keys():
			timePoints = binary_timelapse_interaction_by_chain_dict[chain_id][interaction]
			binarizedTime = binarizeTime(timePoints, totalFrames)
			binary_timelapse_interaction_by_chain_dict[chain_id][interaction] = binarizedTime
	return binary_timelapse_interaction_by_chain_dict


# Create the frequency map either in time point form, or binaryized form with -b flag
def calculateFrequencyMap():
	print("calculateFrequencyMap()")
	binarizeFlagIndex = 6
	stride, TrajectoryPaths, TopologyPath, timelapse_interaction_by_chain_dict, totalFrames = createTimeFrameDict()
	if(len(sys.argv) > binarizeFlagIndex and sys.argv[binarizeFlagIndex] == "-b"):
		binary_timelapse_interaction_by_chain_dict = binarizeTimePoints(timelapse_interaction_by_chain_dict, totalFrames)
		writeFrequencyMap(stride, TrajectoryPaths, TopologyPath, binary_timelapse_interaction_by_chain_dict, totalFrames)
	else:
		writeFrequencyMap(stride, TrajectoryPaths, TopologyPath, timelapse_interaction_by_chain_dict, totalFrames)

if __name__ == "__main__":
	tic = time.clock()
	print ("START MOR_processTimeFrameData() TIME: " + datetime.datetime.now().strftime("%A, %d. %B %Y %I:%M%p"))
	calculateFrequencyMap()
	toc = time.clock()
	print ("END MOR_processTimeFrameData() TIME: " + datetime.datetime.now().strftime("%A, %d. %B %Y %I:%M%p"))
	print("Total Time for MOR_processTimeFrameData: " + str(toc-tic))


