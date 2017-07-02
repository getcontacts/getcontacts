# Author: Anthony Kai Kwang Ma
# Email: anthonyma27@gmail.com, akma327@stanford.edu
# Date: 07/02/17
# output_utils.py

from __future__ import print_function
import os
import sys
import time
import itertools
import mdtraj as md


# Open file descriptor for path and returns error message if invalid path 
def readFileOpener(path):
	f = open(path, 'r')
	if(f == None):
		print("File Opener Invalid For: " + str(path))
		exit(0)
	return f

# Write hbbb, hbsb, hbss to file 
def writeToHBondFile(fhb, frameHBondList):
	for elem in frameHBondList:
		fhb.write(elem + "\n")


# Generate combinations of first order water bridges 
def getWBondCombos(water, resid_list):
	wBondCombos = set()
	for pair in itertools.combinations(resid_list, 2):
		resid1, resid2 = pair 
		if(resid1 != resid2):
			wBondCombos.add((resid1[0], resid2[0]))
	return wBondCombos

# Takes the frameWBondDict which contains key = unique solvent identifier, and 
# value = list of all the residues that the unique solvent has h-bond with 
# at the particular frame and writes out all the combinations 
def writeToWBondFile(fwb, frameWBondDict):
	frameWBonds = set() 
	for water in frameWBondDict:
		resid_list = sorted(frameWBondDict[water])
		frameWBonds = frameWBonds | getWBondCombos(water, resid_list)
	frameWBonds = sorted(list(frameWBonds))
	print(frameWBonds)
	for resid1, resid2 in frameWBonds:
		fwb.write(resid1 + " -- " + resid2 + "\n")

# Generate combinations of first order ligand water bridges
def getLWBondCombos(ligand, resid_ligand_list):
	lwBondCombos = set()
	for resid in resid_ligand_list:
		if("LIG" not in resid[0]): # valid residue
			lwBondCombos.add((ligand, resid[0]))
	return lwBondCombos

# Takes the frameLWBond Dict which contains key = unique solvent identifier, and
# value = list of all residues and ligand that the unique solvent has h-bond with
# at the particular frame and writes out all the valid combinations involving
# ligand -- residue
def writeToLWBondFile(fwb, frameLWBondDict, currFrame=-1):
	def find_ligand(arr):
		"""
			Subroutine to find ligand in list of residue + ligand
		"""
		for elem in arr:
			cand_ligand = elem[0]
			if("LIG" in cand_ligand): return cand_ligand
		return None

	frameLWBonds = set()
	for water in frameLWBondDict:
		resid_ligand_list = sorted(frameLWBondDict[water])
		ligand = find_ligand(resid_ligand_list)
		if(ligand == None): continue
		frameLWBonds = frameLWBonds | getLWBondCombos(ligand, resid_ligand_list)
	frameLWBonds = sorted(list(frameLWBonds))
	print(frameLWBonds)
	for ligand, resid in frameLWBonds:
		fwb.write(ligand + " -- " + resid + "\n")


# Takes two lists of residues and returns cross product combination of the residue pairs
def residueCombinations(water1, water2, resid_list1, resid_list2):	
	return set([(res1[0], res2[0]) for res1 in resid_list1 for res2 in resid_list2]) #doesn't make distinct by water 

def residueLigandCombinations(water1, water2, resid_ligand_list1, resid_ligand_list2):
	rlCombos = set()
	for resid_lig1 in resid_ligand_list1:
		for resid_lig2 in resid_ligand_list2:
			rl1, rl2 = resid_lig1[0], resid_lig2[0]
			if("LIG" in rl1 and "LIG" not in rl2):
				rlCombos.add((rl1, rl2))
			elif("LIG" not in rl1 and "LIG" in rl2):
				rlCombos.add((rl2, rl1))
	return rlCombos




# Considers all pairs hoh1 and hoh2 that have hydrogen bond between them. 
# If residlist1 and residlist2 are the list of residues that exhibit hydrogen 
# bond to hoh1 and hoh2 respectively, then return the cross product combination
# of all resid1,resid2 pairs from the two lists. 
def writeToExtendedWBondFile(fwb, frameWBond2Dicts):
	print("writeToExtendedWBondFile()")
	extendWaterBridges = set()
	waterToWaterDict, waterToResidueDict = frameWBond2Dicts
	print("Length WaterToWater", len(waterToWaterDict), len(waterToResidueDict))
	for water1 in waterToWaterDict:
		water2_list = sorted(waterToWaterDict[water1]) # solvent molecules close to water 1 
		wat1_resid_list = []
		if(water1 in waterToResidueDict):
			wat1_resid_list = waterToResidueDict[water1]
		for i, water2 in enumerate(water2_list):
			wat2_resid_list = []
			if(water2 in waterToResidueDict):
				wat2_resid_list = waterToResidueDict[water2]
			extendWaterBridges = extendWaterBridges | residueCombinations(water1, water2, wat1_resid_list, wat2_resid_list)
	for resid1, resid2 in extendWaterBridges:
		fwb.write(resid1 + " -- " + resid2 + "\n")

# Considers all pairs hoh1 and hoh2 that have hydrogen bond between them. 
# If residlist1 and residlist2 are the list of residues and ligand that exhibit hydrogen 
# bond to hoh1 and hoh2 respectively, then return the cross product combination
# of all resid1,resid2 pairs from the two lists. 
def writeToExtendedLWBondFile(fwb, frameLWBond2Dicts):
	extendLigandWaterBridges = set()
	waterToWaterDict, waterToResidueLigandDict = frameLWBond2Dicts
	print(len(waterToWaterDict), len(waterToResidueLigandDict))
	for water1 in waterToWaterDict:
		water2_list = sorted(waterToWaterDict[water1])
		wat1_resid_ligand_list = []
		if(water1 in waterToResidueLigandDict):
			wat1_resid_ligand_list = waterToResidueLigandDict[water1]
		for i, water2 in enumerate(water2_list):
			wat2_residue_ligand_list = []
			if(water2 in waterToResidueLigandDict):
				wat2_residue_ligand_list = waterToResidueLigandDict[water2]
			extendLigandWaterBridges = extendLigandWaterBridges | residueLigandCombinations(water1, water2, wat1_resid_ligand_list, wat2_residue_ligand_list)
	for ligand, resid in extendLigandWaterBridges:
		fwb.write(ligand + " -- " + resid + "\n")


# Parses the water mediated hydrogen bond frame line by line and writes to the new 
# file that represents hBondType file. It will extract the stride, traj, top, 
# nFrames, compute time information. Every time the frame number updates, then 
# the frameHBondDataStruct written to the new file and then cleared out to prepare
# for the new frame.  
def writeToHeader(f, line, nFrames, currFrame, computingTime, frameHBondDataStruct, hBondType):
	isHeader = True
	if("Stride:" in line or "TrajectoryPath:" in line or "TopologyPath:" in line):
		f.write(line)
	elif("nFrames:" in line):
		nFrames = int(line.split("nFrames:")[1])
	elif('Computing Time:' in line):
		computingTime = float(line.split("Computing Time:")[1])
	elif("Hydrogen Bond-Water Mediated:" in line and "Ligand" not in line):
		top_path = line.split("Hydrogen Bond-Water Mediated:")[1]
		if(hBondType == '-hbbb'): f.write("\nBackbone Backbone Hydrogen Bonds:" + top_path)
		if(hBondType == '-hbsb'): f.write("\nSidechain Backbond Hydrogen Bonds:" + top_path)
		if(hBondType == '-hbss'): f.write("\nSidechain Sidechain Hydrogen Bonds:" + top_path)
		if(hBondType == '-rw'): f.write("\nResidue Water Hydrogen Bonds:" + top_path)
		if(hBondType == '-wb'): f.write("\nWater Bonds:" + top_path)
		if(hBondType == '-wb2'): f.write("\nExtended Water Bonds:" + top_path)

	### Ligand Based Interactions
	elif("Ligand Hydrogen Bond-Water Mediated:" in line):
		top_path = line.split("Ligand Hydrogen Bond-Water Mediated:")[1]
		if(hBondType == "-hlb"): f.write("\nLigand Backbone Hydrogen Bonds:" + top_path)
		if(hBondType == "-hls"): f.write("\nLigand Sidechain Hydrogen Bonds:" + top_path)
		if(hBondType == '-lw'): f.write("\nLigand Water Hydrogen Bonds:" + top_path)
		if(hBondType == '-lwb'): f.write("\nLigand Water Bonds:" + top_path)
		if(hBondType == '-lwb2'): f.write("\nLigand Extended Water Bonds:" + top_path)

	### Protein H-Bond Interactions
	elif("Hydrogen Bond Water Frame:" in line and "Ligand" not in line):
		frame_index = line.split("Hydrogen Bond Water Frame: ")[1].strip()
		if(hBondType == '-wb'):
			writeToWBondFile(f, frameHBondDataStruct)
			frameHBondDataStruct = {}
		elif(hBondType == '-wb2'):
			writeToExtendedWBondFile(f, frameHBondDataStruct)
			frameHBondDataStruct = ({},{})
		elif(hBondType in ['-hbbb', '-hbsb', '-hbss', '-rw']):
			writeToHBondFile(f, frameHBondDataStruct)
			frameHBondDataStruct = []
		currFrame = int(frame_index)
		if(currFrame %100 == 0): print("Finished: " + str(currFrame) + " frame.")
		if(hBondType == '-hbbb'): f.write("\nBackbone Backbone Hydrogen Bond Frame:" + str(frame_index) +"\n")
		if(hBondType == '-hbsb'): f.write("\nSidechain Backbond Hydrogen Bond Frame:" + str(frame_index) +"\n")
		if(hBondType == '-hbss'): f.write("\nSidechain Sidechain Hydrogen Bond Frame:" + str(frame_index) +"\n")
		if(hBondType == '-rw'): f.write("\nResidue Water Hydrogen Bond Frame:" + str(frame_index) +"\n")
		if(hBondType == '-wb'): f.write("\nWater Bond Frame:" + str(frame_index) +"\n")
		if(hBondType == '-wb2'): f.write("\nExtended Water Bond Frame:" + str(frame_index) + "\n")
	
	### Ligand H-Bond Interactions
	elif("Ligand Hydrogen Bond Water Frame:" in line):
		frame_index = line.split("Ligand Hydrogen Bond Water Frame: ")[1].strip()
		if(hBondType == '-lwb'):
			writeToLWBondFile(f, frameHBondDataStruct, currFrame)
			frameHBondDataStruct = {}
		elif(hBondType == '-lwb2'):
			writeToExtendedLWBondFile(f, frameHBondDataStruct)
			frameHBondDataStruct = ({},{})
		elif(hBondType in ['-hlb', '-hls', '-lw']):
			writeToHBondFile(f, frameHBondDataStruct)
			frameHBondDataStruct = []
		
		currFrame = int(frame_index)
		if(hBondType == '-hlb'): f.write("\nLigand Backbone Hydrogen Bond Frame:" + str(frame_index) + "\n")
		if(hBondType == '-hls'): f.write("\nLigand Sidechain Hydrogen Bond Frame:" + str(frame_index) + "\n")
		if(hBondType == '-lw'): f.write("\nLigand Water Hydrogen Bond Frame:" + str(frame_index) +"\n")
		if(hBondType == '-lwb'): f.write("\nLigand Water Bond Frame:" + str(frame_index) +"\n")
		if(hBondType == '-lwb2'): f.write("\nLigand Extended Water Bond Frame:" + str(frame_index) + "\n")


	elif(line !="\n"): isHeader = False
	return isHeader, nFrames, currFrame, computingTime, frameHBondDataStruct
