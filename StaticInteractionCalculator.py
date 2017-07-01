# Author: Anthony Kai Kwang Ma
# Email: anthonyma27@gmail.com, akma327@stanford.edu
# StaticInteractionCalculator.py

from __future__ import print_function
import os
import sys
import time
import itertools
import mdtraj as md
from waterMediatedHBonds import *
from hbond import baker_hubbard
from saltbridge import *
from vanderwaal import *
from hydrogenbond import *
from hydrophobe import *
from pication import *
from aromatic import *
from vmdMeasureHBonds_noang import * 


USAGE_STR = """ 
# Usage:
# python StaticInteractionCalculator.py <TOP_PATH> <OUTPUT_FILE_PATH> <-pdb flag> <PDB CODE> <-append flag> <optional -solv flag> <solvent> <optional -chain> <chainID> <optional -ligand> <Ligand> <-interlist> <INTERACTION LIST>

# Arguments:
# <TOP_PATH> Absolute path to the topology for calculating the interactions in crystal structure 
# <OUTPUT_FILE_PATH> Absolute path to output file containing the interactions in crystal structure
# <-pdb> Optional flag to denote the specific pdb code for the topology being used. 
# <-append> Optional flag to denote that you are appending to an existing file instead of writing a new file. 
# This is useful for the case of first calculating non hydrogen bonds and then appending hydrogen bond calculations
# <-solv> Optional flag to denote the solvent to use 
# <-chain> Optional flag to denote specific chain value 
# <optional -ligand flag> To denote the resname of ligand in the simulation.
# <Interaction List> List of interactions to compute

# Example 1:
TOP_PATH="/scratch/PI/rondror/akma327/DynamicNetworks/data/StaticNetworks/pdb-files/4EIY_edited_Hadded.pdb"
OUTPUT_FILE_PATH="/scratch/PI/rondror/akma327/DynamicNetworks/data/StaticNetworks/static-interaction-output/4EIY.txt"
python StaticInteractionCalculator.py $TOP_PATH $OUTPUT_FILE_PATH -interlist -all

"""

K_MIN_ARG = 3

solvent = None 
solventId = None 
backbone_atoms = ['N', 'O']
ALL_INTERACTION_TYPES = ['-sb', '-pc', '-ps', '-ts', '-vdw', '-hbbb','-hbsb', '-hbss', '-rw', '-wb', '-wb2']


# Create file descriptor for writing to output file 
def createFileWriter(OUTPUT_FILE_PATH, APPEND_FLAG):
	base_directory = "/".join(OUTPUT_FILE_PATH.split("/")[:-1])
	if not os.path.exists(base_directory):
		os.makedirs(base_directory)
	if(APPEND_FLAG == True):
		f = open(OUTPUT_FILE_PATH, 'a')
	else:
		f = open(OUTPUT_FILE_PATH, 'w')
	return f, OUTPUT_FILE_PATH


# Generate combinations of first order water bridges 
def getWBondCombos(water, resid_list):
	wBondCombos = set()
	for pair in itertools.combinations(resid_list, 2):
		resid1, resid2 = pair 
		if(resid1 != resid2):
			wBondCombos.add((resid1, resid2))
	return wBondCombos

# Takes the frameWBondDict which contains key = unique solvent identifier, and 
# value = list of all the residues that the unique solvent has h-bond with 
# at the particular frame and writes out all the combinations 
def writeToWBondFile(fwb, frameWBondDict):
	print("wt", frameWBondDict)
	frameWBonds = set() 
	for water in frameWBondDict:
		resid_list = sorted(frameWBondDict[water])
		frameWBonds = frameWBonds | getWBondCombos(water, resid_list)
	frameWBonds = sorted(list(frameWBonds))
	print(frameWBonds)
	for resid1, resid2 in frameWBonds:
		fwb.write(resid1 + " -- " + resid2 + "@-wb\n")


# Generate combinations of first order ligand water bridges
def getLWBondCombos(ligand, resid_ligand_list):
	lwBondCombos = set()
	for resid in resid_ligand_list:
		if("LIG" not in resid): # valid residue
			lwBondCombos.add((ligand, resid))
	return lwBondCombos

# Takes the frameLWBond Dict which contains key = unique solvent identifier, and
# value = list of all residues and ligand that the unique solvent has h-bond with
# at the particular frame and writes out all the valid combinations involving
# ligand -- residue
def writeToLWBondFile(fwb, frameLWBondDict):
	def find_ligand(arr):
		"""
			Subroutine to find ligand in list of residue + ligand
		"""
		for elem in arr:
			cand_ligand = elem
			if("LIG" in cand_ligand): return cand_ligand
		return None

	frameLWBonds = set()
	for water in frameLWBondDict:
		resid_ligand_list = sorted(frameLWBondDict[water])
		ligand = find_ligand(resid_ligand_list)
		if(ligand == None): continue
		frameLWBonds = frameLWBonds | getLWBondCombos(ligand, resid_ligand_list)
	frameLWBonds = sorted(list(frameLWBonds))
	for ligand, resid in frameLWBonds:
		fwb.write(ligand + " -- " + resid + "@-lwb\n")


# Takes two lists of residues and returns cross product combination of the residue pairs
def residueCombinations(water1, water2, resid_list1, resid_list2):	
	return set([(res1, res2) for res1 in resid_list1 for res2 in resid_list2]) #doesn't make distinct by water 

def residueLigandCombinations(water1, water2, resid_ligand_list1, resid_ligand_list2):
	rlCombos = set()
	for rl1 in resid_ligand_list1:
		for rl2 in resid_ligand_list2:
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
	extendWaterBridges = set()
	waterToWaterDict, waterToResidueDict = frameWBond2Dicts
	print(len(waterToWaterDict), len(waterToResidueDict))
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
		fwb.write(resid1 + " -- " + resid2 + "@-wb2\n")


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
		fwb.write(ligand + " -- " + resid + "@-lwb2\n")

# Computation for Salt Bridges
def calcSaltBridgeResults(allFrames, f, TOP_PATH):
	print("\n\nSalt Bridges:" + TOP_PATH + "\n")
	tic = time.clock()
	chainDict = initSaltBridgeChainDict(allFrames)
	sbFramePairs = calcSaltBridgeFramePairs(allFrames, chainDict)
	toc = time.clock()
	computingTime = toc - tic
	for index, sbPairs in enumerate(sbFramePairs):
		for pairs in sbPairs:
			atom1, atom2 = pairs[0], pairs[1]
			f.write(str(atom1[0]) + "_" + str(atom1[1]) +  " -- " + str(atom2[0]) + "_" + str(atom2[1]) + "@-sb\n")
	print("\nComputing Time:" + str(computingTime) + "\n")
	return computingTime



#Computation for pi-cation interaction
def calcPiCationResults(allFrames, f, TOP_PATH):
	print("\n\nPi-Cation:" + TOP_PATH + "\n")
	tic = time.clock()
	chainDict = initPiCationChainDict(allFrames)
	pcFramePairs = calcPiCationFramePairs(allFrames, chainDict)
	toc = time.clock()
	computingTime = toc - tic 
	for index, picationPairs in enumerate(pcFramePairs):
		for pair in picationPairs:
			cation, aromatic = pair[0], pair[1]
			cation_str = str(cation[0]) + "_" + str(cation[1])
			aromatic_str = str(aromatic[0][0]).split("-")[0] + "_" + str(aromatic[0][1])
			f.write(cation_str + " -- " + aromatic_str + "@-pc\n")
	print("\nComputing Time:" + str(computingTime) + "\n")
	return computingTime


# Compute pi stacking interactions
def calcFaceToFaceResults(allFrames, f, TOP_PATH):
	print("\n\nPi-Stacking:" + TOP_PATH + "\n")
	nFrames = len(allFrames)
	tic = time.clock()
	chainDict = initFaceFaceAromaticChainDict(allFrames)
	psFramePairs = calcPiStackingFramePairs(allFrames, chainDict)
	toc = time.clock()
	computingTime = toc - tic 
	for index, pi_stackPairs in enumerate(psFramePairs):
		for pair in pi_stackPairs:
			aromatic1, aromatic2 = pair[0], pair[1]
			aromatic_str1 = str(aromatic1[0][0]).split("-")[0] + "_" + str(aromatic1[0][1])
			aromatic_str2 = str(aromatic2[0][0]).split("-")[0] + "_" + str(aromatic2[0][1])
			f.write(aromatic_str1 + " -- " + aromatic_str2 + "@-ps\n")
	print("\nComputing Time:" + str(computingTime) + "\n")
	return computingTime


#Computation for edge to face aromatic detection
def calcEdgeToFaceResults(allFrames, f, TOP_PATH):
	print("\n\nT-Stacking:" + TOP_PATH + "\n")
	nFrames = len(allFrames)
	tic = time.clock()
	tsFramePairs = []
	chainDict = initFaceEdgeAromaticChainDict(allFrames)
	tsFramePairs = calcTStackingFramePairs(allFrames, chainDict)
	toc = time.clock()
	computingTime = toc - tic
	for index, t_stackPairs in enumerate(tsFramePairs):
		for pair in t_stackPairs:
			aromatic1, aromatic2 = pair[0], pair[1]
			aromatic_str1 = str(aromatic1[0][0]).split("-")[0] + "_" + str(aromatic1[0][1])
			aromatic_str2 = str(aromatic2[0][0]).split("-")[0] + "_" + str(aromatic2[0][1])
			f.write(aromatic_str1 + " -- " + aromatic_str2 + "@-ts\n")
	print("\nComputing Time:" + str(computingTime) + "\n")
	return computingTime

# Calculate van der waals interactions 
def calcVanderwaalsResults(allFrames, f, TOP_PATH):
	print("\n\nVan Der Waals:" + TOP_PATH + "\n")
	tic = time.clock()
	vdwFramePairs = calcVDWFramePairs(allFrames)
	toc = time.clock()
	computingTime = toc - tic
	for index, vanderwaalPairs in enumerate(vdwFramePairs):
		for pair in vanderwaalPairs:
			atom1, atom2 = pair[0], pair[1]
			f.write(str(atom1[0]) + "_" + str(atom1[1]) +  " -- " + str(atom2[0]) + "_" + str(atom2[1]) + "@-vdw\n")
	print("\nComputing Time:" + str(computingTime) + "\n")
	return computingTime


# Calculate raw water mediated hydrogen bond interactions 
def calcHydrogenBondWaterResults(file_writer, TOP_PATH, solvent_id, chainId=None, ligand=None):
	print("\n\nHydrogen Bond-Water Mediated:" + TOP_PATH + "\n")
	tic = time.clock()
	hbwFramePairs = calcVMDHydrogenBondWaterResults(file_writer, 1, TOP_PATH, None, solvent_id, chainId, None, True)
	for atom1, atom2 in hbwFramePairs:
		file_writer.write(atom1 + " -- " + atom2 + "@-hbw\n")
	toc = time.clock()
	computingTime = toc - tic 
	print("\nComputing Time:" + str(computingTime) + "\n")
	return computingTime


def calcLigandHydrogenBondWaterResults(file_writer, TOP_PATH, solvent_id, chainId=None, ligand=None):
	print("\n\nLigand Hydrogen Bond-Water Mediated:" + TOP_PATH + "\n")
	tic = time.clock()
	lhbwFramePairs = calcVMDHydrogenBondWaterResults(file_writer, 1, TOP_PATH, None, solvent_id, chainId, ligand, True)
	for atom1, atom2 in lhbwFramePairs:
		file_writer.write(atom1 + " -- " + atom2 + "@-lhbw\n")
	toc = time.clock()
	computingTime = toc - tic 
	print("\nComputing Time:" + str(computingTime) + "\n")
	return computingTime

# Calculate backbone to backbone interactions
def calcBackBoneHBondResults(file_writer, TOP_PATH, solvent_id, chainId=None, ligand=None):
	print("\n\nBackbone Backbone Hydrogen Bonds:" + TOP_PATH + "\n")
	tic = time.clock()
	hbwFramePairs = calcVMDHydrogenBondWaterResults(file_writer, 1, TOP_PATH, None, solvent_id, chainId, None, crystal =True)
	print(hbwFramePairs)
	for resAtom1, resAtom2 in hbwFramePairs:
		resAtom1, resAtom2 = str(resAtom1), str(resAtom2)
		if("HOH" not in resAtom1 and "HOH" not in resAtom2):
			atom1, atom2 = resAtom1.split("-")[1].strip(), resAtom2.split("-")[1].strip()
			if(atom1 in backbone_atoms and atom2 in backbone_atoms):
				file_writer.write(str(resAtom1) + " -- " + str(resAtom2) + "@-hbbb\n")
	toc = time.clock()
	computingTime = toc - tic 
	print("\nComputing Time:" + str(computingTime) + "\n")
	return computingTime


def calcLigandBackBoneHBondResults(file_writer, TOP_PATH, solvent_id, chainId=None, ligand=None):
	print("\n\nLigand Backbone Hydrogen Bonds:" + TOP_PATH + "\n")
	tic = time.clock()
	lhbwFramePairs = calcVMDHydrogenBondWaterResults(file_writer, 1, TOP_PATH, None, solvent_id, chainId, ligand, crystal =True)
	print(lhbwFramePairs)
	for resAtom1, resAtom2 in lhbwFramePairs:
		resAtom1, resAtom2 = str(resAtom1), str(resAtom2)
		if("LIG" in resAtom1 and ("LIG" not in resAtom2 and "HOH" not in resAtom2)): 
			atom2 = resAtom2.split("-")[1].strip()
			if(atom2 in backbone_atoms):
				file_writer.write(resAtom1 + " -- " + resAtom2 + "@-hlb\n")
		if(("LIG" not in resAtom1 and "HOH" not in resAtom1) and "LIG" in resAtom2):
			atom1 = resAtom1.split("-")[1].strip()
			if(atom1 in backbone_atoms):
				file_writer.write(resAtom1 + " -- " + resAtom2 + "@-hlb\n")
	toc = time.clock()
	computingTime = toc - tic
	print("\nComputing Time:" + str(computingTime) + "\n")
	return computingTime


# Calculate sidechain to backbone interactions 
def calcSideChainBackBoneHBondResults(file_writer, TOP_PATH, solvent_id, chainId=None, ligand=None):
	print("\n\nSidechain Backbond Hydrogen Bonds:" + TOP_PATH + "\n")
	tic = time.clock()
	hbwFramePairs = calcVMDHydrogenBondWaterResults(file_writer, 1, TOP_PATH, None, solvent_id, chainId, None, crystal=True)
	for resAtom1, resAtom2 in hbwFramePairs:
		resAtom1, resAtom2 = str(resAtom1), str(resAtom2)
		if("HOH" not in resAtom1 and "HOH" not in resAtom2):
			atom1, atom2 = resAtom1.split("-")[1].strip(), resAtom2.split("-")[1].strip()
			if((atom1 in backbone_atoms and atom2 not in backbone_atoms) or (atom1 not in backbone_atoms and atom2 in backbone_atoms)):
				file_writer.write(str(resAtom1) + " -- " + str(resAtom2) + "@-hbsb\n")
	toc = time.clock()
	computingTime = toc - tic 
	print("\nComputing Time:" + str(computingTime) + "\n")
	return computingTime


## BUG DIDNT EVEN CHECK IF LIGAND IS PART OF ONE OF THE ATOMS
def calcLigandSideChainHBondResults(file_writer, TOP_PATH, solvent_id, chainId=None, ligand=None):
	print("\n\nLigand Sidechain Backbond Hydrogen Bonds:" + TOP_PATH + "\n")
	tic = time.clock()
	lhbwFramePairs = calcVMDHydrogenBondWaterResults(file_writer, 1, TOP_PATH, None, solvent_id, chainId, ligand, crystal=True)
	for resAtom1, resAtom2 in lhbwFramePairs:
		resAtom1, resAtom2 = str(resAtom1), str(resAtom2)
		if("LIG" in resAtom1 and ("LIG" not in resAtom2 and "HOH" not in resAtom2)): 
			atom2 = resAtom2.split("-")[1].strip()
			if(atom2 not in backbone_atoms):
				file_writer.write(resAtom1 + " -- " + resAtom2 + "@-hls\n")
		if(("LIG" not in resAtom1 and "HOH" not in resAtom1) and "LIG" in resAtom2):
			atom1 = resAtom1.split("-")[1].strip()
			if(atom1 not in backbone_atoms):
				file_writer.write(resAtom1 + " -- " + resAtom2 + "@-hls\n")
	toc = time.clock()
	computingTime = toc - tic 
	print("\nComputing Time:" + str(computingTime) + "\n")
	return computingTime


# Calculate sidechain to sidechain interactions 
def calcSideChainHBondResults(file_writer, TOP_PATH, solvent_id, chainId=None, ligand=None):
	print("\n\nSidechain Sidechain Hydrogen Bonds:" + TOP_PATH + "\n")
	tic = time.clock()
	hbwFramePairs = calcVMDHydrogenBondWaterResults(file_writer, 1, TOP_PATH, None, solvent_id, chainId, None, crystal=True)
	for resAtom1, resAtom2 in hbwFramePairs:
		resAtom1, resAtom2 = str(resAtom1), str(resAtom2)
		if("HOH" not in resAtom1 and "HOH" not in resAtom2):
			atom1, atom2 = resAtom1.split("-")[1].strip(), resAtom2.split("-")[1].strip()
			if(atom1 not in backbone_atoms and atom2 not in backbone_atoms):
				file_writer.write(str(resAtom1) + " -- " + str(resAtom2) + "@-hbss\n")
	toc = time.clock()
	computingTime = toc - tic 
	print("\nComputing Time:" + str(computingTime) + "\n")
	return computingTime

# Calculate water to residue interactions
def calcResidueWaterBondResults(file_writer, TOP_PATH, solvent_id, chainId=None, ligand=None):
	print("\n\nResidue Water Hydrogen Bonds:" + TOP_PATH + "\n")
	tic = time.clock()
	hbwFramePairs = calcVMDHydrogenBondWaterResults(file_writer, 1, TOP_PATH, None, solvent_id, chainId, None, crystal=True)
	for resAtom1, resAtom2 in hbwFramePairs:
		resAtom1, resAtom2 = str(resAtom1), str(resAtom2)
		if(("HOH" in resAtom1 and "HOH" not in resAtom2) or ("HOH" not in resAtom1 and "HOH" in resAtom2)):
			file_writer.write(str(resAtom1) + " -- " + str(resAtom2) + "@-rw\n")
	toc = time.clock()
	computingTime = toc - tic 
	print("\nComputing Time:" + str(computingTime) + "\n")
	return computingTime

def calcLigandToWaterBondResults(file_writer, TOP_PATH, solvent_id, chainId=None, ligand=None):
	print("\n\nLigand Water Hydrogen Bonds:" + TOP_PATH + "\n")
	tic = time.clock()
	lhbwFramePairs = calcVMDHydrogenBondWaterResults(file_writer, 1, TOP_PATH, None, solvent_id, chainId, ligand, crystal=True)
	for resAtom1, resAtom2 in lhbwFramePairs:
		resAtom1, resAtom2 = str(resAtom1), str(resAtom2)
		if(("HOH" in resAtom1 and "LIG" in resAtom2) or ("LIG" in resAtom1 and "HOH" in resAtom2)):
			file_writer.write(resAtom1 + " -- " + resAtom2 + "@-lw\n")
	toc = time.clock()
	computingTime = toc - tic 
	print("\nComputing Time:" + str(computingTime) + "\n")
	return computingTime


# Calculate first degree water bridges
def calcWaterBondResults(file_writer, TOP_PATH, solvent_id, chainId=None, ligand=None):

	print("calcWaterBondResults()", solvent_id, solventId)

	print("\n\nWater Bonds:" + TOP_PATH + "\n")
	frameWBondDict = {} # THIS IS EMPTY THAT'S WHY
	tic = time.clock()
	hbwFramePairs = calcVMDHydrogenBondWaterResults(file_writer, 1, TOP_PATH, None, solvent_id, chainId, None, crystal=True)
	print("cakcWaterBondResults() frame pair", hbwFramePairs)
	for resAtom1, resAtom2 in hbwFramePairs:
		resAtom1, resAtom2 = str(resAtom1), str(resAtom2)
		if('HOH' in resAtom1 and 'HOH' in resAtom2): continue 
		if(not('HOH' not in resAtom1 and 'HOH' not in resAtom2)):
			if('HOH' in resAtom1):
				solvent_atom, resid_atom = resAtom1, resAtom2
			else:
				solvent_atom, resid_atom = resAtom2, resAtom1
			if(solvent_atom not in frameWBondDict):
				frameWBondDict[solvent_atom] = [resid_atom]
			else:
				frameWBondDict[solvent_atom].append(resid_atom)

	print("calcWaterBondResults() Result", frameWBondDict)
	writeToWBondFile(file_writer, frameWBondDict)


def calcLigandWaterBondResults(file_writer, TOP_PATH, solvent_id, chainId=None, ligand=None):
	print("\n\nLigand Water Bonds:" + TOP_PATH + "\n")
	frameLWBondDict = {}
	tic = time.clock()
	lhbwFramePairs = calcVMDHydrogenBondWaterResults(file_writer, 1, TOP_PATH, None, solvent_id, chainId, ligand, crystal=True)
	for resAtom1, resAtom2 in lhbwFramePairs:
		resAtom1, resAtom2 = str(resAtom1), str(resAtom2)
		if('HOH' in resAtom1 and 'HOH' in resAtom2): continue 
		if(not('HOH' not in resAtom1 and 'HOH' not in resAtom2)):
			if('HOH' in resAtom1):
				solvent_atom, resid_ligand_atom = resAtom1, resAtom2
			else:
				solvent_atom, resid_ligand_atom = resAtom2, resAtom1
			if(solvent_atom not in frameLWBondDict):
				frameLWBondDict[solvent_atom] = [resid_ligand_atom]
			else:
				frameLWBondDict[solvent_atom].append(resid_ligand_atom)
	print("waht", frameLWBondDict)
	writeToLWBondFile(file_writer, frameLWBondDict)


# This method computes second degree water bridges. This class of interaction
# is defined between amino acids aa1 and aa2, which have unique hydrogen bond
# interactions with hoh1 and hoh2 respectively. If hoh1-hoh2 also have 
# hydrogen bond interaction then aa1 and aa2 form an extended second degree
# water bridge. 
def calcWaterBond2Results(file_writer, TOP_PATH, solvent_id, chainId=None, ligand=None):
	frameWBond2Dicts = ({}, {})
	hbwFramePairs = calcVMDHydrogenBondWaterResults(file_writer, 1, TOP_PATH, None, solvent_id, chainId, None, crystal=True)
	for resAtom1, resAtom2 in hbwFramePairs:
		if("HOH" in resAtom1 and "HOH" in resAtom2):
			if(resAtom1 not in frameWBond2Dicts[0]):
				frameWBond2Dicts[0][resAtom1] = [resAtom2]
			else:
				frameWBond2Dicts[0][resAtom1].append(resAtom2)
		elif(("HOH" in resAtom1 and "HOH" not in resAtom2) or ("HOH" not in resAtom1 and "HOH" in resAtom2)):
			if("HOH" in resAtom1): solvent_atom, resid_atom = resAtom1, resAtom2
			else: solvent_atom, resid_atom = resAtom2, resAtom1
			if(solvent_atom not in frameWBond2Dicts[1]):
				frameWBond2Dicts[1][solvent_atom] = [resid_atom]
			else: 
				frameWBond2Dicts[1][solvent_atom].append(resid_atom)
	writeToExtendedWBondFile(file_writer, frameWBond2Dicts)


def calcLigandWaterBond2Results(file_writer, TOP_PATH, solvent_id, chainId=None, ligand=None):
	frameLWBond2Dicts = ({}, {})
	lhbwFramePairs = calcVMDHydrogenBondWaterResults(file_writer, 1, TOP_PATH, None, solvent_id, chainId, ligand, crystal=True)
	for resAtom1, resAtom2 in lhbwFramePairs:
		if("HOH" in resAtom1 and "HOH" in resAtom2):
			if(resAtom1 not in frameLWBond2Dicts[0]):
				frameLWBond2Dicts[0][resAtom1] = [resAtom2]
			else:
				frameLWBond2Dicts[0][resAtom1].append(resAtom2)
		elif(("HOH" in resAtom1 and "HOH" not in resAtom2) or ("HOH" not in resAtom1 and "HOH" in resAtom2)):
			if("HOH" in resAtom1): solvent_atom, resid_atom = resAtom1, resAtom2
			else: solvent_atom, resid_atom = resAtom2, resAtom1 
			if(solvent_atom not in frameLWBond2Dicts[1]):
				frameLWBond2Dicts[1][solvent_atom] = [resid_atom]
			else:
				frameLWBond2Dicts[1][solvent_atom].append(resid_atom)
	writeToExtendedLWBondFile(file_writer, frameLWBond2Dicts)



# Calculate static interactions 
def calcStaticInteractions(file_writer, interaction_selection, TOP_PATH, solvent_id, chainId=None, ligand=None):
	print("calcStaticInteractions()", solvent_id)
	if(interaction_selection in ['-sb', '-pc', '-ps', '-ts', '-vdw']):
		traj = md.load(TOP_PATH)
	if(interaction_selection == '-sb'):
		calcSaltBridgeResults(traj, file_writer, TOP_PATH)
	if(interaction_selection == '-pc'):
		calcPiCationResults(traj, file_writer, TOP_PATH)
	if(interaction_selection == '-ps'):
		calcFaceToFaceResults(traj, file_writer, TOP_PATH)
	if(interaction_selection == '-ts'):
		calcEdgeToFaceResults(traj, file_writer, TOP_PATH)
	if(interaction_selection == '-vdw'):
		calcVanderwaalsResults(traj, file_writer, TOP_PATH)
	if(interaction_selection == '-hbw'):
		calcHydrogenBondWaterResults(file_writer, TOP_PATH, solvent_id, chainId, ligand)
	if(interaction_selection == '-hbbb'):
		calcBackBoneHBondResults(file_writer, TOP_PATH, solvent_id, chainId, ligand)
	if(interaction_selection == '-hbsb'):
		calcSideChainBackBoneHBondResults(file_writer, TOP_PATH, solvent_id, chainId, ligand)
	if(interaction_selection == '-hbss'):
		calcSideChainHBondResults(file_writer, TOP_PATH, solvent_id, chainId, ligand)
	if(interaction_selection == '-rw'):
		calcResidueWaterBondResults(file_writer, TOP_PATH, solvent_id, chainId, ligand)
	if(interaction_selection == '-wb'):
		calcWaterBondResults(file_writer, TOP_PATH, solvent_id, chainId, ligand)
	if(interaction_selection == '-wb2'):
		calcWaterBond2Results(file_writer, TOP_PATH, solvent_id, chainId, ligand)

	### Ligand Based Interactions
	if(interaction_selection == "-lhbw"):
		calcLigandHydrogenBondWaterResults(file_writer, TOP_PATH, solvent_id, chainId, ligand)
	if(interaction_selection == "-hlb"):
		calcLigandBackBoneHBondResults(file_writer, TOP_PATH, solvent_id, chainId, ligand)
	if(interaction_selection == "-hls"):
		calcLigandSideChainHBondResults(file_writer, TOP_PATH, solvent_id, chainId, ligand)
	if(interaction_selection == "-lw"):
		calcLigandToWaterBondResults(file_writer, TOP_PATH, solvent_id, chainId, ligand)
	if(interaction_selection == "-lwb"):
		calcLigandWaterBondResults(file_writer, TOP_PATH, solvent_id, chainId, ligand)
	if(interaction_selection == "-lwb2"):
		calcLigandWaterBond2Results(file_writer, TOP_PATH, solvent_id, chainId, ligand)


if __name__ == "__main__":
	if(len(sys.argv) < K_MIN_ARG):
		print(USAGE_STR)
		exit(0)
	TOP_PATH = sys.argv[1]
	OUTPUT_FILE_PATH = sys.argv[2]
	APPEND_FLAG = False 
	PDB_CODE = None
	if('-pdb' in sys.argv): PDB_CODE = sys.argv[sys.argv.index('-pdb') + 1]
	if('-append' in sys.argv): APPEND_FLAG = True 
	solvent_id = None
	if('-solv' in sys.argv):
		solvent_index = sys.argv.index('-solv') + 1
		solvent_id = sys.argv[solvent_index]
	print("main()", solvent_id)
	chainId = None 
	if("-chain" in sys.argv):
		chainId = str(sys.argv[sys.argv.index("-chain") + 1])
	ligand = None
	if("-ligand" in sys.argv):
		ligand = str(sys.argv[sys.argv.index("-ligand") + 1])
	inter_sel_index = sys.argv.index('-interlist') + 1
	interaction_selections = sys.argv[inter_sel_index:]
	if("-all" in interaction_selections):
		interaction_selections = ALL_INTERACTION_TYPES
	file_writer, outputPath = createFileWriter(OUTPUT_FILE_PATH, APPEND_FLAG)
	if(PDB_CODE != None): file_writer.write("PDB:" + PDB_CODE + "\n")
	for interaction_selection in interaction_selections:
		print("Processing Interaction " + interaction_selection, solvent_id)
		calcStaticInteractions(file_writer, interaction_selection, TOP_PATH, solvent_id, chainId, ligand)

