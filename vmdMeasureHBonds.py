# Author: Anthony Kai Kwang Ma
# Email: akma327@stanford.edu

from __future__ import division

import vmd, molecule
from vmd import *
import sys
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import ast 
import time 

USAGE_STR = """
# Purpose:
# This script utilizes VMD's measure command to extract all the hydrogen bonds within a specified trajectory 

# Usage: python vmdMeasureHBonds-ExternalLookup.py  <TRAJECTORY_PATH> <TOPOLOGY_PATH> <OUTPUT_FILE> <Optional -stride flag> <Stride> <Optional -chain flag> <Chain> <Optional -solv flag> <solvent_id> 

# Arguments:
# <TRAJECTORY_PATH> Absolute path to trajectory to compute hydrogen bond contacts for 
# <TOPOLOGY_PATH> Absolute path to topology file 
# <OUTPUT_FILE> Absolute path to output file 
# <-stride> Optional flag to denote stride value other than default of 1
# <-chain> Optional flag to denote specific chain to compute hydrogen bonds for 
# <-solv> Optional flag to denote solvent other than default of TIP3

# Example:
# TRAJ_PATH="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutput/B2AR-active-Gs-BI-science2015/condition-22/rep_1/Prod_0/Prod_0_rewrapped.dcd"
# TOP_PATH="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutput/B2AR-active-Gs-BI-science2015/condition-22/step5_assembly.mae"
# OUTPUT_FILE="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutput/B2AR-active-Gs-BI-science2015/condition-22/rep_1/Prod_0/HB_3.5A_70D/hydrogen_bond_water_result.txt"
# python vmdMeasureHBonds-ExternalLookup.py $TRAJ_PATH $TOP_PATH $OUTPUT_FILE -chain R 

# Example VMD:
# import sys
# TRAJ_PATH="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutput/B2AR-active-Gs-BI-science2015/condition-22/rep_1/Prod_0/Prod_0_rewrapped.dcd"
# TOP_PATH="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutput/B2AR-active-Gs-BI-science2015/condition-22/step5_assembly.mae"
# OUTPUT_FILE="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutput/B2AR-active-Gs-BI-science2015/condition-22/rep_1/Prod_0/HB_3.5A_70D/hydrogen_bond_water_result.txt"
# sys.argv = ['a', TRAJ_PATH, TOP_PATH, OUTPUT_FILE, '-chain', 'R']
# execfile('vmdMeasureHBonds-ExternalLookup.py')

"""

K_MIN_ARG = 4
start=0
stop=-1
stride=1
smoothing=0
trajid = None
solvent_id = "TIP3"
lookup_table = {}


def createFileWriter(OutputPath):
	"""
		Generate file descriptor for path 
	"""
	filename = OutputPath
	print("Output File Path: " + filename)
	if not os.path.exists(os.path.dirname(filename)):
		os.makedirs(os.path.dirname(filename))
	f = open(filename, 'w')
	return f


def debugRepeats(frame_index, dist_cutoff, angle_cutoff):
	"""
		Script that checks whether there are any duplicates in the interaction key output 
	"""
	donors, acceptors = calcDonorAcceptors(frame_index, None, dist_cutoff, angle_cutoff)
	donors, acceptors = removeDuplicateIndices(donors, acceptors)
	key_pairs = []
	for i, d in enumerate(donors):
		key_pairs.append((donors[i], acceptors[i]))
	key_dict = {}
	for a in key_pairs:
		if(a not in key_dict):
			key_dict[a] = 1
		else:
			key_dict[a] += 1
	for k in key_dict:
		v = key_dict[k]
		if(v > 1):
			print(k,v)
	print(len(key_dict))


def removeDuplicateIndices(donors, acceptors):
	"""
		Removes duplicate indice pairs that represent the same interaction key
	"""
	pairs = []
	for i, d in enumerate(donors):
		pairs.append((donors[i], acceptors[i]))
	pairs = sorted(list(set(pairs)))
	newDonors, newAcceptors = [], []
	for p in pairs:
		newDonors.append(p[0])
		newAcceptors.append(p[1])
	return newDonors, newAcceptors



def getFileType(fileName):
	"""
		Extract suffix file type 
	"""

	file_type = fileName.split(".")[-1].strip()
	if(file_type == "nc"): file_type = 'netcdf'
	return file_type


def load_traj(top_path, traj_path):
	"""
		Load trajectory and topology into VMD display
	"""
	top_file_type = getFileType(top_path)
	trajid = molecule.load(top_file_type, top_path)
	if(traj_path != None):
		traj_file_type = getFileType(traj_path)
		molecule.read(trajid, traj_file_type, traj_path, beg=start, end=stop, skip=stride, waitfor=-1)


def indexToLabel(frame_index, hbond_index_list, crystal=False):
	"""
		Convert list of atom indices for a frame to atom label that includes information 
		regarding resname, resid, atom_name, and chain in the format of 
		<resname><resid>-<atom_name>_<chain> (ie GLU142-NZ_A)
	"""
	global lookup_table
	print("beg indexToLabel", len(lookup_table))
	label_list = []
	for display_index, atom_index in enumerate(hbond_index_list):
		if(atom_index in lookup_table):
			val = lookup_table[atom_index]
			label_list.append(val)
		else:
			if(display_index % 400 == 0): print("Adding key " + str(atom_index) + " to lookup_table")
			key = int(atom_index)
			aselect = evaltcl("set sel [atomselect top \"index " + str(atom_index) + "\"]")
			resname = evaltcl('$sel get resname')
			resid = evaltcl('$sel get resid')
			atom_name = evaltcl('$sel get name')
			chain = evaltcl('$sel get chain')
			value = resname + resid + "-" + atom_name + "_" + chain
			evaltcl('$sel delete')
			lookup_table[key] = value 
			label_list.append(lookup_table[atom_index])
	return label_list


# Generate atomselection command 
def genHBondsAtomSelection(frame_index, dist_cutoff, angle_cutoff, chainId, ligand, solvent_id):
	"""
		Query for hydrogen bonds involving water near the small molecule peptide and GPCR
	"""

	print("genHBondsAtomSelection() " + str(ligand))
	if(ligand == None):
		if(solvent_id == "IP3"):
			if(chainId == None):
				print("try this01")
				x = evaltcl("measure hbonds " + str(dist_cutoff) + " " + str(angle_cutoff) + " [atomselect top \"(resname IP3 and within 3.5 of protein) or protein and not lipid and not carbon\" frame " + str(frame_index) + "]")
			else:
				print("try this02")
				x = evaltcl("measure hbonds " + str(dist_cutoff) + " " + str(angle_cutoff) + " [atomselect top \"(resname IP3 and within 3.5 of protein and chain " + str(chainId) + ") or protein and chain " + str(chainId) + " and not lipid and not carbon\" frame " + str(frame_index) + "]")
		else:
			if(chainId == None):
				print("try this03")
				x = evaltcl("measure hbonds " + str(dist_cutoff) + " " + str(angle_cutoff) + " [atomselect top \"(water and within 3.5 of protein) or protein and not lipid and not carbon\" frame " + str(frame_index) + "]")
			else:
				print("try this04")
				x = evaltcl("measure hbonds " + str(dist_cutoff) + " " + str(angle_cutoff) + " [atomselect top \"(water and within 3.5 of protein and chain " + str(chainId) + ") or protein and chain " + str(chainId) + " and not lipid and not carbon\" frame " + str(frame_index) + "]")	
	else:
		if(solvent_id == "IP3"):
			if(chainId == None):
				print("try this1")
				x = evaltcl("measure hbonds " + str(dist_cutoff) + " " + str(angle_cutoff) + " [atomselect top \"(resname IP3 and within 12 of resname " + str(ligand) + ") or (not carbon and not sulfur and protein within 12 of resname " + str(ligand) + ") or (not carbon and not sulfur and resname " + str(ligand) + ") and (not lipid)\" frame " + str(frame_index) + "]")
			else:
				print("try this2")
				x = evaltcl("measure hbonds " + str(dist_cutoff) + " " + str(angle_cutoff) + " [atomselect top \"(resname IP3 and within 12 of resname " + str(ligand) + ") or (not carbon and not sulfur and protein and chain " + str(chainId) + " and within 12 of resname " + str(ligand) + ") or (not carbon and not sulfur and resname " + str(ligand) + ") and (not lipid)\" frame " + str(frame_index) + "]")
		else:
			if(chainId == None):
				print("try this3")
				x = evaltcl("measure hbonds " + str(dist_cutoff) + " " + str(angle_cutoff) + " [atomselect top \"(water within 12 of resname " + str(ligand) + ") or (not carbon and not sulfur and protein within 12 of resname " + str(ligand) + ") or (not carbon and not sulfur and resname " + str(ligand) + ") and (not lipid)\" frame " + str(frame_index) + "]")
			else:
				print("try this4")
				x = evaltcl("measure hbonds " + str(dist_cutoff) + " " + str(angle_cutoff) + " [atomselect top \"(water within 12 of resname " + str(ligand) + ") or (not carbon and not sulfur and protein and chain " + str(chainId) + " and within 12 of resname " + str(ligand) + ") or (not carbon and not sulfur and resname " + str(ligand) + ") and (not lipid)\" frame " + str(frame_index) + "]")
	return x 


def genProteinLigandHBondsAtomSelection(frame_index, dist_cutoff, angle_cutoff, chainId, ligand_chain, solvent_id):
	"""
		Query for hydrogen bonds involving water near the peptide ligand and GPCR
	"""
	x = evaltcl("measure hbonds " + str(dist_cutoff) + " " + str(angle_cutoff) + " [atomselect top \"(water within 12 of (protein within 12 of chain " + str(ligand_chain) + ")) or (not carbon and not sulfur and protein and chain " + str(chainId) + " and within 12 of chain " + str(ligand_chain) + ") or (not carbon and not sulfur and chain " + str(ligand_chain) + ") and (not lipid)\" frame " + str(frame_index) + "]")
	return x 


def calcDonorAcceptors(frame_index, dist_cutoff, angle_cutoff, chainId=None, ligand=None):
	"""
		Calculate donor and acceptor atom pairs for hydrogen bonds in terms of numeric indices 
	"""

	print("calcDonorAcceptors() " + str(ligand))

	### For all cases except where you have peptide ligand
	if(ligand != None and solvent_id != "IP3" and "CHAIN" in ligand):
		ligand_chain = ligand.strip().split("_")[1]
		x = genProteinLigandHBondsAtomSelection(frame_index, dist_cutoff, angle_cutoff, chainId, ligand_chain, solvent_id)
	
	### Retrieve hydrogen bond pair indices for peptide ligand
	else:
		x = genHBondsAtomSelection(frame_index, dist_cutoff, angle_cutoff, chainId, ligand, solvent_id)
	
	### Parse atom indices 

	lists = x.split("}")
	list1 = lists[0].split("{")[1].split(" ")
	donors = [int(a) for a in list1]
	list2 = lists[1].split("{")[1].split(" ")
	acceptors = [int(a) for a in list2]
	return donors, acceptors


def calcHBondFrame(f, frame_index, dist_cutoff, angle_cutoff, lookup_table, chainId=None, ligand=None, crystal=False):
	"""
		Calculate the hydrogen bonds for a single frame of simulation OR for crystal structure. 
	"""

	crystalHBond = []
	if(frame_index %20 == 0): print("calcHBondFrame() for frame_index: " + str(frame_index))

	### Retrieve donor and acceptor pairs based on chainId and ligand query
	if(ligand == None):
		donors, acceptors = calcDonorAcceptors(frame_index, dist_cutoff, angle_cutoff, chainId) #indices 
	else:
		donors, acceptors = calcDonorAcceptors(frame_index, dist_cutoff, angle_cutoff, chainId, ligand) #indices 
	donors, acceptors = removeDuplicateIndices(donors, acceptors) # Remove duplicates 
	donor_label = indexToLabel(frame_index, donors, crystal)
	acceptor_label = indexToLabel(frame_index, acceptors, crystal)


	### Write file header 
	if(crystal == False):
		if(ligand == None):
			f.write("\nHydrogen Bond Water Frame: " + str(frame_index-1) + "\n")
		else:
			f.write("\nLigand Hydrogen Bond Water Frame: " + str(frame_index-1) + "\n")

	### Get rid of chain information if no peptide ligand 
	if(ligand == None or "CHAIN" not in ligand):
		donor_label = [dl.split("_")[0] for dl in donor_label]
		acceptor_label = [al.split("_")[0] for al in acceptor_label]

	### For every donor acceptor pair
	for index, v in enumerate(donor_label):
		print("Index: " + str(index) + " " + donor_label[index] + " " + acceptor_label[index])
		donor_atom = donor_label[index].replace(solvent_id, "HOH")
		acceptor_atom = acceptor_label[index].replace(solvent_id, "HOH")
		if(ligand != None):

			### Normal small ligand drug molecule
			if("CHAIN" not in ligand):
				if(ligand in donor_atom):
					donor_atom = "LIG-" + donor_atom.split("-")[1]
				if(ligand in acceptor_atom): 
					acceptor_atom = "LIG-" + acceptor_atom.split("-")[1]
			
			### Peptide ligand case: ONLY PLACE WHERE CHAIN INFO STILL EXISTS AND IS USED
			else: 
				ligand_chain = ligand.strip().split("_")[1]
				donor_atom, donor_chain = donor_atom.strip().split("_")
				acceptor_atom, acceptor_chain = acceptor_atom.strip().split("_")
				if(donor_chain == ligand_chain):
					donor_atom = "LIG-" + donor_atom.split("-")[1]
				if(acceptor_chain == ligand_chain):
					acceptor_atom = "LIG-" + acceptor_atom.split("-")[1]

		### Append to crystal hydrogen bond output
		if(crystal == True):
			crystalHBond.append((donor_atom, acceptor_atom))

		### Simulation hydrogen bond, write to file 
		else:
			hbond_label = donor_atom + " -- " + acceptor_atom
			f.write(hbond_label + " -- " + str(donors[index]) + " -- " + str(acceptors[index]) + "\n") #include index
	
	### Return crystal hydrogen bond output 
	if(crystal == True): return crystalHBond


def calcVMDHydrogenBondWaterResults(f, stride_val, TOP_PATH, TRAJ_PATH, solventId, chainId=None, ligand=None, crystal=False, dist_cutoff=3.5, angle_cutoff=70):
	"""
		Calculate and write all hydrogen bonds for entire trajectory 
	"""
	print("blah", solventId)
	global solvent_id
	global stride
	if(solventId != None and solventId != 'None'):
		solvent_id = solventId
	stride = stride_val 
	tic = time.clock()
	load_traj(TOP_PATH, TRAJ_PATH)
	numFrames = int(evaltcl('molinfo top get numframes'))
	if(crystal == True):
		crystalHBond = list(set(calcHBondFrame(f, 0, dist_cutoff, angle_cutoff, lookup_table, chainId, ligand, crystal)))
		return crystalHBond
	else:
		for frame_index in range(1, numFrames):
			calcHBondFrame(f, frame_index, dist_cutoff, angle_cutoff, lookup_table, chainId, ligand)
		toc = time.clock()
		computingTime = toc - tic 
		print("calcVMDHydrogenBondWaterResults() Computing Time: " + str(computingTime))
		print(len(lookup_table))
		return numFrames - 1, computingTime


# Interactive module testing 
if __name__ == "__main__":
	if(len(sys.argv) < K_MIN_ARG):
		print(USAGE_STR)
		exit(1)
	TRAJ_PATH = sys.argv[1]
	TOP_PATH = sys.argv[2]
	OUTPUT_FILE = sys.argv[3]
	chainId = None 
	stride = 1
	if("-stride" in sys.argv):
		stride = int(sys.argv[sys.argv.index("-stride")+1])
	if("-chain" in sys.argv):
		chainId = sys.argv[sys.argv.index("-chain") + 1]
	if("-solv" in sys.argv):
		solvent_id = sys.argv[sys.argv.index("-solv") + 1]
	f = createFileWriter(OUTPUT_FILE)
	calcVMDHydrogenBondWaterResults(f, stride, TOP_PATH, TRAJ_PATH, solvent_id, chainId)


