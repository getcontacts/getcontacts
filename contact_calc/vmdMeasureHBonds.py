# Author: Anthony Kai Kwang Ma
# Email: anthonyma27@gmail.com, akma327@stanford.edu
# Date: 07/02/17
# vmdMeasureHBonds.py

### Molecular Dynamics Trajectory Simulation - HBond Detection ###

from __future__ import division

import vmd, molecule
from vmd import *
import sys
import os
import subprocess
import numpy as np
import time 

### Global Variables
start=0
stop=-1
stride=1
smoothing=0
trajid = None
solvent_id = "TIP3"
lookup_table = {}


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


def indexToLabel(frame_index, hbond_index_list, crystal=False):
	"""
		Convert list of atom indices for a frame to atom label that includes information 
		regarding resname, resid, atom_name, and chain in the format of 
		<resname><resid>-<atom_name>_<chain> (ie GLU142-NZ_A)
	"""
	global lookup_table
	label_list = []
	for display_index, atom_index in enumerate(hbond_index_list):
		if(atom_index in lookup_table):
			val = lookup_table[atom_index]
			label_list.append(val)
		else:
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

	if(ligand == None):
		if(solvent_id == "IP3"):
			if(chainId == None):
				x = evaltcl("measure hbonds " + str(dist_cutoff) + " " + str(angle_cutoff) + " [atomselect top \"(resname IP3 and within 3.5 of protein) or protein and not lipid and not carbon\" frame " + str(frame_index) + "]")
			else:
				x = evaltcl("measure hbonds " + str(dist_cutoff) + " " + str(angle_cutoff) + " [atomselect top \"(resname IP3 and within 3.5 of protein and chain " + str(chainId) + ") or protein and chain " + str(chainId) + " and not lipid and not carbon\" frame " + str(frame_index) + "]")
		else:
			if(chainId == None):
				x = evaltcl("measure hbonds " + str(dist_cutoff) + " " + str(angle_cutoff) + " [atomselect top \"(water and within 3.5 of protein) or protein and not lipid and not carbon\" frame " + str(frame_index) + "]")
			else:
				x = evaltcl("measure hbonds " + str(dist_cutoff) + " " + str(angle_cutoff) + " [atomselect top \"(water and within 3.5 of protein and chain " + str(chainId) + ") or protein and chain " + str(chainId) + " and not lipid and not carbon\" frame " + str(frame_index) + "]")	
	else:
		if(solvent_id == "IP3"):
			if(chainId == None):
				x = evaltcl("measure hbonds " + str(dist_cutoff) + " " + str(angle_cutoff) + " [atomselect top \"(resname IP3 and within 12 of resname " + str(ligand) + ") or (not carbon and not sulfur and protein within 12 of resname " + str(ligand) + ") or (not carbon and not sulfur and resname " + str(ligand) + ") and (not lipid)\" frame " + str(frame_index) + "]")
			else:
				x = evaltcl("measure hbonds " + str(dist_cutoff) + " " + str(angle_cutoff) + " [atomselect top \"(resname IP3 and within 12 of resname " + str(ligand) + ") or (not carbon and not sulfur and protein and chain " + str(chainId) + " and within 12 of resname " + str(ligand) + ") or (not carbon and not sulfur and resname " + str(ligand) + ") and (not lipid)\" frame " + str(frame_index) + "]")
		else:
			if(chainId == None):
				x = evaltcl("measure hbonds " + str(dist_cutoff) + " " + str(angle_cutoff) + " [atomselect top \"(water within 12 of resname " + str(ligand) + ") or (not carbon and not sulfur and protein within 12 of resname " + str(ligand) + ") or (not carbon and not sulfur and resname " + str(ligand) + ") and (not lipid)\" frame " + str(frame_index) + "]")
			else:
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
		donor_atom = donor_label[index].replace(solvent_id, "HOH")
		acceptor_atom = acceptor_label[index].replace(solvent_id, "HOH")
		if(ligand != None):

			### Normal small ligand drug molecule
			if("CHAIN" not in ligand):
				if(ligand in donor_atom):
					donor_atom = "LIG-" + donor_atom.split("-")[1]
				if(ligand in acceptor_atom): 
					acceptor_atom = "LIG-" + acceptor_atom.split("-")[1]
			
			### Peptide ligand case:
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
		return numFrames - 1, computingTime
