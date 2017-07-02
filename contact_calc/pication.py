# Author: Anthony Kai Kwang Ma
# Email: anthonyma27@gmail.com, akma327@stanford.edu
# Date: 07/02/17
# pication.py

### Molecular Dynamics Trajectory Simulation - pi-Cationic Interaction Detection ###

from __future__ import print_function, division
import math
import time
from itertools import product
import numpy as np
import mdtraj as md
from mdtraj.utils import ensure_type
from mdtraj.geometry import compute_distances, compute_angles
from mdtraj.geometry import _geometry
from contact_utils import *
import math

__all__ = ['pication_detect', 'initPiCationChainDict', 'calcPiCationFramePairs']

SOFT_CUTOFF_DISTANCE = 0.65 # First pass 6.5 Angstrom distance cutof
CUTOFF_DISTANCE = 0.6 #99% occur within 6 Angstroms
CUTOFF_ANGLE = 60 #45 #angle that cation can deviate from the norm vector given by aromatic ring

CATION_DONOR = ['LYS', 'ARG']
AROMATIC_DONOR = ['PHE', 'TYR', 'TRP']

LYS_FILTER = ['NZ']
ARG_FILTER = ['NH']

PHE_FILTER = ['CG', 'CE1', 'CE2']
TRP_FILTER = ['CD2', 'CZ2', 'CZ3']
TYR_FILTER = ['CG', 'CE1', 'CE2']

def fillTimeGaps(timePoints, filtered_candidate_pairs):
	"""
		Fill in data for missing time points 
	"""
	missingTimes = set(timePoints) - set(filtered_candidate_pairs.keys())
	for time in missingTimes:
		filtered_candidate_pairs[time] = {}
	return filtered_candidate_pairs

def cationFilter(atom, chain_index, cation_list):
	"""
		Identify residue indices that are LYS or ARG cations
	"""
	indicator, residueID = atomInfo(atom)
	if(residueID in CATION_DONOR):
		if(residueID == 'LYS'):
			appendAminoAcidToList(cation_list, (atom, chain_index), LYS_FILTER)
		if(residueID == 'ARG'):
			appendAminoAcidToList(cation_list, (atom, chain_index), ARG_FILTER)


def aromaticFilter(atom, chain_index, cand_aromatic_dict):
	"""
		Identify residue indices that are aromatics PHE, TRP, TYR
	"""
	indicator, residueID = atomInfo(atom)
	if(residueID in AROMATIC_DONOR):
		if(residueID == 'PHE'):
			appendAminoAcidToDict(cand_aromatic_dict, (atom, chain_index), PHE_FILTER)
		if(residueID == 'TRP'):
			appendAminoAcidToDict(cand_aromatic_dict, (atom, chain_index), TRP_FILTER)
		if(residueID == 'TYR'):
			appendAminoAcidToDict(cand_aromatic_dict, (atom, chain_index), TYR_FILTER)


def filter_candidate_pairs_by_chain(traj, cation_list, cand_aromatic_dict):
	"""
		Apply distance and angle geometric criterion to filter out cation-aromatic residues
		that can form pi-cation interactions. 
	"""
	filtered_candidate_pairs = {} # outer dictionary key = time, value = inner dictionary with key = aromatic key, value = list of cations
	pairKeys = []
	atomPairs = []
	# generate all candidate pairs of cation and aromatic atoms
	for k1 in cand_aromatic_dict.keys():
		atoms_in_aromatic_resid = cand_aromatic_dict[k1]
		if(len(atoms_in_aromatic_resid) == 3):
			for cation in cation_list:
				pairKeys.append((cation, (k1, atoms_in_aromatic_resid[0])))
				atomPairs.append([int(cation[0].index), int(atoms_in_aromatic_resid[0][0].index)])
				pairKeys.append((cation, (k1, atoms_in_aromatic_resid[1])))
				atomPairs.append([int(cation[0].index), int(atoms_in_aromatic_resid[1][0].index)])
				pairKeys.append((cation, (k1, atoms_in_aromatic_resid[2])))
				atomPairs.append([int(cation[0].index), int(atoms_in_aromatic_resid[2][0].index)])
	# calculate atom pair distances
	atomPairs = np.array(atomPairs)
	pairDistances = md.compute_distances(traj, atomPairs)

	# Filter to get smaller candidate set 
	for time in range(len(pairDistances)):
		t_distances = pairDistances[time]
		for i in range(0, len(t_distances), 3):
			cation_arom_distances = t_distances[i:i+3]
			if(min(cation_arom_distances) <= SOFT_CUTOFF_DISTANCE):
				cation_identifier = pairKeys[i][0]
				aromatic_key = pairKeys[i][1][0]
				if(time not in filtered_candidate_pairs.keys()):
					filtered_candidate_pairs[time] = {aromatic_key:[cation_identifier]}
				else:
					if(aromatic_key not in filtered_candidate_pairs[time].keys()):
						filtered_candidate_pairs[time][aromatic_key] = [cation_identifier]
					else:
						filtered_candidate_pairs[time][aromatic_key].append(cation_identifier)
	filtered_candidate_pairs = fillTimeGaps(range(len(pairDistances)), filtered_candidate_pairs)
	return filtered_candidate_pairs


def initPiCationChainDict(initFrame):
	"""
		Compute cation list and candidate aromatic residue list for all chains in protein. 
		This provides a starting point to identifying pairs of cationic-aromatic residues 
		that can form pi-cation interactions.
	"""
	topology = initFrame.topology
	chainDict = {}
	for chain_index in range(topology.n_chains):
		cation_list = []
		cand_aromatic_dict = {}
		for index, atom in enumerate(topology.chain(chain_index).atoms):
			cationFilter(atom, chain_index, cation_list)
			aromaticFilter(atom, chain_index, cand_aromatic_dict)
		for key in cand_aromatic_dict.keys():
			cand_aromatic_dict[key] = list(cand_aromatic_dict[key])
		chainDict[chain_index] = (cation_list, cand_aromatic_dict)
	return chainDict


def calcPiCationFramePairs(traj, chainDict):
	"""
		Compute all pi-cation interactions in protein throughout MD simulation 
	"""
	pcFramePairs = {}
	for chain_index in chainDict.keys():
		cation_list, cand_aromatic_dict = chainDict[chain_index]
		if(len(cation_list) == 0 or len(cand_aromatic_dict) == 0):
			continue
		filtered_candidate_pairs = filter_candidate_pairs_by_chain(traj, cation_list, cand_aromatic_dict)
		for time, arom_cation_dict in filtered_candidate_pairs.iteritems():
			frame_pairs = []
			for arom_key in arom_cation_dict.keys():
				cation_list = arom_cation_dict[arom_key]
				atoms_in_aromatic_resid = cand_aromatic_dict[arom_key]
				if(len(atoms_in_aromatic_resid) == 3):
					aromatic_center = calcCentroid(atoms_in_aromatic_resid, traj, time)
					aromatic_norm_vec, planeCoord = calcNormVec(atoms_in_aromatic_resid, traj, time)
					for cation in cation_list:
						cation_loc = tuple(traj.xyz[time, cation[0].index, :])
						cation_aromatic_distance = distBetweenTwoPoints(cation_loc, aromatic_center)
						if(cation_aromatic_distance <= CUTOFF_DISTANCE):
							centerToCationVec = pointsToVec(aromatic_center, cation_loc)
							deviationAngle = angleBtwnVec(aromatic_norm_vec, centerToCationVec)
							deviationAngle = min(math.fabs(deviationAngle - 0), math.fabs(deviationAngle - 180))
							if(deviationAngle < CUTOFF_ANGLE):
								pairInfo = [cation, [atoms_in_aromatic_resid[0], atoms_in_aromatic_resid[1], atoms_in_aromatic_resid[2]], aromatic_center, deviationAngle]
								frame_pairs.append(pairInfo)
			if (time not in pcFramePairs):
				pcFramePairs[time] = frame_pairs
			else:
				pcFramePairs[time] += frame_pairs
	pcFramePairs = dictToList(pcFramePairs)
	return pcFramePairs



def calcPiCationResults(traj, f, PROTEIN_CODE):
	"""
		Driver for computing pi-cation interactions
	"""
	print("\n\nPi-Cation:" + PROTEIN_CODE + "\n")
	tic = time.clock()
	chainDict = initPiCationChainDict(traj)
	pcFramePairs = calcPiCationFramePairs(traj, chainDict)
	toc = time.clock()
	computingTime = toc - tic 
	f.write("nFrames:" + str(len(pcFramePairs)) + "\n")
	f.write("\n\nPi-Cation:" + PROTEIN_CODE + "\n")
	for index, picationPairs in enumerate(pcFramePairs):
		f.write("Pi_Cation Frame: " + str(index) + "\n")
		for pair in picationPairs:
			cation, aromatic = pair[0], pair[1]
			cation_str = str(cation[0]) + "_" + str(cation[1])
			aromatic_str = str(aromatic[0][0]).split("-")[0] + "_" + str(aromatic[0][1])
			f.write(cation_str + " -- " + aromatic_str + "\n")
	f.write("\nComputing Time:" + str(computingTime) + "\n")
	print("\nComputing Time:" + str(computingTime) + "\n")
	return computingTime

