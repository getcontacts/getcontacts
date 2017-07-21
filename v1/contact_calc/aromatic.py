# Author: Anthony Kai Kwang Ma
# Email: anthonyma27@gmail.com, akma327@stanford.edu
# Date: 07/02/17
# aromatic.py

### Molecular Dynamics Trajectory Simulation - Pi-stacking and T-stacking Interaction Detection ###

from __future__ import print_function, division
from itertools import product
import math
import numpy as np
import time
import mdtraj as md
from mdtraj.utils import ensure_type
from mdtraj.geometry import compute_distances, compute_angles
from mdtraj.geometry import _geometry

from contact_utils import *

PI_PI_SOFT_CUTOFF_DIST = .8
PI_PI_INTERACTING_DIST_CUTOFF = .7#7 Ang
PI_STACKING_ANGLE_TOLERANCE = 30 #20 #15
PI_PI_PSI_ANGLE_CUTOFF = 45

T_SOFT_CUTOFF_DIST = 0.6
T_STACKING_DIST_CUTOFF = 0.5 #0.4
T_STACKING_ANGLE_TOLERANCE = 30
T_STACK_PSI_ANGLE_CUTOFF = 45

AROMATIC_DONOR = ['PHE', 'TYR', 'TRP']

PHE_FILTER = ['CG', 'CE1', 'CE2']
TRP_FILTER = ['CD2', 'CZ2', 'CZ3']
TYR_FILTER = ['CG', 'CE1', 'CE2']


def fillTimeGaps(timePoints, filtered_candidate_pairs):
	"""
		Fill in the time points that should be empty but didn't make soft cutoff
	"""
	missingTimes = set(timePoints) - set(filtered_candidate_pairs.keys())
	for time in missingTimes:
		filtered_candidate_pairs[time] = {}
	return filtered_candidate_pairs

def calcGamma(xnorm, ynorm, xcent, ycent):
	"""
		Takes in the candidate residues and do the geometry to see 
		if any pairs would lead to pi stacking
	"""
	gamma = angleBtwnVec(xnorm, pointsToVec(xcent, ycent))
	psi = angleBtwnVec(ynorm, pointsToVec(ycent, xcent))
	return min(gamma, psi)

def atomInfo(atom):
	"""
		Return atom properties
	"""
	return (atom.name, atom.residue.name)

def aromaticFilter(atom, chain_index, cand_aromatic_dict):
	"""
		Find residue indices of aromatic residues
	"""
	indicator, residueID = atomInfo(atom)
	if(residueID in AROMATIC_DONOR):
		if(residueID == 'PHE'):
			appendAminoAcidToDict(cand_aromatic_dict, (atom, chain_index), PHE_FILTER)
		if(residueID == 'TRP'):
			appendAminoAcidToDict(cand_aromatic_dict, (atom, chain_index), TRP_FILTER)
		if(residueID == 'TYR'):
			appendAminoAcidToDict(cand_aromatic_dict, (atom, chain_index), TYR_FILTER)


def generatePairKeys(pairKeys, atomPairs, k1, aromatic1, k2, aromatic2):
	"""
		Generate pairs of aromatic interactions 
	"""
	pairKeys.append(((k1, aromatic1[0]), (k2, aromatic2[0])))
	atomPairs.append([int(aromatic1[0][0].index), int(aromatic2[0][0].index)])
	pairKeys.append(((k1, aromatic1[0]), (k2, aromatic2[1])))
	atomPairs.append([int(aromatic1[0][0].index), int(aromatic2[1][0].index)])
	pairKeys.append(((k1, aromatic1[0]), (k2, aromatic2[2])))
	atomPairs.append([int(aromatic1[0][0].index), int(aromatic2[2][0].index)])

	pairKeys.append(((k1, aromatic1[1]), (k2, aromatic2[0])))
	atomPairs.append([int(aromatic1[1][0].index), int(aromatic2[0][0].index)])
	pairKeys.append(((k1, aromatic1[1]), (k2, aromatic2[1])))
	atomPairs.append([int(aromatic1[1][0].index), int(aromatic2[1][0].index)])
	pairKeys.append(((k1, aromatic1[1]), (k2, aromatic2[2])))
	atomPairs.append([int(aromatic1[1][0].index), int(aromatic2[2][0].index)])

	pairKeys.append(((k1, aromatic1[2]), (k2, aromatic2[0])))
	atomPairs.append([int(aromatic1[2][0].index), int(aromatic2[0][0].index)])
	pairKeys.append(((k1, aromatic1[2]), (k2, aromatic2[1])))
	atomPairs.append([int(aromatic1[2][0].index), int(aromatic2[1][0].index)])
	pairKeys.append(((k1, aromatic1[2]), (k2, aromatic2[2])))
	atomPairs.append([int(aromatic1[2][0].index), int(aromatic2[2][0].index)])

def filtered_candidate_pairs_by_chain(traj, cand_aromatic_dict, soft_cutoff_dist):
	"""
		Perform soft distance cutoff
	"""
	filtered_candidate_pairs = {}
	pairKeys = []
	atomPairs = []
	tempDict = {}
	for k1 in cand_aromatic_dict.keys():
		for k2 in cand_aromatic_dict.keys():
			if(k1 != k2 and (k1, k2) not in tempDict and (k2, k1) not in tempDict):
				tempDict[(k1, k2)] = 1
				aromatic1 = cand_aromatic_dict[k1]
				aromatic2 = cand_aromatic_dict[k2]
				generatePairKeys(pairKeys, atomPairs, k1, aromatic1, k2, aromatic2)
	# calculate atom pair distances
	atomPairs = np.array(atomPairs)
	pairDistances = md.compute_distances(traj, atomPairs)
	for time in range(len(pairDistances)):
		t_distances = pairDistances[time]
		for i in range(0, len(t_distances), 9):
			arom_arom_distances = t_distances[i:i+9]
			if(min(arom_arom_distances) <= soft_cutoff_dist):
				arom1_key = pairKeys[i][0][0]
				arom2_key = pairKeys[i][1][0]
				if(time not in filtered_candidate_pairs.keys()):
					filtered_candidate_pairs[time] = {arom1_key:[arom2_key]}
				else:
					if(arom1_key not in filtered_candidate_pairs[time].keys()):
						filtered_candidate_pairs[time][arom1_key] = [arom2_key]
					else:
						filtered_candidate_pairs[time][arom1_key].append(arom2_key)
	filtered_candidate_pairs = fillTimeGaps(range(len(pairDistances)), filtered_candidate_pairs)
	return filtered_candidate_pairs


def calcPiStackingFramePairs(traj, chainDict):
	"""
		Apply distance and angle geometric criterion for pi-stacking interactions 
	"""
	psFramePairs = {}
	for chain_index in chainDict.keys():
		cand_aromatic_dict = chainDict[chain_index]
		if(len(cand_aromatic_dict) == 0):
			continue
		filtered_candidate_pairs = filtered_candidate_pairs_by_chain(traj, cand_aromatic_dict, PI_PI_SOFT_CUTOFF_DIST)
		for time, arom_arom_dict in filtered_candidate_pairs.iteritems():
			frame_pairs = []
			for arom1_key in arom_arom_dict.keys():
				arom2_key_list = arom_arom_dict[arom1_key]
				aromatic1 = cand_aromatic_dict[arom1_key]
				if(len(aromatic1) == 3):
					arom_center1 = calcCentroid(aromatic1, traj, time)
					for arom2_key in arom2_key_list:
						aromatic2 = cand_aromatic_dict[arom2_key]
						if(len(aromatic2) == 3):
							arom_center2 = calcCentroid(aromatic2, traj, time)
							centroidDist = distBetweenTwoPoints(arom_center1, arom_center2)
							if(centroidDist <= PI_PI_INTERACTING_DIST_CUTOFF):
								normVec1, planeCoord1 = calcNormVec(aromatic1, traj, time)
								normVec2, planeCoord2 = calcNormVec(aromatic2, traj, time)
								alignAngle = angleBtwnVec(normVec1, normVec2) #whether norm vec point in same direction
								alignAngle = min(math.fabs(alignAngle - 0), math.fabs(alignAngle - 180))
								if(alignAngle <=PI_STACKING_ANGLE_TOLERANCE):
									psiAngle1 = calcPsiAngle(arom_center1, arom_center2, normVec1)
									psiAngle2 = calcPsiAngle(arom_center2, arom_center1, normVec2)
									psiAngle = min(psiAngle1, psiAngle2)
									if(psiAngle < PI_PI_PSI_ANGLE_CUTOFF):
										pairInfo = [aromatic1, aromatic2, centroidDist, alignAngle, psiAngle]
										frame_pairs.append(pairInfo)
			if(time not in psFramePairs):
				psFramePairs[time] = frame_pairs
			else:
				psFramePairs[time] += frame_pairs
	psFramePairs = dictToList(psFramePairs)
	return psFramePairs


def initFaceFaceAromaticChainDict(initFrame):
	"""
		Determine candidate aromatic pairs in each chain of protein
	"""
	topology = initFrame.topology
	chainDict = {}
	for chain_index in range(topology.n_chains):
		cand_aromatic_dict = {}
		for index, atom in enumerate(topology.chain(chain_index).atoms):
			aromaticFilter(atom, chain_index, cand_aromatic_dict)
		for key in cand_aromatic_dict.keys():
			cand_aromatic_dict[key] = list(cand_aromatic_dict[key])
		chainDict[chain_index] = cand_aromatic_dict
	return chainDict


def calcTStackingFramePairs(traj, chainDict):
	"""
		Apply distance and angle geometric criterion for T-stacking interactions 
	"""
	tsFramePairs = {}
	for chain_index in chainDict.keys():
		cand_aromatic_dict = chainDict[chain_index]
		if(len(cand_aromatic_dict) == 0):
			continue
		filtered_candidate_pairs = filtered_candidate_pairs_by_chain(traj, cand_aromatic_dict, T_SOFT_CUTOFF_DIST)
		for time, arom_arom_dict in filtered_candidate_pairs.iteritems():
			frame_pairs = []
			for arom1_key in arom_arom_dict.keys():
				arom2_key_list = arom_arom_dict[arom1_key]
				aromatic1 = cand_aromatic_dict[arom1_key]
				if(len(aromatic1) == 3):
					arom_center1 = calcCentroid(aromatic1, traj, time)
					for arom2_key in arom2_key_list:
						aromatic2 = cand_aromatic_dict[arom2_key]
						if(len(aromatic2) == 3):
							arom_center2 = calcCentroid(aromatic2, traj, time)
							centroidDist = distBetweenTwoPoints(arom_center1, arom_center2)
							if(centroidDist <= T_STACKING_DIST_CUTOFF):
								normVec1, planeCoord1 = calcNormVec(aromatic1, traj, time)
								normVec2, planeCoord2 = calcNormVec(aromatic2, traj, time)
								perpAngleDeviation = math.fabs(angleBtwnVec(normVec1, normVec2) - 90)
								if(perpAngleDeviation < T_STACKING_ANGLE_TOLERANCE):
									psiAngle1 = calcPsiAngle(arom_center1, arom_center2, normVec1)
									psiAngle2 = calcPsiAngle(arom_center2, arom_center1, normVec2)
									psiAngle = min(psiAngle1, psiAngle2)
									if(psiAngle < T_STACK_PSI_ANGLE_CUTOFF):
										pairInfo = [aromatic1, aromatic2, centroidDist, perpAngleDeviation, psiAngle]
										frame_pairs.append(pairInfo)
			if(time not in tsFramePairs):
				tsFramePairs[time] = frame_pairs
			else:
				tsFramePairs[time] += frame_pairs
	tsFramePairs = dictToList(tsFramePairs)
	return tsFramePairs


def initFaceEdgeAromaticChainDict(initFrame):
	"""
		Determine candidate aromatic pairs in each chain of protein
	"""
	topology = initFrame.topology
	chainDict = {}
	for chain_index in range(topology.n_chains):
		cand_aromatic_dict = {}
		for index, atom in enumerate(topology.chain(chain_index).atoms):
			aromaticFilter(atom, chain_index, cand_aromatic_dict)
		for key in cand_aromatic_dict.keys():
			cand_aromatic_dict[key] = list(cand_aromatic_dict[key])
		chainDict[chain_index] = cand_aromatic_dict
	return chainDict


def calcPiStackingResults(traj, f, PROTEIN_CODE):
	"""
		Driver for Pi-stacking interactions 
	"""
	print("\n\nPi-Stacking:" + PROTEIN_CODE + "\n")
	nFrames = len(traj)
	tic = time.clock()
	chainDict = initFaceFaceAromaticChainDict(traj)
	psFramePairs = calcPiStackingFramePairs(traj, chainDict)
	toc = time.clock()
	computingTime = toc - tic 
	f.write("nFrames:" + str(len(psFramePairs)) + "\n")
	f.write("\n\nPi-Stacking:" + PROTEIN_CODE + "\n")
	for index, pi_stackPairs in enumerate(psFramePairs):
		f.write("Pi_Stacking Frame: " + str(index) + "\n")
		for pair in pi_stackPairs:
			aromatic1, aromatic2 = pair[0], pair[1]
			aromatic_str1 = str(aromatic1[0][0]).split("-")[0] + "_" + str(aromatic1[0][1])
			aromatic_str2 = str(aromatic2[0][0]).split("-")[0] + "_" + str(aromatic2[0][1])
			f.write(aromatic_str1 + " -- " + aromatic_str2 + "\n")
	f.write("\nComputing Time:" + str(computingTime) + "\n")
	print("\nComputing Time:" + str(computingTime) + "\n")
	return computingTime


def calcTStackingResults(traj, f, PROTEIN_CODE):
	"""
		Driver for T-stacking interactions 
	"""
	print("\n\nT-Stacking:" + PROTEIN_CODE + "\n")
	nFrames = len(traj)
	tic = time.clock()
	tsFramePairs = []
	chainDict = initFaceEdgeAromaticChainDict(traj)
	tsFramePairs = calcTStackingFramePairs(traj, chainDict)
	toc = time.clock()
	computingTime = toc - tic
	f.write("nFrames:" + str(len(tsFramePairs)) + "\n")
	f.write("\n\nT-Stacking:" + PROTEIN_CODE + "\n")
	for index, t_stackPairs in enumerate(tsFramePairs):
		f.write("T_Stacking Frame: " + str(index) + "\n")
		for pair in t_stackPairs:
			aromatic1, aromatic2 = pair[0], pair[1]
			aromatic_str1 = str(aromatic1[0][0]).split("-")[0] + "_" + str(aromatic1[0][1])
			aromatic_str2 = str(aromatic2[0][0]).split("-")[0] + "_" + str(aromatic2[0][1])
			f.write(aromatic_str1 + " -- " + aromatic_str2 + "\n")
	f.write("\nComputing Time:" + str(computingTime) + "\n")
	print("\nComputing Time:" + str(computingTime) + "\n")
	return computingTime

