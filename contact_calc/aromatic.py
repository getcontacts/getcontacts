# Author: Anthony Kai Kwang Ma
# Email: akma327@stanford.edu
# Date: October 6, 2015

# Molecular Dynamics Trajectory Simulation - pi-Cationic Interaction Detection

from __future__ import print_function, division
import math
from itertools import product
import numpy as np
import mdtraj as md
from mdtraj.utils import ensure_type
from mdtraj.geometry import compute_distances, compute_angles
from mdtraj.geometry import _geometry
import math
from noncovalentInteractionUtils import *

__all__ = ['face_face_aromatic_detect', 'face_edge_aromatic_detect', 'initFaceFaceAromaticChainDict', 'initFaceEdgeAromaticChainDict', 'calcPiStackingFramePairs', 'calcTStackingFramePairs']

PI_PI_SOFT_CUTOFF_DIST = .8
PI_PI_INTERACTING_DIST_CUTOFF = .7#7 Ang
PI_STACKING_ANGLE_TOLERANCE = 30 #20 #15
PI_PI_PSI_ANGLE_CUTOFF = 45

T_SOFT_CUTOFF_DIST = 0.6
T_STACKING_DIST_CUTOFF = 0.5 #0.4
T_STACKING_ANGLE_TOLERANCE = 30
T_STACK_PSI_ANGLE_CUTOFF = 45

# PI_PI_CENTER_DIST_TOLERANCE = 0.35 #probably should be more like .25
# T_CENTER_DIST_TOLERANCE = 0.35


AROMATIC_DONOR = ['PHE', 'TYR', 'TRP']

PHE_FILTER = ['CG', 'CE1', 'CE2']
TRP_FILTER = ['CD2', 'CZ2', 'CZ3']
TYR_FILTER = ['CG', 'CE1', 'CE2']

# 12/13/15 Version 2.0 

# Fill in the time points that should be empty but didn't make soft cutoff 
def fillTimeGaps(timePoints, filtered_candidate_pairs):
	missingTimes = set(timePoints) - set(filtered_candidate_pairs.keys())
	for time in missingTimes:
		filtered_candidate_pairs[time] = {}
	return filtered_candidate_pairs

#takes in the candidate residues and do the geometry to see if any pairs would lead to pi stacking
def calcGamma(xnorm, ynorm, xcent, ycent):
	gamma = angleBtwnVec(xnorm, pointsToVec(xcent, ycent))
	psi = angleBtwnVec(ynorm, pointsToVec(ycent, xcent))
	return min(gamma, psi)

def atomInfo(atom):
	return (atom.name, atom.residue.name)

def aromaticFilter(atom, chain_index, cand_aromatic_dict):
	indicator, residueID = atomInfo(atom)
	if(residueID in AROMATIC_DONOR):
		if(residueID == 'PHE'):
			appendAminoAcidToDict(cand_aromatic_dict, (atom, chain_index), PHE_FILTER)
		if(residueID == 'TRP'):
			appendAminoAcidToDict(cand_aromatic_dict, (atom, chain_index), TRP_FILTER)
		if(residueID == 'TYR'):
			appendAminoAcidToDict(cand_aromatic_dict, (atom, chain_index), TYR_FILTER)


def generatePairKeys(pairKeys, atomPairs, k1, aromatic1, k2, aromatic2):
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Edge Face Aromatic Interactions ~~~~~~~~~~~~~~~~~~~~~~~~~~


def calcTStackingFramePairs(traj, chainDict):
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



# Version 1.0
def assignFaceFacePairs(time, traj, chainDict, faceFacePairs):
	for chain_index in chainDict.keys():
		cand_aromatic_dict = chainDict[chain_index]
		tempDict = {}
		for k1 in cand_aromatic_dict.keys():
			for k2 in cand_aromatic_dict.keys():
				if(k1 != k2 and (k1, k2) not in tempDict and (k2,k1) not in tempDict):
					tempDict[(k1, k2)] = 1
					aromatic1 = cand_aromatic_dict[k1]
					aromatic2 = cand_aromatic_dict[k2]
					if(len(aromatic1) == 3 and len(aromatic2) == 3):
						arom_center1 = calcCentroid(aromatic1, traj, time)
						arom_center2 = calcCentroid(aromatic2, traj, time)
						centroidDist = distBetweenTwoPoints(arom_center1, arom_center2)
						if(centroidDist <= PI_PI_INTERACTING_DIST_CUTOFF):
							normVec1, planeCoord1 = calcNormVec(aromatic1, traj, time)
							normVec2, planeCoord2 = calcNormVec(aromatic2, traj, time)
							alignAngle = angleBtwnVec(normVec1, normVec2) #whether norm vec point in same direction
							alignAngle = min(math.fabs(alignAngle - 0), math.fabs(alignAngle - 180))
							if(alignAngle <=PI_STACKING_ANGLE_TOLERANCE):
								#method 5: treat arom_center2 as the cation and compute angle between normvec1
								#and the vector(center2, center1)
								psiAngle1 = calcPsiAngle(arom_center1, arom_center2, normVec1)
								psiAngle2 = calcPsiAngle(arom_center2, arom_center1, normVec2)
								psiAngle = min(psiAngle1, psiAngle2)
								if(psiAngle < PI_PI_PSI_ANGLE_CUTOFF):
									pairInfo = [aromatic1, aromatic2, centroidDist, alignAngle, psiAngle]
									faceFacePairs.append(pairInfo)


#Detects face to face aromatic stacking or interaction of pi bonds
def face_face_aromatic_detect(time, traj, initChainDict):
	faceFacePairs = []
	assignFaceFacePairs(time, traj, initChainDict, faceFacePairs)
	return faceFacePairs


def assignFaceEdgePairs(time, traj, chainDict, faceEdgePairs):
	for chain_index in chainDict.keys():
		cand_aromatic_dict = chainDict[chain_index]
		tempDict = {}
		for k1 in cand_aromatic_dict.keys():
			for k2 in cand_aromatic_dict.keys():
				if (k1 != k2 and (k1, k2) not in tempDict and (k2, k1) not in tempDict):
					tempDict[(k1, k2)] = 1
					aromatic1 = cand_aromatic_dict[k1]
					aromatic2 = cand_aromatic_dict[k2]
					if(len(aromatic1) == 3 and len(aromatic2) == 3):
						arom_center1 = calcCentroid(aromatic1, traj, time)
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
									faceEdgePairs.append(pairInfo)


def face_edge_aromatic_detect(time, traj, initChainDict):
	faceEdgePairs = []
	assignFaceEdgePairs(time, traj, initChainDict, faceEdgePairs)
	return faceEdgePairs

