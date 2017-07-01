# Molecular Dynamics Trajectory Simulation - Hydrophobic Interaction Detection
# 10/01/15 - Anthony Ma 

from __future__ import print_function, division
import math
from itertools import product
import numpy as np
from mdtraj.utils import ensure_type
from mdtraj.geometry import compute_distances, compute_angles
from mdtraj.geometry import _geometry
from noncovalentInteractionUtils import *

__all__ = ['hydrophobe_detect']

CUTOFF_DISTANCE = .45 #A as defined by the paper 
#Filtering for significant hydrophobic atoms in side chain
GLY_FILTER = ['CA', 'HA2', 'HA3']
ALA_FILTER = ['CB']
VAL_FILTER = ['CB', 'CG1', 'CG2']
LEU_FILTER = ['CB', 'CG', 'CD1', 'CD2']
ILE_FILTER = ['CB', 'CG1', 'CG2', 'CD1']
MET_FILTER = ['CB', 'CG', 'SD', 'CE']
PHE_FILTER = ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']
TYR_FILTER = ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']
TRP_FILTER = ['CB', 'CG', 'CD1', 'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3']

polarAAList = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TYR', 'TRP']

# GLY_FILTER = ['CA']
# ALA_FILTER = ['CB']
# VAL_FILTER = ['CG1']
# LEU_FILTER = ['CD1']
# ILE_FILTER = ['CD1']
# MET_FILTER = ['SD', 'CE']
# PHE_FILTER = ['CZ']
# TYR_FILTER = ['CZ']
# TRP_FILTER = ['CZ3']


def candidatesToPair(trajectoryFrame, chainDict, hydrophobePairs):
	for chain_index in chainDict.keys():
		cand_Dict = chainDict[chain_index]
		for k1 in cand_Dict.keys():
			for k2 in cand_Dict.keys():
				if(k1 != k2):
					for atom1 in cand_Dict[k1]:
						for atom2 in cand_Dict[k2]:
							atomIndex1 = atom1[0].index
							atomIndex2 = atom2[0].index
							pairVec = tuple(trajectoryFrame.xyz[0, atomIndex1,:] - trajectoryFrame.xyz[0, atomIndex2,:])
							pairDist = length(pairVec)
							if(pairDist <=CUTOFF_DISTANCE):
								pairInfo = [atom1, atom2, pairDist]
								hydrophobePairs.append(pairInfo)

def hydrophobicFilter(atom, chain_index, cand_Dict):
	correspondingResidue = atom.residue.name
	if(correspondingResidue in polarAAList):
		if(correspondingResidue == 'GLY'):
			appendAminoAcidToDict(cand_Dict, (atom, chain_index), GLY_FILTER)
		if(correspondingResidue == 'ALA'):
			appendAminoAcidToDict(cand_Dict, (atom, chain_index), ALA_FILTER)
		if(correspondingResidue == 'VAL'):
			appendAminoAcidToDict(cand_Dict, (atom, chain_index), VAL_FILTER)
		if(correspondingResidue == 'LEU'):
			appendAminoAcidToDict(cand_Dict, (atom, chain_index), LEU_FILTER)
		if(correspondingResidue == 'ILE'):
			appendAminoAcidToDict(cand_Dict, (atom, chain_index), ILE_FILTER)
		if(correspondingResidue == 'MET'):
			appendAminoAcidToDict(cand_Dict, (atom, chain_index), MET_FILTER)
		if(correspondingResidue == 'PHE'):
			appendAminoAcidToDict(cand_Dict, (atom, chain_index), PHE_FILTER)
		if(correspondingResidue == 'TYR'):
			appendAminoAcidToDict(cand_Dict, (atom, chain_index), TYR_FILTER)
		if(correspondingResidue == 'TRP'):
			appendAminoAcidToDict(cand_Dict, (atom, chain_index), TRP_FILTER)

def hydrophobe_detect(trajectoryFrame):
	hydrophobePairs = []
	topology = trajectoryFrame.topology
	chainDict = {} #key: chainID, value: cand_Dict
	for chain_index in range(topology.n_chains):
		cand_Dict = {} #dictionary containing list of the atoms in a particular residue that has hydrophobe properties
		for index, atom in enumerate(topology.chain(chain_index).atoms):
			hydrophobicFilter(atom, chain_index, cand_Dict)
		for key in cand_Dict.keys():
			cand_Dict[key] = list(cand_Dict[key])
		chainDict[chain_index] = cand_Dict
	candidatesToPair(trajectoryFrame, chainDict, hydrophobePairs)
	return hydrophobePairs
