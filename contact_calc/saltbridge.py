# Author: Anthony Kai Kwang Ma
# Email: anthonyma27@gmail.com, akma327@stanford.edu
# Date: 07/02/17
# saltbridge.py 

### Molecular Dynamics Trajectory Simulation - Salt Bridge Detection ###

from __future__ import print_function, division
import itertools
from itertools import product
import math
import numpy as np
import time
import mdtraj as md
from mdtraj.utils import ensure_type
from mdtraj.geometry import compute_distances, compute_angles
from mdtraj.geometry import _geometry
from contact_utils import *

ANIONIC_LIST = ['ASP', 'GLU']
CATIONIC_LIST = ['LYS', 'ARG', 'HIS']

ASP_FILTER = ['OD'] 
GLU_FILTER = ['OE'] 
LYS_FILTER = ['NZ']
ARG_FILTER = ['NH']
HIS_FILTER = ['NE', 'ND']

sel_basic = "(resname ARG or resname LYS or resname HIS) and (name NH* or name NZ or name NE or name ND)"
sel_acidic = "(resname ASP or resname GLU) and (name OE* or name OD*)"

CUTOFF_DISTANCE = .40 # 4 A as defined by the paper


def initSaltBridgeChainDict(initFrame):
	"""
		Determine list of anions and cations in each chain of the protein
	"""
	topology = initFrame.topology
	chainDict = {}
	for chain_index in range(topology.n_chains):
		O_List = []
		N_List = []
		for index, atom in enumerate(topology.chain(chain_index).atoms):
			name = atom.name
			residueName = atom.residue.name
			if(residueName in ANIONIC_LIST):
				if(residueName == 'ASP'):
					appendAminoAcidToList(O_List, (atom, chain_index), ASP_FILTER)
				if(residueName == 'GLU'):
					appendAminoAcidToList(O_List, (atom, chain_index), GLU_FILTER)
			if(residueName in CATIONIC_LIST):
				if(residueName == 'LYS'):
					appendAminoAcidToList(N_List, (atom, chain_index), LYS_FILTER)
				if(residueName == 'ARG'):
					appendAminoAcidToList(N_List, (atom, chain_index), ARG_FILTER)
				if(residueName == 'HIS'):
					appendAminoAcidToList(N_List, (atom, chain_index), HIS_FILTER)
		chainDict[chain_index] = (O_List, N_List)
	return chainDict


def calcSaltBridgeFramePairs(traj, chainDict):
	"""
		Compute all salt bridges in protein throughout the MD simulation
	"""
	sbFramePairs = {}
	for chain_index in chainDict.keys():
		O_List, N_List = chainDict[chain_index]
		if(len(O_List) == 0 or len(N_List) == 0):
			continue
		pairKeys = []
		atompairs = []
		for bridge_donor in N_List:
			for bridge_acceptor in O_List:
				pairKeys.append([bridge_donor, bridge_acceptor])
				atompairs.append([int(bridge_donor[0].index), int(bridge_acceptor[0].index)])
		atompairs = np.array(atompairs)
		pairDistances = md.compute_distances(traj, atompairs)
		for time in range(len(pairDistances)):
			t_distances = pairDistances[time]
			sb_indices = [i for i in range(len(t_distances)) if t_distances[i] <= CUTOFF_DISTANCE]
			if(time not in sbFramePairs):
				sbFramePairs[time] = [pairKeys[i] for i in sb_indices]
			else:
				sbFramePairs[time] += [pairKeys[i] for i in sb_indices]
	sbFramePairs = dictToList(sbFramePairs)
	return sbFramePairs


def calcSaltBridgeResults(traj, f, PROTEIN_CODE):
	"""
		Driver for computing salt bridges 
	"""
	print("\n\nSalt Bridges:" + PROTEIN_CODE + "\n")
	tic = time.clock()
	chainDict = initSaltBridgeChainDict(traj)
	sbFramePairs = calcSaltBridgeFramePairs(traj, chainDict)
	toc = time.clock()
	computingTime = toc - tic
	f.write("nFrames:" + str(len(sbFramePairs)) + "\n")
	f.write("\n\nSalt Bridges:" + PROTEIN_CODE + "\n")
	for index, sbPairs in enumerate(sbFramePairs):
		f.write("Salt Bridge Frame: " + str(index) + "\n")
		for pairs in sbPairs:
			atom1, atom2 = pairs[0], pairs[1]
			f.write(str(atom1[0]) + "_" + str(atom1[1]) +  " -- " + str(atom2[0]) + "_" + str(atom2[1]) + "\n")
	f.write("\nComputing Time:" + str(computingTime) + "\n")
	print("\nComputing Time:" + str(computingTime) + "\n")
	return computingTime


