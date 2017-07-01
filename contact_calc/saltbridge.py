# Author: Anthony Kai Kwang Ma
# Email: akma327@stanford.edu
# Date: September 29, 2015

# Molecular Dynamics Trajectory Simulation - Salt Bridge Detection

from __future__ import print_function, division
import math
from itertools import product
import itertools
import numpy as np
import mdtraj as md
from mdtraj.utils import ensure_type
from mdtraj.geometry import compute_distances, compute_angles
from mdtraj.geometry import _geometry
from noncovalentInteractionUtils import *

__all__ = ['initSaltBridgeChainDict', 'calcSaltBridgeFramePairs']

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

def calcSaltBridgeFramePairs(traj, chainDict):
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


def initSaltBridgeChainDict(initFrame):
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
			#Handle the nitrogens
			if(residueName in CATIONIC_LIST):
				if(residueName == 'LYS'):
					appendAminoAcidToList(N_List, (atom, chain_index), LYS_FILTER)
				if(residueName == 'ARG'):
					appendAminoAcidToList(N_List, (atom, chain_index), ARG_FILTER)
				if(residueName == 'HIS'):
					appendAminoAcidToList(N_List, (atom, chain_index), HIS_FILTER)
		chainDict[chain_index] = (O_List, N_List)
	return chainDict



