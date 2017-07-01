# Author: Anthony Kai Kwang Ma
# Email: akma327@stanford.edu
# Date: December 14, 2015

# Molecular Dynamics Trajectory Simulation - Van Der Waals Interaction Detection

from __future__ import print_function, division
import math
from itertools import product
import time
import numpy as np
import mdtraj as md
from mdtraj.utils import ensure_type
from mdtraj.geometry import compute_distances, compute_angles
from mdtraj.geometry import _geometry
from noncovalentInteractionUtils import *

__all__ = ['calcVDWFramePairs']

ALPHA_CARBON_DIST_CUTOFF = 1 #10 angstrom
VDW_EPSILON = .025 #.25 Angstrom

AMINO_ACIDS = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE','LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

#http://crystalmaker.com/support/tutorials/crystalmaker/atomic-radii/index.html
ATOM_RADIUS = {'hydrogen': .120, 
				'carbon': .170,
				'nitrogen': .155, 
				'oxygen': .152,
				'sulfur': .180}

def fillTimeGaps(timePoints, filtered_candidate_pairs):
	missingTimes = set(timePoints) - set(filtered_candidate_pairs.keys())
	for time in missingTimes:
		filtered_candidate_pairs[time] = []
	return filtered_candidate_pairs

#given two atoms determine the vanderwaal cutoff as radius1 + radius2 + epsilon
def VDWCutoff(atom1, atom2):
	return ATOM_RADIUS[str(atom1.element)] + ATOM_RADIUS[str(atom2.element)] + VDW_EPSILON


# Input: the topology for a single chain
# Output: Utilize the md.compute_distance to find distance between all alpha carbons, filter out
# the pairs that are too far apart and return a dictionary with key = time, and value being a list
# of tuples of alpha carbon atoms that are sufficiently close together 
def filter_alpha_carbon(traj, chain_top):
	filtered_candidate_pairs = {} # outer dictionary key = time, value = list of tuples (ca_1, ca_2)
	ca_keys, ca_indices = [], []
	for ca_atom in chain_top.atoms_by_name('CA'):
		ca_keys.append(ca_atom)
		ca_indices.append(ca_atom.index)
	ca_pair_keys, ca_pair_indices = [], []
	n = len(ca_indices)
	if(n == 0): return 
	for i1 in range(n):
		for i2 in range(i1 + 1, n):
			ca_pair_keys.append((ca_keys[i1], ca_keys[i2]))
			ca_pair_indices.append([ca_indices[i1], ca_indices[i2]])
	ca_pair_indices = np.array(ca_pair_indices)
	if(len(ca_pair_indices) == 0): return 
	pairDistances = md.compute_distances(traj, ca_pair_indices)
	for time in range(len(pairDistances)):
		t_distances = pairDistances[time]
		potential_indices = [i for i in range(len(t_distances)) if t_distances[i] <= ALPHA_CARBON_DIST_CUTOFF]
		filtered_candidate_pairs[time] = [ca_pair_keys[i] for i in potential_indices]
	filtered_candidate_pairs = fillTimeGaps(range(len(pairDistances)), filtered_candidate_pairs)
	return filtered_candidate_pairs


# Verify whether two atoms fulfill van der waals criterion 
# Can't have vdw between two bonded CA atoms because that is a covalent bond
# Must be within vanderwaal radii 
# Using the brute force method 
def fulfillVDWCriterion(atom1, atom2, traj, time):
	if(atom1.name == 'CA' and atom2.name == 'CA'): # don't count covalent interactions
		return False
	vdwCutoff = VDWCutoff(atom1, atom2)
	p1 = traj.xyz[time, atom1.index, :]
	p2 = traj.xyz[time, atom2.index, :]
	distance = distBetweenTwoPoints(p1, p2)
	if(distance < vdwCutoff):
		return True
	return False

# Iterate through all pairs of alpha carbon, and then do another distance filter for vdw
# interaction between the individual non alpha carbon atoms within the entire residue 
def calcFramePairs(traj, topology, chain_index, time, ca_pair_list):
	frame_pairs = []
	for ca_1, ca_2 in ca_pair_list:
		resid1_atoms = ca_1.residue.atoms
		resid2_atoms = ca_2.residue.atoms
		for atom1 in resid1_atoms:
			for atom2 in resid2_atoms:
				if fulfillVDWCriterion(atom1, atom2, traj, time):
					pairKey = [[(atom1, chain_index), (atom2, chain_index)]]
					frame_pairs += pairKey
	return frame_pairs

# Calculate all the van der waal pairs for all time points in trajectory 
def calcVDWFramePairs(traj):
	vdwFramePairs = {}
	topology = traj.topology
	for chain_index in range(topology.n_chains):
		filtered_CA = filter_alpha_carbon(traj, topology.chain(chain_index))
		if(filtered_CA == None):
			continue
		for time, ca_pair_list in filtered_CA.iteritems():
			frame_pairs = calcFramePairs(traj, topology, chain_index, time, ca_pair_list)
			if(time not in vdwFramePairs):
				vdwFramePairs[time] = frame_pairs
			else:
				vdwFramePairs[time] += frame_pairs
	vdwFramePairs = dictToList(vdwFramePairs)
	return vdwFramePairs




########################### Other Code #############################

# Verify whether two atoms fulfill van der waals criterion 
# Can't have vdw between two bonded CA atoms because that is a covalent bond
# Must be within vanderwaal radii 
def fulfillVDWCriterion2(pairKey, distance):
	atom1, atom2 = pairKey[0][0], pairKey[1][0]
	if(atom1.name == 'CA' and atom2.name == 'CA'): # don't count covalent interactions
		return False
	vdwCutoff = VDWCutoff(atom1, atom2)
	if(distance < vdwCutoff):
		return True
	return False

# Iterate through all pairs of alpha carbon, and then do another distance filter for vdw
# interaction between the individual non alpha carbon atoms within the entire residue 
def calcFramePairs2(traj, topology, chain_index, time, ca_pair_list):
	frame_pairs = []
	for ca_1, ca_2 in ca_pair_list: # might be faster to do the dumb way
		resid1 = ca_1.residue.index
		resid2 = ca_2.residue.index
		resid1_atoms = topology.chain(chain_index).residue(resid1).atoms
		resid2_atoms = topology.chain(chain_index).residue(resid2).atoms
		pairKeys, atomPairs = [], []
		for atom1 in resid1_atoms:
			for atom2 in resid2_atoms:
				pairKeys.append([(atom1, chain_index), (atom2, chain_index)])
				atomPairs.append([atom1.index, atom2.index])
		atomPairs = np.array(atomPairs)
		pairDistances = md.compute_distances(traj, atomPairs) # will compute for all time points
		t_distances = pairDistances[time] #only care about this specific time point
		vdw_indices = [i for i in range(len(t_distances)) if fulfillVDWCriterion(pairKeys[i], t_distances[i]) == True]
		frame_pairs += [pairKeys[i] for i in vdw_indices]
	return frame_pairs

