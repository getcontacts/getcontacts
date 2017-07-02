# Author: Anthony Kai Kwang Ma
# Email: anthonyma27@gmail.com, akma327@stanford.edu
# Date: 07/02/17
# vanderwaal.py

### Molecular Dynamics Trajectory Simulation - Van Der Waals Interaction Detection ###

from __future__ import print_function, division
import math
from itertools import product
import time
import numpy as np
import mdtraj as md
from mdtraj.utils import ensure_type
from mdtraj.geometry import compute_distances, compute_angles
from mdtraj.geometry import _geometry
from contact_utils import *

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
	"""
		Fill in data for missing time points 
	"""
	missingTimes = set(timePoints) - set(filtered_candidate_pairs.keys())
	for time in missingTimes:
		filtered_candidate_pairs[time] = []
	return filtered_candidate_pairs

def VDWCutoff(atom1, atom2):
	"""
		Compute the van-der-waals distance between two atoms
	"""
	return ATOM_RADIUS[str(atom1.element)] + ATOM_RADIUS[str(atom2.element)] + VDW_EPSILON


def filter_alpha_carbon(traj, chain_top):
	"""
		Input: Topology for a single chain
		Output: Compute distances between all alpha carbons in the protein.
		Filter out the pairs that are too far apart and return a dictionary 
		mapping time points to a list of tuples representing pairs of alpha
		carbons that are sufficiently close together. 
	"""
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



def fulfillVDWCriterion(atom1, atom2, traj, time):
	"""
		Determine whether two atoms that are not covalently bonded fulfill the 
		van der waals distance criterion.
	"""
	if(atom1.name == 'CA' and atom2.name == 'CA'): # don't count covalent interactions
		return False
	vdwCutoff = VDWCutoff(atom1, atom2)
	p1 = traj.xyz[time, atom1.index, :]
	p2 = traj.xyz[time, atom2.index, :]
	distance = distBetweenTwoPoints(p1, p2)
	if(distance < vdwCutoff):
		return True
	return False


def calcFramePairs(traj, topology, chain_index, time, ca_pair_list):
	"""
		Iterate through all pairs of alpha carbon, and apply a distance filter to identify 
		van der waal interactions between non-alpha carbon atoms within the entire residues.
	"""
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

def calcVDWFramePairs(traj):
	"""
		Calculate all Van-der-waals interactions in protein throughout all frames of simulation
	"""
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

def calcVanderwaalsResults(traj, f, PROTEIN_CODE):
	"""
		Driver for computing Van Der Waals interactions
	"""
	print("\n\nVan Der Waals:" + PROTEIN_CODE + "\n")
	tic = time.clock()
	vdwFramePairs = calcVDWFramePairs(traj)
	toc = time.clock()
	computingTime = toc - tic
	f.write("nFrames:" + str(len(vdwFramePairs)) + "\n")
	f.write("\n\nVan Der Waals:" + PROTEIN_CODE + "\n")
	for index, vanderwaalPairs in enumerate(vdwFramePairs):
		f.write("Vanderwaal Frame: " + str(index) + "\n")
		for pair in vanderwaalPairs:
			atom1, atom2 = pair[0], pair[1]
			f.write(str(atom1[0]) + "_" + str(atom1[1]) +  " -- " + str(atom2[0]) + "_" + str(atom2[1]) + "\n")
	f.write("\nComputing Time:" + str(computingTime) + "\n")
	print("\nComputing Time:" + str(computingTime) + "\n")
	return computingTime

