#!/share/PI/rondror/software/miniconda/bin/python

import mdtraj as md
import operator
from collections import defaultdict
import time

#Finds out hbonds and their frames in a given trajectory, with certain boundary coordinates. 
#Analyzing frame by frame reduces the memory used. 

__all__ = ['calcHydrogenBondFramePairs']

MIN_X = -2
MAX_X = 2
MIN_Y = -2
MAX_Y = 2
MIN_Z = -1.3
MAX_Z = 1.3

#determine which waters to consider
def atoms_in_bound(frameCoord, boundIndices, allAtoms):
	for index in allAtoms:
		if MIN_X <= frameCoord[index][0] and frameCoord[index][0] <= MAX_X:
			if MIN_Y <= frameCoord[index][1] and frameCoord[index][1] <= MAX_Y:
				if MIN_Z <= frameCoord[index][2] and frameCoord[index][2] <= MAX_Z:

					boundIndices.append(index)
	print("Number of points in bound indices: " + str(len(boundIndices)))

def atom_label(traj, hbond):

	make_label = lambda hbond : '%s -- %s' % (traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2]))
	return make_label(hbond)

#finds the hydrogen bonds in the trajectory
def findHbonds(traj, boundIndices, results, frameIndex, resultsFrame):
	print("Computing Hydrogen Bond at Index: " + str(frameIndex))
	tic = time.clock()
	hbonds = md.baker_hubbard(traj, freq = 0.55, exclude_water = True, periodic = False)
	for hbond in hbonds:
		if (hbond[0] in boundIndices) and (hbond[2] in boundIndices):	
			results[tuple(hbond)] += 1
			resultsFrame[tuple(hbond)].append(frameIndex)
	toc = time.clock()
	print("Finished computing for ", str(frameIndex), " in time: ", str(toc - tic))


def calcHydrogenBondFramePairs(traj):
	print("Start calcHydrogenBondFramePairs")
	print("Number of frames in this trajectory is ", str(traj.n_frames))
	coord = traj.xyz
	allAtoms = traj.topology.select('all')
	results = defaultdict(int)
	resultsFrame = defaultdict(list)
	for i in range(traj.n_frames):
		boundIndices = []
		atoms_in_bound(coord[i], boundIndices, allAtoms)
		boundIndicesSet = frozenset(boundIndices)
		findHbonds(traj[i], boundIndices, results, i, resultsFrame)
	frameDict = {}
	for result in results:
		for tp in resultsFrame[result]:
			if(tp not in frameDict.keys()):
				frameDict[tp] = [atom_label(traj, list(result))]
			else:
				frameDict[tp].append(atom_label(traj, list(result)))
	hbFramePairs = []
	for time, framePair in frameDict.iteritems():
		print("Number of h-bond pairs at time: ", time, " = ", len(framePair))
		hbFramePairs.append(framePair)
	return hbFramePairs


