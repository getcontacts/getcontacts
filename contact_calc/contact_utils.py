# Author: Anthony Kai Kwang Ma
# Email: anthonyma27@gmail.com, akma327@stanford.edu
# Date: 07/02/17
# contact_utils.py

### Molecular Dynamics Trajectory Simulation - Noncovalent Interaction Utilities ###


from __future__ import print_function, division
from itertools import product
import copy
import math
import numpy as np
from mdtraj.utils import ensure_type
from mdtraj.geometry import compute_distances, compute_angles
from mdtraj.geometry import _geometry


def dictToList(timeDict):
	"""
		Input: Dictionary mapping key to value 
		Output: List of sorted values 
	"""
	framePairs = []
	for time in timeDict.iterkeys():
		frameInteractions = timeDict[time]
		framePairs.append(frameInteractions)
	return framePairs

#take set of two points and return vector pointing from a to b
def pointsToVec(a, b):
	"""
		Input: point a and point b 
		Output: Vector ab
	"""
	returnVec = []
	for index, val in enumerate(a):
		returnVec.append(b[index] - a[index])
	return np.array(returnVec)

#distance between two points
def distBetweenTwoPoints(x, y):
	"""
		Return euclidean distance between two points 
	"""
	dist = 0
	for index in range(len(x)):
		dist += (x[index] - y[index])**2
	return math.sqrt(dist)

def retrieveAtomProperties(atom):
	"""
		Retrieve residue ID and the indicator
	"""
	props = str(atom).split("-")
	residID = props[0]
	atomIndicator = props[1]
	return(residID, atomIndicator)

def atomInfo(atom):
	"""
		Return atom name and residue 
	"""
	return (atom.name, atom.residue.name)


#provided trajectory frame and pair of atoms, returns distance
def pairAtomDistance(trajectoryFrame, atom1, atom2):
	"""
		Input: Trajectory frame and pair of atoms
		Output: Distance between the pair of atom in that particular frame
	"""
	atom1_loc = trajectoryFrame.xyz[0, atom1.index, :]
	atom2_loc = trajectoryFrame.xyz[0, atom2.index, :]
	pairVec = tuple(atom2_loc - atom1_loc)
	pairDist = length(pairVec)
	return pairDist

def appendAminoAcidToList(cand_list, atom, filterArray):
	"""
		Determine whether atom should be appended to candidate atom list. 
		(ie if filtering for OE then accepts OE1 and OE2 for distinct salt bridges)
	"""
	residID, atomIndicator = retrieveAtomProperties(atom[0])
	for elem in filterArray:
		if (elem in atomIndicator):
			cand_list.append(atom)
			break

def appendAminoAcidToDict(cand_Dict, atom, filterArray):
	"""
		Determine whether atom should be appended to candidate dictionary
	"""
	residID, atomIndicator = retrieveAtomProperties(atom[0])
	if(atomIndicator in filterArray):
		if residID not in cand_Dict.keys():
			cand_Dict[residID] = set()
			cand_Dict[residID].add(atom)
		else:
			cand_Dict[residID].add(atom)


def cross(a,b):
	"""
		Computes cross product of two vectors
	"""
	c = np.array([a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]) #used to not be np.array
	return c

def dotproduct(v1, v2):
	"""
		Computes dot product of two vectors
	"""
	return sum((a*b) for a, b in zip(v1, v2))

def length(v):
	"""
		Compute length of vector
	"""
	return math.sqrt(dotproduct(v,v))

def angleBtwnVec(v1,v2):
	"""
		Calculates angle between two vectors
	"""
	radBtwnVec = math.acos(dotproduct(v1,v2)/(length(v1)* length(v2)))
	return math.degrees(radBtwnVec)

def getTripletCoord(arom_resid_triplet, traj, time):
	"""
		Determine three equally spaced out points on an aromatic ring
	"""
	coord1 = tuple(traj.xyz[time, arom_resid_triplet[0][0].index, :])
	coord2 = tuple(traj.xyz[time, arom_resid_triplet[1][0].index, :])
	coord3 = tuple(traj.xyz[time, arom_resid_triplet[2][0].index, :])
	p1 = np.array([coord1[0], coord1[1], coord1[2]])
	p2 = np.array([coord2[0], coord2[1], coord2[2]])
	p3 = np.array([coord3[0], coord3[1], coord3[2]])
	return p1, p2, p3 

def calcCentroid(arom_resid_triplet, traj, time):
	"""
		Calculates centroid between three points
	"""
	p1, p2, p3 = getTripletCoord(arom_resid_triplet, traj, time)
	x = float(p1[0] + p2[0] + p3[0])/3
	y = float(p1[1] + p2[1] + p3[1])/3
	z = float(p1[2] + p2[2] + p3[2])/3
	centroid = np.array([x,y,z])
	return centroid 

def calcNormVec(arom_resid_triplet, traj, time):
	"""
		Calculates the normal vector for a plane
	"""
	p1, p2, p3 = getTripletCoord(arom_resid_triplet, traj, time)
	v1 = p3-p1
	v2 = p2-p1
	cp = np.cross(v1, v2)
	a, b, c = cp
	d = np.dot(cp, p3)
	planeCoord = np.array([a, b, c, d])
	return cp, planeCoord

def projectPointOntoPlane(planeCoord, point):
	"""
		Performs a projection of point onto 3D plane
	"""
	a, b, c, d = planeCoord[0], planeCoord[1], planeCoord[2], planeCoord[3]
	s, u, v, = point[0], point[1], point[2]
	t = (d - a*s - b*u - c*v) / (a*a + b*b + c*c)
	x = s + a*t
	y = u + b*t
	z = v + c*t
	return np.array([x,y,z])

def projectedCenterDist(planeCoord, center1, center2):
	"""
		Projects aromatic center2 onto plane of aromatic1 and determines
		the distance between the projected point and aromatic center1
	"""
	projectedCenter = projectPointOntoPlane(planeCoord, center2)
	projCenterDist = distBetweenTwoPoints(center1, projectedCenter)
	return projCenterDist

def removeDuplicate(cand_dict):
	"""
		Removes duplicate keys in dictionary
	"""
	for key in cand_dict.keys():
		atomStrSet = set()
		atomSet = set()
		val = cand_dict[key]
		for atom in val:
			if(str(atom) not in atomStrSet):
				atomStrSet.add(str(atom))
				atomSet.add(atom)
		cand_dict[key] = list(atomSet)

def calcPsiAngle(arom_center1, arom_center2, normVec1):
	"""
		Input: Aromatic center1 and center2 
		Output: Angle between normal_vector1 and vector(center2 - center1)
	"""
	centerToCenterVec = pointsToVec(arom_center2, arom_center1)
	psiAngle = angleBtwnVec(normVec1, centerToCenterVec)
	psiAngle = min(math.fabs(psiAngle - 0), math.fabs(psiAngle - 180))
	return psiAngle
