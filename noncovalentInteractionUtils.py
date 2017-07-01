# Molecular Dynamics Trajectory Simulation - Noncovalent Interaction Utilities
# 10/12/15 - Anthony Ma 

from __future__ import print_function, division
import math
from itertools import product
import numpy as np
from mdtraj.utils import ensure_type
from mdtraj.geometry import compute_distances, compute_angles
from mdtraj.geometry import _geometry
import copy

# Convert dictionary with key and value as list to a list of these values
def dictToList(timeDict):
	print("Converting dict to list")
	sbFramePairs = []
	for time in timeDict.iterkeys():
		frameInteractions = timeDict[time]
		sbFramePairs.append(frameInteractions)
	return sbFramePairs

#take set of two points and return vector pointing from a to b
def pointsToVec(a, b):
	returnVec = []
	for index, val in enumerate(a):
		returnVec.append(b[index] - a[index])
	return np.array(returnVec)

#distance between two points
def distBetweenTwoPoints(x, y):
	dist = 0
	for index in range(len(x)):
		dist += (x[index] - y[index])**2
	return math.sqrt(dist)

#retrieve residue ID and the indicator
def retrieveAtomProperties(atom):
	props = str(atom).split("-")
	residID = props[0]
	atomIndicator = props[1]
	return(residID, atomIndicator)

def atomInfo(atom):
	return (atom.name, atom.residue.name)


#provided trajectory frame and pair of atoms, returns distance
def pairAtomDistance(trajectoryFrame, atom1, atom2):
	atom1_loc = trajectoryFrame.xyz[0, atom1.index, :]
	atom2_loc = trajectoryFrame.xyz[0, atom2.index, :]
	pairVec = tuple(atom2_loc - atom1_loc)
	pairDist = length(pairVec)
	return pairDist

# Check whether atom is should be appended to candidate list
# If filter is for OE then accept OE1 and OE2 for distinct salt bridges
def appendAminoAcidToList(cand_list, atom, filterArray):
	residID, atomIndicator = retrieveAtomProperties(atom[0])
	for elem in filterArray:
		if (elem in atomIndicator):
			cand_list.append(atom)
			break

# Check whether atom should be appended to candidate dictionary
# Purpose is for storing multiple crucial atoms for given residue grouping
def appendAminoAcidToDict(cand_Dict, atom, filterArray):
	residID, atomIndicator = retrieveAtomProperties(atom[0])
	if(atomIndicator in filterArray):
		if residID not in cand_Dict.keys():
			cand_Dict[residID] = set()
			cand_Dict[residID].add(atom)
		else:
			cand_Dict[residID].add(atom)

# cross product of two vectors
def cross(a,b):
	c = np.array([a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]) #used to not be np.array
	return c

def dotproduct(v1, v2):
	return sum((a*b) for a, b in zip(v1, v2))

def length(v):
	return math.sqrt(dotproduct(v,v))

def angleBtwnVec(v1,v2):
	radBtwnVec = math.acos(dotproduct(v1,v2)/(length(v1)* length(v2)))
	return math.degrees(radBtwnVec)


#returns the coordinate of the triplet of atoms 
# def getTripletCoord(arom_resid_triplet, trajectoryFrame):
# 	coord1 = tuple(trajectoryFrame.xyz[0, arom_resid_triplet[0][0].index, :])
# 	coord2 = tuple(trajectoryFrame.xyz[0, arom_resid_triplet[1][0].index, :])
# 	coord3 = tuple(trajectoryFrame.xyz[0, arom_resid_triplet[2][0].index, :])
# 	p1 = np.array([coord1[0], coord1[1], coord1[2]])
# 	p2 = np.array([coord2[0], coord2[1], coord2[2]])
# 	p3 = np.array([coord3[0], coord3[1], coord3[2]])
# 	return p1, p2, p3 


# def calcNormVec(arom_resid_triplet, trajectoryFrame):
# 	p1, p2, p3 = getTripletCoord(arom_resid_triplet, trajectoryFrame) #np arrays
# 	v1 = p3-p1
# 	v2 = p2-p1
# 	cp = np.cross(v1, v2)
# 	a, b, c = cp
# 	d = np.dot(cp, p3)
# 	planeCoord = np.array([a, b, c, d]) #used to not be np.array
# 	return cp, planeCoord


#calculate centroid of aromatic ring 
# def calcCentroid(arom_resid_triplet, trajectoryFrame):
# 	p1, p2, p3 = getTripletCoord(arom_resid_triplet, trajectoryFrame)
# 	x = float(p1[0] + p2[0] + p3[0])/3
# 	y = float(p1[1] + p2[1] + p3[1])/3
# 	z = float(p1[2] + p2[2] + p3[2])/3
# 	centroid = np.array([x, y, z])
# 	return centroid

def getTripletCoord(arom_resid_triplet, traj, time):
	coord1 = tuple(traj.xyz[time, arom_resid_triplet[0][0].index, :])
	coord2 = tuple(traj.xyz[time, arom_resid_triplet[1][0].index, :])
	coord3 = tuple(traj.xyz[time, arom_resid_triplet[2][0].index, :])
	p1 = np.array([coord1[0], coord1[1], coord1[2]])
	p2 = np.array([coord2[0], coord2[1], coord2[2]])
	p3 = np.array([coord3[0], coord3[1], coord3[2]])
	return p1, p2, p3 

def calcCentroid(arom_resid_triplet, traj, time):
	p1, p2, p3 = getTripletCoord(arom_resid_triplet, traj, time)
	x = float(p1[0] + p2[0] + p3[0])/3
	y = float(p1[1] + p2[1] + p3[1])/3
	z = float(p1[2] + p2[2] + p3[2])/3
	centroid = np.array([x,y,z])
	return centroid 

def calcNormVec(arom_resid_triplet, traj, time):
	p1, p2, p3 = getTripletCoord(arom_resid_triplet, traj, time) #np arrays
	v1 = p3-p1
	v2 = p2-p1
	cp = np.cross(v1, v2)
	a, b, c = cp
	d = np.dot(cp, p3)
	planeCoord = np.array([a, b, c, d]) #used to not be np.array
	return cp, planeCoord



def projectPointOntoPlane(planeCoord, point):
	a, b, c, d = planeCoord[0], planeCoord[1], planeCoord[2], planeCoord[3]
	s, u, v, = point[0], point[1], point[2]
	t = (d - a*s - b*u - c*v) / (a*a + b*b + c*c)
	x = s + a*t
	y = u + b*t
	z = v + c*t
	return np.array([x,y,z]) #used to not be np.array


#projects the aromatic center 2 onto aromatic 1 plane and finds the
#distance between projection and aromatic center 1 
def projectedCenterDist(planeCoord, center1, center2):
	projectedCenter = projectPointOntoPlane(planeCoord, center2)
	projCenterDist = distBetweenTwoPoints(center1, projectedCenter)
	return projCenterDist


#remove duplicates from dictionary
def removeDuplicate(cand_dict):
	for key in cand_dict.keys():
		atomStrSet = set()
		atomSet = set()
		val = cand_dict[key]
		for atom in val:
			if(str(atom) not in atomStrSet):
				atomStrSet.add(str(atom))
				atomSet.add(atom)
		cand_dict[key] = list(atomSet)


#method takes in aromatic center 1 and 2 and will compute the angle
#between normvector1 and the vector(center2 - center1)
def calcPsiAngle(arom_center1, arom_center2, normVec1):
	centerToCenterVec = pointsToVec(arom_center2, arom_center1)
	psiAngle = angleBtwnVec(normVec1, centerToCenterVec)
	psiAngle = min(math.fabs(psiAngle - 0), math.fabs(psiAngle - 180))
	return psiAngle

# other methods for gamma

#method 5: treat one aromatic center as cation and use same method as in pi cation
#similar in results to method 4 but more restrictive. 
# psiAngle1 = calcPsiAngle(arom_center1, arom_center2, normVec1)
# psiAngle2 = calcPsiAngle(arom_center2, arom_center1, normVec2)
# psiAngle = min(psiAngle1, psiAngle2)
# print("psi angles: ", psiAngle1, psiAngle2)
# if(psiAngle < T_STACK_PSI_ANGLE_CUTOFF):
# 	pairInfo = [aromatic1, aromatic2, centroidDist, perpAngleDeviation, psiAngle]
# 	faceEdgePairs.append(pairInfo)

#method 4: project centers onto first plane and see if close to other
# projCenterDist1 = projectedCenterDist(planeCoord1, arom_center1, arom_center2)
# projCenterDist2 = projectedCenterDist(planeCoord2, arom_center2, arom_center1)
# projCenterDist = min(projCenterDist1, projCenterDist2)
# print("projected center dist: ", projCenterDist1, projCenterDist2)
# if(projCenterDist < PI_PI_CENTER_DIST_TOLERANCE):
# 	pairInfo = [aromatic1, aromatic2, centroidDist, alignAngle, projCenterDist]
# 	faceFacePairs.append(pairInfo)

# takes in the centroid and plane coordinates of one aromatic ring, and projects
# the points of the second ring onto this plane. Return average distance from
# the center
# def calcProjPointsDeviation(planeCoord, arom_center, arom_resid_triplet, trajectoryFrame):
# 	p1, p2, p3 = getTripletCoord(arom_resid_triplet, trajectoryFrame)
# 	totDist = 0
# 	for point in [p1, p2, p3]:
# 		totDist += distBetweenTwoPoints(point, arom_center)
# 	return float(totDist)/3

# Finds the angle between center1 to center2 vector and the center1 projection
# of center2 onto aromatic plane1 vector. This would give you the angle of elevation
# that planes are corresponding to each other. 90 degree elevation
# would mean one is above the other and 0 degree means side by side
# def calcGammaAngle(planeCoord, arom_center1, arom_center2):
# 	projectedPoint = projectPointOntoPlane(planeCoord, arom_center2)
# 	v1 = pointsToVec(arom_center1, projectedPoint)
# 	v2 = pointsToVec(arom_center1, arom_center2)
# 	return angleBtwnVec(v1, v2)

#method 3
# gammaAngle = calcGamma(normVec1, normVec2, arom_center1, arom_center2)
# if(gammaAngle <= CUTOFF_GAMMA):
# 	tempDict[(k1, k2)] = 1
# 	pairInfo = [aromatic1, aromatic2, centroidDist, alignAngle, gammaAngle]
# 	faceFacePairs.append(pa

