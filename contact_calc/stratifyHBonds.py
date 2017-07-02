# Author: Anthony Kai Kwang Ma
# Email: anthonyma27@gmail.com, akma327@stanford.edu
# Date: 07/02/17
# stratifyHBonds.py

from output_utils import *

### After computing hydrogen bonds for a fragment of simulation, these functions allow user to 
### Stratify output into hbond subtypes (backbone-backbone, backbone-sidechain, sidechain-sidechain
### residue-water, water-mediated, extended-water mediated)


def calcBackBoneHBondResults(f, post_process_file, solventId = 'HOH'):
	"""
		Compute interactions between two residues that are part of backbone
	"""
	fp = readFileOpener(post_process_file)
	currFrame, nFrames, computingTime, frameHBondList = -1, -1, -1, []
	backbone_atoms = ['N', 'O']
	for line in fp:
		isHeader, nFrames, currFrame, computingTime, frameHBondList = writeToHeader(f, line, nFrames, currFrame, computingTime, frameHBondList, '-hbbb')
		if(isHeader == False):
			resAtoms = line.split(" -- ")
			resAtom1, resAtom2 = resAtoms[0].strip(), resAtoms[1].strip()
			if(solventId not in resAtom1 and solventId not in resAtom2):
				atom1, atom2 = resAtom1.split("-")[1].strip(), resAtom2.split("-")[1].strip()
				if(atom1 in backbone_atoms and atom2 in backbone_atoms):
					frameHBondList.append(resAtom1 + " -- " + resAtom2)
	writeToHBondFile(f, frameHBondList)
	f.write("nFrames:" + str(nFrames) + "\n")
	f.write("Computing Time:" + str(computingTime) + "\n")
	frameHBondList = []


def calcSideChainBackBoneHBondResults(f, post_process_file, solventId='HOH'):
	"""
		Processes full water mediated hydrogen file and sees if both atoms are part of backbone
	"""
	fp = readFileOpener(post_process_file)
	currFrame, nFrames, computingTime, frameHBondList = -1, -1, -1, []
	backbone_atoms = ['N', 'O']
	for line in fp:
		isHeader, nFrames, currFrame, computingTime, frameHBondList = writeToHeader(f, line, nFrames, currFrame, computingTime, frameHBondList, '-hbsb')
		if(isHeader == False):
			resAtoms = line.split(" -- ")
			resAtom1, resAtom2 = resAtoms[0].strip(), resAtoms[1].strip()
			if(solventId not in resAtom1 and solventId not in resAtom2):
				atom1, atom2 = resAtom1.split("-")[1].strip(), resAtom2.split("-")[1].strip()
				if((atom1 in backbone_atoms and atom2 not in backbone_atoms) or (atom1 not in backbone_atoms and atom2 in backbone_atoms)):
					frameHBondList.append(resAtom1 + " -- " + resAtom2)
	writeToHBondFile(f, frameHBondList)
	f.write("nFrames:" + str(nFrames) + "\n")
	f.write("Computing Time:" + str(computingTime) + "\n")
	frameHBondList = []


def calcSideChainHBondResults(f, post_process_file, solventId='HOH'):
	"""
		Processes full water mediated hydrogen file and sees if exactly one 
		atom part of side chain and other atom part of back bone 
	"""
	fp = readFileOpener(post_process_file)
	currFrame, nFrames, computingTime, frameHBondList = -1, -1, -1, []
	backbone_atoms = ['N', 'O']
	for line in fp:
		isHeader, nFrames, currFrame, computingTime, frameHBondList = writeToHeader(f, line, nFrames, currFrame, computingTime, frameHBondList, '-hbss')
		if(isHeader == False):
			resAtoms = line.split(" -- ")
			resAtom1, resAtom2 = resAtoms[0].strip(), resAtoms[1].strip()
			if(solventId not in resAtom1 and solventId not in resAtom2):
				atom1, atom2 = resAtom1.split("-")[1].strip(), resAtom2.split("-")[1].strip()
				if(atom1 not in backbone_atoms and atom2 not in backbone_atoms):
					frameHBondList.append(resAtom1 + " -- " + resAtom2)
	writeToHBondFile(f, frameHBondList)
	f.write("nFrames:" + str(nFrames) + "\n")
	f.write("Computing Time:" + str(computingTime) + "\n")
	frameHBondList = []		


def calcResidueWaterBondResults(f, post_process_file, solventId='HOH'):
	"""
		Processes full water mediated hydrogen file and sees if exactly one 
		of the atoms is water and the other is part of residue 
	"""
	fp = readFileOpener(post_process_file)
	currFrame, nFrames, computingTime, frameHBondList = -1, -1, -1, []
	backbone_atoms = ['N', 'O']
	for line in fp:
		isHeader, nFrames, currFrame, computingTime, frameHBondList = writeToHeader(f, line, nFrames, currFrame, computingTime, frameHBondList, '-rw')
		if(isHeader == False):
			resAtoms = line.split(" -- ")
			resAtom1, resAtom2 = resAtoms[0].strip(), resAtoms[1].strip()
			if((solventId in resAtom1 and solventId not in resAtom2) or (solventId not in resAtom1 and solventId in resAtom2)):
				frameHBondList.append(resAtom1 + " -- " + resAtom2)
	writeToHBondFile(f, frameHBondList)
	f.write("nFrames:" + str(nFrames) + "\n")
	f.write("Computing Time:" + str(computingTime) + "\n")
	frameHBondList = []		


def calcWaterBondResults(f, post_process_file, solventId='HOH'):
	"""
		This method computes first degree water bridges which means there 
		is at least one water that has hydrogen bond interaction to both of 
		the amino acids simultaneously in a single frame. 
	"""
	fp = readFileOpener(post_process_file)
	currFrame, nFrames, computingTime, frameWBondDict = -1, -1, -1, {}
	backbone_atoms = ['N', 'O']
	for line in fp:
		isHeader, nFrames, currFrame, computingTime, frameWBondDict = writeToHeader(f, line, nFrames, currFrame, computingTime, frameWBondDict, '-wb')
		if(isHeader == False):
			resAtoms = line.split(" -- ")
			atomLabel1, atomLabel2, atomIndex1, atomIndex2 = resAtoms[0].strip(), resAtoms[1].strip(), resAtoms[2].strip(), resAtoms[3].strip()
			resAtom1 = (atomLabel1, atomIndex1)
			resAtom2 = (atomLabel2, atomIndex2)
			if(solventId in resAtom1[0] and solventId in resAtom2[0]): #skip if both atoms are water 
				continue 
			if(not(solventId not in resAtom1[0] and solventId not in resAtom2[0])): #at least one of atoms is water
				if(solventId in resAtom1[0]):
					solvent_atom, resid_atom = resAtom1, resAtom2
				else:
					solvent_atom, resid_atom = resAtom2, resAtom1
				if(solvent_atom not in frameWBondDict):
					frameWBondDict[solvent_atom] = [resid_atom]
				else:
					frameWBondDict[solvent_atom].append(resid_atom)
	writeToWBondFile(f, frameWBondDict)
	f.write("nFrames:" + str(nFrames) + "\n")
	f.write("Computing Time:" + str(computingTime) + "\n")
	frameWBondDict = {}


def calcWaterBond2Results(f, post_process_file, solventId='HOH'):
	"""
		This method computes second degree water bridges. This class of interaction
		is defined between amino acids aa1 and aa2, which have unique hydrogen bond
		interactions with hoh1 and hoh2 respectively. If hoh1-hoh2 also have hydrogen 
		bond interaction then aa1 and aa2 form an extended second degree water bridge. 
	"""
	fp = readFileOpener(post_process_file)
	currFrame, nFrames, computingTime = -1, -1, -1
	frameWBond2Dicts = ({},{}) 
	backbone_atoms = ['N', 'O']
	line_count = 0
	for line in fp: 
		if(line_count %100 == 0): print("Line count: " + str(line_count))
		isHeader, nFrames, currFrame, computingTime, frameWBond2Dicts = writeToHeader(f, line, nFrames, currFrame, computingTime, frameWBond2Dicts, '-wb2')
		if(isHeader == False):
			resAtoms = line.split(" -- ")
			atomLabel1, atomLabel2, atomIndex1, atomIndex2 = resAtoms[0].strip(), resAtoms[1].strip(), resAtoms[2].strip(), resAtoms[3].strip()
			resAtom1 = (atomLabel1, atomIndex1)
			resAtom2 = (atomLabel2, atomIndex2)
			if(solventId in resAtom1[0] and solventId in resAtom2[0]):
				if(resAtom1 not in frameWBond2Dicts[0]):
					frameWBond2Dicts[0][resAtom1] = [resAtom2]
				else:
					frameWBond2Dicts[0][resAtom1].append(resAtom2)
			elif((solventId in resAtom1[0] and solventId not in resAtom2[0]) or (solventId not in resAtom1[0] and solventId in resAtom2[0])):
				if(solventId in resAtom1[0]): solvent_atom, resid_atom = resAtom1, resAtom2
				else: solvent_atom, resid_atom = resAtom2, resAtom1 
				if(solvent_atom not in frameWBond2Dicts[1]):
					frameWBond2Dicts[1][solvent_atom] = [resid_atom]
				else:
					frameWBond2Dicts[1][solvent_atom].append(resid_atom)
		line_count += 1
	writeToExtendedWBondFile(f, frameWBond2Dicts)
	f.write("nFrames:" + str(nFrames) + "\n")
	f.write("Computing Time:" + str(computingTime) + "\n")
	frameWBond2Dicts = ({},{})



