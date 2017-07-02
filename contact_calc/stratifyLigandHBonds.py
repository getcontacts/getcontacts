# Author: Anthony Kai Kwang Ma
# Email: anthonyma27@gmail.com, akma327@stanford.edu
# Date: 07/02/17
# stratifyLigandHBonds.py

### After computing hydrogen bonds between ligand and binding pocket residues for a fragment of 
### simulation, these functions allow user to stratify output into hbond subtypes 
### (ligand-backbone, ligand-sidechain, ligand-water, ligand water-mediated, ligand extended-water mediated)


from output_utils import *

def calcLigandBackBoneHBondResults(f, post_process_file, ligand, solventId = "HOH"):
	"""
		Compute interaction between ligand to backbone of residue
	"""
	fp = readFileOpener(post_process_file)
	currFrame, nFrames, computingTime, frameHBondList = -1, -1, -1, []
	backbone_atoms = ['N', 'O']
	for line in fp:
		isHeader, nFrames, currFrame, computingTime, frameHBondList = writeToHeader(f, line, nFrames, currFrame, computingTime, frameHBondList, '-hlb')
		if(isHeader == False):
			resAtoms = line.split(" -- ")
			resAtom1, resAtom2 = resAtoms[0].strip(), resAtoms[1].strip()
			if("LIG" in resAtom1 and ("LIG" not in resAtom2 and solventId not in resAtom2)): 
				atom2 = resAtom2.split("-")[1].strip()
				if(atom2 in backbone_atoms):
					frameHBondList.append(resAtom1 + " -- " + resAtom2)
			if(("LIG" not in resAtom1 and solventId not in resAtom1) and "LIG" in resAtom2):
				atom1 = resAtom1.split("-")[1].strip()
				if(atom1 in backbone_atoms):
					frameHBondList.append(resAtom1 + " -- " + resAtom2)
	writeToHBondFile(f, frameHBondList)
	f.write("nFrames:" + str(nFrames) + "\n")
	f.write("Computing Time:" + str(computingTime) + "\n")
	frameHBondList = []
			

def calcLigandSideChainHBondResults(f, post_process_file, ligand, solventId="HOH"):
	"""
		Compute interactions between ligand to sidechain of residue
	"""
	fp = readFileOpener(post_process_file)
	currFrame, nFrames, computingTime, frameHBondList = -1, -1, -1, []
	backbone_atoms = ['N', 'O']
	for line in fp:
		isHeader, nFrames, currFrame, computingTime, frameHBondList = writeToHeader(f, line, nFrames, currFrame, computingTime, frameHBondList, '-hls')
		if(isHeader == False):
			resAtoms = line.split(" -- ")
			resAtom1, resAtom2 = resAtoms[0].strip(), resAtoms[1].strip()
			if("LIG" in resAtom1 and ("LIG" not in resAtom2 and solventId not in resAtom2)): 
				atom2 = resAtom2.split("-")[1].strip()
				if(atom2 not in backbone_atoms):
					frameHBondList.append(resAtom1 + " -- " + resAtom2)
			if(("LIG" not in resAtom1 and solventId not in resAtom1) and "LIG" in resAtom2):
				atom1 = resAtom1.split("-")[1].strip()
				if(atom1 not in backbone_atoms):
					frameHBondList.append(resAtom1 + " -- " + resAtom2)
	writeToHBondFile(f, frameHBondList)
	f.write("nFrames:" + str(nFrames) + "\n")
	f.write("Computing Time:" + str(computingTime) + "\n")
	frameHBondList = []



def calcLigandToWaterBondResults(f, post_process_file, ligand, solventId="HOH"):
	"""
		Compute ligand to water interactions
	"""
	fp = readFileOpener(post_process_file)
	currFrame, nFrames, computingTime, frameHBondList = -1, -1, -1, []
	backbone_atoms = ['N', 'O']
	for line in fp:
		isHeader, nFrames, currFrame, computingTime, frameHBondList = writeToHeader(f, line, nFrames, currFrame, computingTime, frameHBondList, '-lw')
		if(isHeader == False):
			resAtoms = line.split(" -- ")
			resAtom1, resAtom2 = resAtoms[0].strip(), resAtoms[1].strip()
			if((solventId in resAtom1 and "LIG" in resAtom2) or ("LIG" in resAtom1 and solventId in resAtom2)):
				frameHBondList.append(resAtom1 + " -- " + resAtom2)
	writeToHBondFile(f, frameHBondList)
	f.write("nFrames:" + str(nFrames) + "\n")
	f.write("Computing Time:" + str(computingTime) + "\n")
	frameHBondList = []		


def calcLigandWaterBondResults(f, post_process_file, ligand, solventId="HOH"):
	"""
		Compute water-mediated interactions between ligand and residues in binding pocket
	"""
	fp = readFileOpener(post_process_file)
	currFrame, nFrames, computingTime, frameLWBondDict = -1, -1, -1, {}
	backbone_atoms = ['N', 'O']
	for line in fp:
		isHeader, nFrames, currFrame, computingTime, frameLWBondDict = writeToHeader(f, line, nFrames, currFrame, computingTime, frameLWBondDict, '-lwb')
		if(isHeader == False):
			resAtoms = line.split(" -- ")
			atomLabel1, atomLabel2, atomIndex1, atomIndex2 = resAtoms[0].strip(), resAtoms[1].strip(), resAtoms[2].strip(), resAtoms[3].strip()
			resAtom1 = (atomLabel1, atomIndex1)
			resAtom2 = (atomLabel2, atomIndex2)
			if(solventId in resAtom1[0] and solventId in resAtom2[0]): #skip if both atoms are water 
				continue 
			if(not(solventId not in resAtom1[0] and solventId not in resAtom2[0])): #at least one of atoms is water
				if(solventId in resAtom1[0]):
					solvent_atom, resid_ligand_atom = resAtom1, resAtom2
				else:
					solvent_atom, resid_ligand_atom = resAtom2, resAtom1
				if(solvent_atom not in frameLWBondDict):
					frameLWBondDict[solvent_atom] = [resid_ligand_atom]
				else:
					frameLWBondDict[solvent_atom].append(resid_ligand_atom)
	writeToLWBondFile(f, frameLWBondDict)
	f.write("nFrames:" + str(nFrames) + "\n")
	f.write("Computing Time:" + str(computingTime) + "\n")
	frameLWBondDict = {}


def calcLigandWaterBond2Results(f, post_process_file, ligand, solventId="HOH"):
	"""
		This method computes second degree water bridges. This class of interaction
		is defined between ligand and aa, which have unique hydrogen bond interactions 
		with hoh1 and hoh2 respectively. If hoh1-hoh2 also have hydrogen bond interaction 
		then ligand and aa form an extended second degree water bridge. 
	"""
	fp = readFileOpener(post_process_file)
	currFrame, nFrames, computingTime = -1, -1, -1
	frameLWBond2Dicts = ({},{}) 
	backbone_atoms = ['N', 'O']
	line_count = 0
	for line in fp: 
		isHeader, nFrames, currFrame, computingTime, frameLWBond2Dicts = writeToHeader(f, line, nFrames, currFrame, computingTime, frameLWBond2Dicts, '-lwb2')
		if(isHeader == False):
			resAtoms = line.split(" -- ")
			atomLabel1, atomLabel2, atomIndex1, atomIndex2 = resAtoms[0].strip(), resAtoms[1].strip(), resAtoms[2].strip(), resAtoms[3].strip()
			resAtom1 = (atomLabel1, atomIndex1)
			resAtom2 = (atomLabel2, atomIndex2)
			if(solventId in resAtom1[0] and solventId in resAtom2[0]):
				if(resAtom1 not in frameLWBond2Dicts[0]):
					frameLWBond2Dicts[0][resAtom1] = [resAtom2]
				else:
					frameLWBond2Dicts[0][resAtom1].append(resAtom2)
			elif((solventId in resAtom1[0] and solventId not in resAtom2[0]) or (solventId not in resAtom1[0] and solventId in resAtom2[0])):
				if(solventId in resAtom1[0]): solvent_atom, resid_atom = resAtom1, resAtom2
				else: solvent_atom, resid_atom = resAtom2, resAtom1 
				if(solvent_atom not in frameLWBond2Dicts[1]):
					frameLWBond2Dicts[1][solvent_atom] = [resid_atom]
				else:
					frameLWBond2Dicts[1][solvent_atom].append(resid_atom)
		line_count += 1
	writeToExtendedLWBondFile(f, frameLWBond2Dicts)
	f.write("nFrames:" + str(nFrames) + "\n")
	f.write("Computing Time:" + str(computingTime) + "\n")
	frameLWBond2Dicts = ({},{})

