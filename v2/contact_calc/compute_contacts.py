##############################################################################
# MDContactNetworks: A Python Library for computing non-covalent contacts
#                    throughout Molecular Dynamics Trajectories. 
# Copyright 2016-2017 Stanford University and the Authors
#
# Authors: Anthony Kai Kwang Ma
# Email: anthony.ma@yale.edu, anthonyma27@gmail.com, akma327@stanford.edu
##############################################################################

##############################################################################
# Imports
##############################################################################

from vmd import *
import molecule
from multiprocessing import *
import MDAnalysis as mda
from contact_utils import *
from hbonds import *
from salt_bridges import *
from pi_cation import *
from vanderwaals import *

##############################################################################
# Global Variables
##############################################################################

TRAJ_FRAG_SIZE = 1000

##############################################################################
# Functions
##############################################################################


def compute_frame_contacts(traj_frag_molid, frame_idx, ITYPES, solvent_resn, chain_id, ligand, index_to_label):
	"""
	Computes each of the specified non-covalent interaction type for a single frame

	Parameters
	----------
	traj_frag_molid: int
		Identifier to simulation fragment in VMD
	frame_idx: int
		Frame number to query
	ITYPES: list
		Denotes the list of non-covalent interaction types to compute contacts for 
	solvent_resn: string, default = TIP3
		Denotes the resname of solvent in simulation
	chain_id: string, default = None
		Specify chain of protein to perform computation on 
	ligand: string, default = None
		Include ligand resname if computing contacts between ligand and binding pocket residues
	index_to_label: dict 
		Maps VMD atom index to label "chain:resname:resid:name:index"
		{11205: "A:ASP:114:CA:11205, ...}

	Returns
	-------
	frame_contacts: list of tuples, [(frame_index, atom1_label, atom2_label, itype), ...]

	"""

	### NOTE
	# Eventually we will pass in 
	# prep_computation = {"index_to_label": index_to_label, "salt_bridge": (anion_list, cation_list), etc}
	# And pull feed these pre computed information into the calculator. This avoids repeating the 
	# same calculation throughout simulation. 

	frame_contacts = []

	if("-sb" in ITYPES):
		anion_list, cation_list = prep_salt_bridge_computation(traj_frag_molid, frame_idx, chain_id)
		frame_contacts += compute_salt_bridges(traj_frag_molid, frame_idx, anion_list, cation_list)
	if("-pc" in ITYPES):
		frame_contacts += compute_pi_cation(traj_frag_molid, frame_idx, index_to_label, chain_id)
	# if("-ps" in ITYPES):
	# 	frame_contacts += compute_pi_stacking(traj_frag_molid, frame_idx, index_to_label, chain_id)
	# if("-ts" in ITYPES):
	# 	frame_contacts += compute_t_stacking(traj_frag_molid, frame_idx, index_to_label, chain_id)
	if("-vdw" in ITYPES):
		frame_contacts += compute_vanderwaals(traj_frag_molid, frame_idx, index_to_label, chain_id)
	if("-hb" in ITYPES):
		frame_contacts += compute_hydrogen_bonds(traj_frag_molid, frame_idx, index_to_label, solvent_resn, chain_id)
	if("-lhb" in ITYPES):
		frame_contacts += compute_hydrogen_bonds(traj_frag_molid, frame_idx, index_to_label, solvent_resn, chain_id, ligand)

	# print("Finished computing contacts")
	return frame_contacts

def compute_fragment_contacts(frag_idx, beg_frame, end_frame, TOP, TRAJ, ITYPES, stride, solvent_resn, chain_id, ligand, index_to_label):
	""" 
	Reads in a single trajectory fragment and calls compute_frame_contacts on each frame

	Parameters
	----------
	frag_idx: int
		Trajectory fragment index for worker to keep track of
	beg_frame: int
		Start frame of trajectory fragment
	end_frame: int
		End frame of trajectory fragment
	TOP: Topology
		In .pdb or .mae format
	TRAJ: Trajectory
		In .nc or .dcd format
	ITYPES: list
		Denotes the list of non-covalent interaction types to compute contacts for 
	stride: int, default = 1
		Frequency to skip frames in trajectory
	solvent_resn: string, default = TIP3
		Denotes the resname of solvent in simulation
	chain_id: string, default = None
		Specify chain of protein to perform computation on 
	ligand: string, default = None
		Include ligand resname if computing contacts between ligand and binding pocket residues
	index_to_label: dict 
		Maps VMD atom index to label "chain:resname:resid:name:index"
		{11205: "A:ASP:114:CA:11205, ...}

	Return
	------
	frag_idx: int 

	fragment_contacts: list of tuples, [(frame_index, atom1_label, atom2_label, itype), ...]

	"""

	traj_frag_molid = load_traj(TOP, TRAJ, beg_frame, end_frame, stride)
	fragment_contacts = []

	### Compute contacts for each frame
	num_frag_frames = molecule.numframes(traj_frag_molid)
	for frame_idx in range(1, num_frag_frames):
		# if(frame_idx > 1): break
		fragment_contacts += compute_frame_contacts(traj_frag_molid, frame_idx, ITYPES, solvent_resn, chain_id, ligand, index_to_label)

	### Delete trajectory fragment to clear memory
	molecule.delete(traj_frag_molid)
	return (frag_idx, fragment_contacts)

def compute_fragment_contacts_helper(args):
	return compute_fragment_contacts(*args)

def compute_dynamic_contacts(TOP, TRAJ, OUTPUT_DIR, ITYPES, cores, stride, solvent_resn, chain_id, ligand):
	""" 
	Computes non-covalent contacts across the entire trajectory and writes to output.

	Parameters
	----------
	TOP: Topology
		In .pdb or .mae format
	TRAJ: Trajectory
		In .nc or .dcd format
	OUTPUT_DIR: string
		Absolute path to output directory 
	ITYPES: list
		Denotes the list of non-covalent interaction types to compute contacts for 
	cores: int, default = 6
		Number of CPU cores to parallelize over
	stride: int, default = 1
		Frequency to skip frames in trajectory
	solvent_resn: string, default = TIP3
		Denotes the resname of solvent in simulation
	chain_id: string, default = None
		Specify chain of protein to perform computation on 
	ligand: string, default = None
		Include ligand resname if computing contacts between ligand and binding pocket residues

	"""

	index_to_label = gen_index_to_atom_label(TOP, TRAJ)

	sim_length = len(mda.Universe(TOP, TRAJ).trajectory) # Temporary command. Doesn't work for mae 
	input_args = []

	### Serial 
	# output = []
	# for frag_idx, beg_frame in enumerate(range(0, sim_length, TRAJ_FRAG_SIZE)):
	# 	if(frag_idx > 0): break
	# 	end_frame = beg_frame + TRAJ_FRAG_SIZE
	# 	output += compute_fragment_contacts(frag_idx, beg_frame, end_frame, TOP, TRAJ, ITYPES, stride, solvent_resn, chain_id, ligand, index_to_label)[1]

	# for o in output: print o
	# print(len(output))

	### Generate input arguments for each trajectory piece
	for frag_idx, beg_frame in enumerate(range(0, sim_length, TRAJ_FRAG_SIZE)):
		# if(frag_idx > 0): break
		end_frame = beg_frame + TRAJ_FRAG_SIZE
		input_args.append((frag_idx, beg_frame, end_frame, TOP, TRAJ, ITYPES, stride, solvent_resn, chain_id, ligand, index_to_label))

	pool = Pool(processes=cores)
	output = pool.map(compute_fragment_contacts_helper, input_args)
	pool.close()
	pool.join()


	### Sort output by trajectory fragments
	output = sorted(output, key = lambda x: (x[0]))
	
	### Assign absolute frame indices 
	# print("Assigning absolute frame indices")
	contact_types = set()
	full_output = []
	num_frames = 0
	for frag_idx, contacts in output:
		for c in contacts:
			c[0] = num_frames + c[0] - 1
			full_output.append(c)
			contact_types.add(c[-1]) ### Keep track of all itypes
		frag_len = len(set([c[0] for c in contacts]))
		num_frames += frag_len

	### Writing output to seperate files, one for each itype
	# print("Writing to output")
	if not os.path.exists(OUTPUT_DIR):
		os.makedirs(OUTPUT_DIR)

	print(contact_types)

	fd_map = {itype: open(OUTPUT_DIR + "/" + itype.strip("-") + ".txt", 'w') for itype in contact_types}
	for contact in full_output:
		itype_key = contact[-1]
		fd_map[itype_key].write("\t".join(map(str, contact)) + "\n")

