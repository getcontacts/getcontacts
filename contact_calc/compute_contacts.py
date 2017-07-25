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
import os
import molecule
import datetime
import glob
from multiprocessing import *
from contact_utils import *
from aromatics import *
from hbonds import *
from salt_bridges import *
from pi_cation import *
from vanderwaals import *

##############################################################################
# Global Variables
##############################################################################
TRAJ_FRAG_SIZE = 100

##############################################################################
# Functions
##############################################################################

### Note: Trying to write directly to output and doing postprocessing after
### May save a lot of memory instead of having massive arrays, but there is 
### also the IO time (which would happen anyways). Figure out best output
### format and most efficient way to write to disk. 

def compute_frame_contacts(traj_frag_molid, frag_idx, frame_idx, ITYPES, geom_criterion_values, solvent_resn, sele_id, ligand, index_to_label):
	"""
	Computes each of the specified non-covalent interaction type for a single frame

	Parameters
	----------
	traj_frag_molid: int
		Identifier to simulation fragment in VMD
	frag_idx: int
		Trajectory fragment index for worker to keep track of
	frame_idx: int
		Frame number to query
	ITYPES: list
		Denotes the list of non-covalent interaction types to compute contacts for 
	geom_criterion_values: dict
		Dictionary containing the cutoff values for all geometric criteria
	solvent_resn: string, default = TIP3
		Denotes the resname of solvent in simulation
	sele_id: string, default = None
		Compute contacts on subset of atom selection based on VMD query
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
	tic = datetime.datetime.now()

	### Extract geometric criterion 
	SALT_BRIDGE_CUTOFF_DISTANCE = geom_criterion_values['SALT_BRIDGE_CUTOFF_DISTANCE']
	PI_CATION_CUTOFF_DISTANCE = geom_criterion_values['PI_CATION_CUTOFF_DISTANCE']
	PI_CATION_CUTOFF_ANGLE = geom_criterion_values['PI_CATION_CUTOFF_ANGLE']
	PI_STACK_CUTOFF_DISTANCE = geom_criterion_values['PI_STACK_CUTOFF_DISTANCE']
	PI_STACK_CUTOFF_ANGLE = geom_criterion_values['PI_STACK_CUTOFF_ANGLE']
	PI_STACK_PSI_ANGLE = geom_criterion_values['PI_STACK_PSI_ANGLE']
	T_STACK_CUTOFF_DISTANCE = geom_criterion_values['T_STACK_CUTOFF_DISTANCE']
	T_STACK_CUTOFF_ANGLE = geom_criterion_values['T_STACK_CUTOFF_ANGLE']
	T_STACK_PSI_ANGLE = geom_criterion_values['T_STACK_PSI_ANGLE']
	HBOND_CUTOFF_DISTANCE = geom_criterion_values['HBOND_CUTOFF_DISTANCE']
	HBOND_CUTOFF_ANGLE = geom_criterion_values['HBOND_CUTOFF_ANGLE']
	VDW_EPSILON = geom_criterion_values['VDW_EPSILON']
	

	frame_contacts = []
	if("-sb" in ITYPES):
		frame_contacts += compute_salt_bridges(traj_frag_molid, frame_idx, sele_id, SALT_BRIDGE_CUTOFF_DISTANCE)
	if("-pc" in ITYPES):
		frame_contacts += compute_pi_cation(traj_frag_molid, frame_idx, index_to_label, sele_id, PI_CATION_CUTOFF_DISTANCE, PI_CATION_CUTOFF_ANGLE)
	if("-ps" in ITYPES):
		frame_contacts += compute_pi_stacking(traj_frag_molid, frame_idx, index_to_label, sele_id, PI_STACK_CUTOFF_DISTANCE, PI_STACK_CUTOFF_ANGLE, PI_STACK_PSI_ANGLE)
	if("-ts" in ITYPES):
		frame_contacts += compute_t_stacking(traj_frag_molid, frame_idx, index_to_label, sele_id, T_STACK_CUTOFF_DISTANCE, T_STACK_CUTOFF_ANGLE, T_STACK_PSI_ANGLE)
	if("-vdw" in ITYPES):
		frame_contacts += compute_vanderwaals(traj_frag_molid, frame_idx, index_to_label, sele_id, VDW_EPSILON)
	if("-hb" in ITYPES):
		frame_contacts += compute_hydrogen_bonds(traj_frag_molid, frame_idx, index_to_label, solvent_resn, sele_id, None, HBOND_CUTOFF_DISTANCE, HBOND_CUTOFF_ANGLE)
	if("-lhb" in ITYPES):
		frame_contacts += compute_hydrogen_bonds(traj_frag_molid, frame_idx, index_to_label, solvent_resn, sele_id, ligand, HBOND_CUTOFF_DISTANCE, HBOND_CUTOFF_ANGLE)

	toc = datetime.datetime.now()
	print("Finished computing contacts for frag %s frame %s in %s s" % (frag_idx, frame_idx, (toc-tic).total_seconds()))
	return frame_contacts

def compute_fragment_contacts(frag_idx, beg_frame, end_frame, TOP, TRAJ, OUTPUT_DIR, contact_types, ITYPES, geom_criterion_values, stride, solvent_resn, sele_id, ligand, index_to_label):
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
	geom_criterion_values: dict
		Dictionary containing the cutoff values for all geometric criteria
	stride: int, default = 1
		Frequency to skip frames in trajectory
	solvent_resn: string, default = TIP3
		Denotes the resname of solvent in simulation
	sele_id: string, default = None
		Compute contacts on subset of atom selection based on VMD query
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
		fragment_contacts += compute_frame_contacts(traj_frag_molid, frag_idx, frame_idx, ITYPES, geom_criterion_values, solvent_resn, sele_id, ligand, index_to_label)

	### Delete trajectory fragment to clear memory
	molecule.delete(traj_frag_molid)

	### Write directly out to temporary output 
	print("Writing output to seperate files, one for each itype ...")
	
	fd_map = {itype: open(OUTPUT_DIR + "/" + itype.strip("-") + "_frag_" + str(frag_idx) + ".txt", 'w') for itype in contact_types}
	for contact in fragment_contacts:
		itype_key = contact[-1]
		output_string = str(frag_idx) + "\t" + "\t".join(map(str, contact)) + "\n"
		fd_map[itype_key].write(output_string)

	for itype in fd_map:
		fd_map[itype].close()
	return (frag_idx, num_frag_frames - 1)

def compute_fragment_contacts_helper(args):
	return compute_fragment_contacts(*args)

def stitch_fragment_contacts(itype, OUTPUT_DIR, frag_contact_files, frag_idx_to_length):
	"""
	Stitch together multiple fragments of non-covalent contacts into
	single file and delete individual fragment files. 

	Parameters
	----------
	itype: Type of non-covalent contact
	OUTPUT_DIR: string
		Absolute path to output directory 
	frag_contact_files: list of strings
		List of paths to fragment contact files
	frag_idx_to_length: dict from int to int 
		Map the fragment index to length of fragment 
	"""
	print("Stitching %s ..." % (itype))
	fo = open(OUTPUT_DIR + "/" + itype + ".txt", 'w')

	num_frames = 0
	for frag_contact_file in frag_contact_files:
		frag_idx = int(frag_contact_file.split("/")[-1].strip(".txt").split("_")[2])
		ffrag = open(frag_contact_file, 'r')		
		for line in ffrag:
			linfo = line.split("\t")
			frame_idx = int(linfo[1])
			new_frame_idx = num_frames + frame_idx - 1
			output_string = str(new_frame_idx) + "\t" + "\t".join(linfo[2:])
			fo.write(output_string)
		num_frames += frag_idx_to_length[frag_idx]
		ffrag.close()
		os.remove(frag_contact_file)

	fo.close()


def compute_dynamic_contacts(TOP, TRAJ, OUTPUT_DIR, ITYPES, geom_criterion_values, cores, stride, solvent_resn, sele_id, ligand):
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
	geom_criterion_values: dict
		Dictionary containing the cutoff values for all geometric criteria
	cores: int, default = 6
		Number of CPU cores to parallelize over
	stride: int, default = 1
		Frequency to skip frames in trajectory
	solvent_resn: string, default = TIP3
		Denotes the resname of solvent in simulation
	sele_id: string, default = None
		Compute contacts on subset of atom selection based on VMD query
	chain_id: string, default = None
		Specify chain of protein to perform computation on 
	ligand: string, default = None
		Include ligand resname if computing contacts between ligand and binding pocket residues

	"""
	### Set up file descriptors for writing output
	if not os.path.exists(OUTPUT_DIR):
		os.makedirs(OUTPUT_DIR)

	contact_types = []
	for itype in ITYPES:
		if(itype == "-hb"):
			contact_types += ["hbbb", "hbsb", "hbss", "wb", "wb2"]
		elif(itype == "-lhb"):
			contact_types += ["hls", "hlb", "lwb", "lwb2"]
		else:
			contact_types += [itype.strip("-")]


	index_to_label = gen_index_to_atom_label(TOP, TRAJ)
	sim_length = estimate_simulation_length(TOP, TRAJ)
	input_args = []

	### Serial 
	# output = []
	# for frag_idx, beg_frame in enumerate(range(0, sim_length, TRAJ_FRAG_SIZE)):
	# 	if(frag_idx > 0): break
	# 	end_frame = beg_frame + TRAJ_FRAG_SIZE
	# 	frag_idx, frag_length = compute_fragment_contacts(frag_idx, beg_frame, end_frame, TOP, TRAJ, OUTPUT_DIR, contact_types, ITYPES, geom_criterion_values, stride, solvent_resn, sele_id, ligand, index_to_label)
	# 	output.append((frag_idx, frag_length))

	### Generate input arguments for each trajectory piece
	print("MDContactNetworks processing TRAJ: %s with %s total frames with stride %s ..." % (TRAJ, str(sim_length), str(stride)))
	for frag_idx, beg_frame in enumerate(range(0, sim_length, TRAJ_FRAG_SIZE)):
		# if(frag_idx > 0): break
		end_frame = beg_frame + TRAJ_FRAG_SIZE
		print("Processing fragment %s beg_frame %s end_frame %s" % (frag_idx, beg_frame, end_frame))
		input_args.append((frag_idx, beg_frame, end_frame, TOP, TRAJ, OUTPUT_DIR, contact_types, ITYPES, geom_criterion_values, stride, solvent_resn, sele_id, ligand, index_to_label))

	pool = Pool(processes=cores)
	output = pool.map(compute_fragment_contacts_helper, input_args)
	pool.close()
	pool.join()

	### Sort output by trajectory fragments
	print("Map fragments to length")
	frag_idx_to_length = {}
	output = sorted(output, key = lambda x: (x[0]))
	for frag_idx, frag_length in output:
		frag_idx_to_length[frag_idx] = frag_length
		print(frag_idx, frag_length)

	### Combine fragments to single large stitched files
	for itype in contact_types:
		frag_contact_files = glob.glob(OUTPUT_DIR + "/" + itype + "_frag*")
		frag_contact_files.sort(key=natural_keys)
		stitch_fragment_contacts(itype, OUTPUT_DIR, frag_contact_files, frag_idx_to_length)


	

