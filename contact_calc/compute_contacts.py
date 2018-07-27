############################################################################
# Copyright 2018 Anthony Ma & Stanford University                          #
#                                                                          #
# Licensed under the Apache License, Version 2.0 (the "License");          #
# you may not use this file except in compliance with the License.         #
# You may obtain a copy of the License at                                  #
#                                                                          #
#     http://www.apache.org/licenses/LICENSE-2.0                           #
#                                                                          #
# Unless required by applicable law or agreed to in writing, software      #
# distributed under the License is distributed on an "AS IS" BASIS,        #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. #
# See the License for the specific language governing permissions and      #
# limitations under the License.                                           #
############################################################################

##############################################################################
# Imports
##############################################################################

import datetime
from multiprocessing import Process, Queue
from vmd import *  # Loads the static `molecule` object

from .contact_utils import *
from .aromatics import *
from .hbonds import *
from .salt_bridges import *
from .pi_cation import *
from .vanderwaals import *
from .hydrophobics import compute_hydrophobics

##############################################################################
# Global Variables
##############################################################################
TRAJ_FRAG_SIZE = 100

##############################################################################
# Functions
##############################################################################


def compute_frame_contacts(traj_frag_molid, frame_idx, ITYPES, geom_criterion_values, solvent_resn, sele1, sele2,
                           sele1_atoms, sele2_atoms, index_to_atom):
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
    sele1: string, default = None
        Compute contacts on subset of atom selection based on VMD query
    sele2: string, default = None
        If second VMD query is specified, then compute contacts between atom selection 1 and 2 
    sele1_atoms: list 
        List of atom label indices for all atoms in selection 1
    sele2_atoms: list 
        List of atom label indices for all atoms in selection 2
    chain_id: string, default = None
        Specify chain of protein to perform computation on 
    index_to_atom: dict
        Maps VMD atom index to Atom

    Returns
    -------
    frame_contacts: list of tuples, [(frame_index, atom1_label, atom2_label, itype), ...]

    """

    # Extract geometric criteria
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
    HBOND_RES_DIFF = geom_criterion_values['HBOND_RES_DIFF']
    VDW_EPSILON = geom_criterion_values['VDW_EPSILON']
    VDW_RES_DIFF = geom_criterion_values['VDW_RES_DIFF']

    frame_contacts = []
    if "sb" in ITYPES:
        frame_contacts += compute_salt_bridges(traj_frag_molid, frame_idx, index_to_atom, sele1, sele2,
                                               sele1_atoms, sele2_atoms, geom_criterion_values)
    if "pc" in ITYPES:
        frame_contacts += compute_pi_cation(traj_frag_molid, frame_idx, index_to_atom, sele1, sele2,
                                            sele1_atoms, sele2_atoms, PI_CATION_CUTOFF_DISTANCE, PI_CATION_CUTOFF_ANGLE)
    if "ps" in ITYPES:
        frame_contacts += compute_pi_stacking(traj_frag_molid, frame_idx, index_to_atom, sele1, sele2,
                                              sele1_atoms, sele2_atoms, PI_STACK_CUTOFF_DISTANCE,
                                              PI_STACK_CUTOFF_ANGLE, PI_STACK_PSI_ANGLE)
    if "ts" in ITYPES:
        frame_contacts += compute_t_stacking(traj_frag_molid, frame_idx, index_to_atom, sele1, sele2,
                                             sele1_atoms, sele2_atoms, T_STACK_CUTOFF_DISTANCE, T_STACK_CUTOFF_ANGLE,
                                             T_STACK_PSI_ANGLE)
    if "vdw" in ITYPES:
        frame_contacts += compute_vanderwaals(traj_frag_molid, frame_idx, index_to_atom, sele1, sele2,
                                              VDW_EPSILON, VDW_RES_DIFF)
    if "hb" in ITYPES:
        frame_contacts += compute_hydrogen_bonds(traj_frag_molid, frame_idx, index_to_atom, solvent_resn,
                                                 sele1, sele2, sele1_atoms, sele2_atoms, None,
                                                 HBOND_CUTOFF_DISTANCE, HBOND_CUTOFF_ANGLE, HBOND_RES_DIFF)
    if "lhb" in ITYPES:
        frame_contacts += compute_hydrogen_bonds(traj_frag_molid, frame_idx, index_to_atom, solvent_resn,
                                                 sele1, sele2, sele1_atoms, sele2_atoms,
                                                 HBOND_CUTOFF_DISTANCE, HBOND_CUTOFF_ANGLE)
    if "hp" in ITYPES:
        frame_contacts += compute_hydrophobics(traj_frag_molid, frame_idx, index_to_atom, sele1, sele2,
                                               VDW_EPSILON, VDW_RES_DIFF)

    return frame_contacts


def compute_fragment_contacts(frag_idx, beg_frame, end_frame, top, traj, itypes, geom_criterion_values, stride,
                              solvent_resn, sele1, sele2, sele1_atoms, sele2_atoms, index_to_atom):
    """ 
    Reads in a single trajectory fragment and calls compute_frame_contacts on each frame

    Parameters
    ----------
    frag_idx: int
        Trajectory fragment index for worker to keep track of
    beg_frame: int
        Start frame of trajectory fragment
    end_frame: int
        Last frame of trajectory fragment
    top: str
        Topology in .pdb or .mae format
    traj: str
        Trajectory in .nc or .dcd format
    itypes: list
        Denotes the list of non-covalent interaction types to compute contacts for 
    geom_criterion_values: dict
        Dictionary containing the cutoff values for all geometric criteria
    stride: int, default = 1
        Frequency to skip frames in trajectory
    solvent_resn: string, default = TIP3
        Denotes the resname of solvent in simulation
    sele1: string, default = None
        Compute contacts on subset of atom selection based on VMD query
    sele2: string, default = None
        If second VMD query is specified, then compute contacts between atom selection 1 and 2 
    sele1_atoms: list 
        List of atom label indices for all atoms in selection 1
    sele2_atoms: list 
        List of atom label indices for all atoms in selection 2
    index_to_atom: dict
        Maps VMD atom index to Atom

    Return
    ------
    frag_idx: int 

    fragment_contacts: list of tuples, [(frame_index, atom1_label, atom2_label, itype), ...]
    """
    tic = datetime.datetime.now()
    traj_frag_molid = load_traj(top, traj, beg_frame, end_frame, stride)
    fragment_contacts = []

    # Compute contacts for each frame
    num_frag_frames = molecule.numframes(traj_frag_molid)
    for frame_idx in range(num_frag_frames):
        # if frame_idx > 1: break
        fragment_contacts += compute_frame_contacts(traj_frag_molid, frame_idx, itypes, geom_criterion_values,
                                                    solvent_resn, sele1, sele2, sele1_atoms, sele2_atoms, index_to_atom)

    # Delete trajectory fragment to clear memory
    molecule.delete(traj_frag_molid)

    # Update frame-number so it's not relative to beg_frame
    for fc in fragment_contacts:
        fc[0] = beg_frame + (fc[0] * stride)

    toc = datetime.datetime.now()
    print("Finished computing contacts for fragment %d: %d frames from %d to %d in strides of %d taking %s s" %
          (frag_idx, num_frag_frames, beg_frame, beg_frame + num_frag_frames * stride - 1, stride, (toc-tic).total_seconds()))

    return fragment_contacts


# def compute_fragment_contacts_helper(args):
#     return compute_fragment_contacts(*args)
#

def compute_contacts(top, traj, output, itypes, geom_criterion_values, cores,
                     beg, end, stride, solvent_resn, lipid_resn, sele1, sele2):
    """
    Computes non-covalent contacts across the entire trajectory and writes them to `output`.

    Parameters
    ----------
    top: Topology
        In .pdb or .mae format
    traj: Trajectory
        In .nc or .dcd format
    output: string
        Absolute path to output file
    itypes: list
        Denotes the list of non-covalent interaction types to compute contacts for 
    geom_criterion_values: dict
        Dictionary containing the cutoff values for all geometric criteria
    cores: int, default
        Number of CPU cores to parallelize over
    beg: int
        First frame to read
    end: int
        Last frame to read
    stride: int, default
        The number of frames to increment after each read frame
    solvent_resn: string
        Denotes the resname of solvent in simulation
    lipid_resn: string
        Denotes the resnames of lipids in simulation
    sele1: string, default = None
        Compute contacts on subset of atom selection based on VMD query
    sele2: string, default = None
        If second VMD query is specified, then compute contacts between atom selection 1 and 2 
    """

    contact_types = []
    for itype in itypes:
        if itype == "hb":
            contact_types += ["hbbb", "hbsb", "hbss", "wb", "wb2"]
        elif itype == "lhb":
            contact_types += ["hbls", "hblb", "lwb", "lwb2"]
        else:
            contact_types += [itype]

    index_to_atom = gen_index_to_atom(top, traj)
    sim_length = simulation_length(top, traj)
    solvent_resn = configure_solv(top, traj, solvent_resn)
    configure_lipid(lipid_resn)

    trajid = load_traj(top, traj, 0, 1, 1)
    sele1_atoms = get_selection_indices(trajid, 0, sele1)
    sele2_atoms = get_selection_indices(trajid, 0, sele2)
    molecule.delete(trajid)

    beg = max(min(beg, sim_length - 1), 0)
    end = min(max(end, beg), sim_length - 1)
    stride = max(1, stride)
    num_fragments = math.ceil((end - beg + 1) / (TRAJ_FRAG_SIZE * stride))

    # Generate input arguments for each trajectory piece
    inputqueue = Queue()
    print("Processing %s (frame %d to %d with stride of %d) as %d fragments\n" %
          (traj, beg, end, stride, num_fragments))

    for frag_idx, beg_frame in enumerate(range(beg, end + 1, TRAJ_FRAG_SIZE * stride)):
        end_frame = beg_frame + (TRAJ_FRAG_SIZE * stride) - 1
        # print(frag_idx, beg_frame, end_frame, stride)
        inputqueue.put((frag_idx, beg_frame, end_frame, top, traj, itypes, geom_criterion_values,
                        stride, solvent_resn, sele1, sele2, sele1_atoms, sele2_atoms, index_to_atom))

    # Set up result queue for workers to transfer results to the consumer
    resultsqueue = Queue()

    # Set up and start worker processes
    num_workers = max(1, cores)
    for _ in range(num_workers):
        inputqueue.put("DONE")  # Stops each worker process

    if num_workers == 1:  # Run everything in series (in addition to being slow it will consume memory)
        contact_worker(inputqueue, resultsqueue)
        resultsqueue.put("DONE")
        output_fd = open(output, "w")
        contact_consumer(resultsqueue, output_fd, itypes, beg, end, stride)
        output_fd.close()

    else:
        workers = [Process(target=contact_worker, args=(inputqueue, resultsqueue)) for _ in range(num_workers)]
        for w in workers:
            w.start()

        # Set up and start consumer process which takes contact results and saves them to output
        output_fd = open(output, "w")
        consumer = Process(target=contact_consumer, args=(resultsqueue, output_fd, itypes, beg, end, stride))
        consumer.start()

        # Wait for everyone to finish
        for w in workers:
            w.join()
        resultsqueue.put("DONE")
        consumer.join()
        output_fd.close()


def contact_worker(inputqueue, resultsqueue):
    while True:
        args = inputqueue.get()
        # print("====================== producer ====================== ")
        # print(args)
        if args == "DONE":
            return
        contacts = compute_fragment_contacts(*args)
        # print(args)
        frag_idx = args[0]
        resultsqueue.put((frag_idx, contacts))


def contact_consumer(resultsqueue, output_fd, itypes, beg, end, stride):
    import heapq

    total_frames = math.ceil((end - beg + 1) / stride)
    output_fd.write("# total_frames:%d beg:%d end:%d stride:%d interaction_types:%s\n" %
                    (total_frames, beg, end, stride, ",".join(itypes)))
    output_fd.write("# Columns: frame, interaction_type, atom_1, atom_2[, atom_3[, atom_4]]\n")

    resultheap = []
    waiting_for_fragment = 0

    while True:
        result = resultsqueue.get()
        if result == "DONE":
            assert not resultheap, "resultheap not empty .. unwritten contacts on termination"
            output_fd.close()
            break

        fragment = result[0]
        contacts = result[1]
        heapq.heappush(resultheap, (fragment, contacts))

        # If this is the contact frame we're waiting for, go ahead and write it
        while resultheap and waiting_for_fragment == resultheap[0][0]:
            fragment, contacts = heapq.heappop(resultheap)
            waiting_for_fragment = fragment + 1

            # Strip VMD id, order atom 1 and 2
            for interaction in contacts:
                for a in range(2, len(interaction)):
                    atom_str = interaction[a]
                    interaction[a] = atom_str[0:atom_str.rfind(":")]

                # Sort atom 1 and 2 lexicographically
                if interaction[3] < interaction[2]:
                    interaction[2], interaction[3] = interaction[3], interaction[2]  # Swap

            contact_line_hash = set()
            for interaction in contacts:
                # Write to file
                # waiting_for_frame = max(waiting_for_frame, interaction[0] + stride)
                interaction[0] = str(interaction[0])
                contact_line = "\t".join(interaction)
                if contact_line not in contact_line_hash:
                    output_fd.write(contact_line)
                    output_fd.write("\n")
                    contact_line_hash.add(contact_line)
