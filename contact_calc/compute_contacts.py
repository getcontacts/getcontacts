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
try:
    from vmd import *  # Loads the static `molecule` object
except ModuleNotFoundError:
    import sys
    sys.stderr.write("Error: The 'vmd-python' dependency could not be located. Please follow instructions "
                     "at https://github.com/getcontacts/getcontacts to install.\n")
    sys.exit(1)


from .contact_utils import *
from .aromatics import *
from .hbonds import *
from .salt_bridges import *
from .pi_cation import *
from .vanderwaals import *
from .hydrophobics import *

##############################################################################
# Global Variables
##############################################################################
TRAJ_FRAG_SIZE = 100

##############################################################################
# Functions
##############################################################################


def compute_frame_contacts(molid, frame, itypes, geom_criteria, sele1, sele2,
                           sele1_atoms, sele2_atoms, index_to_atom, ligand_anions, ligand_cations):
    """
    Computes each of the specified non-covalent interaction type for a single frame

    Parameters
    ----------
    molid: int
        Identifier to simulation fragment in VMD
    frame: int
        Frame number to query
    itypes: list
        Denotes the list of non-covalent interaction types to compute contacts for
    geom_criteria: dict
        Dictionary containing the cutoff values for all geometric criteria
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

    Returns
    -------
    frame_contacts: list of tuples, [(frame_index, atom1_label, atom2_label, itype), ...]

    """

    frame_contacts = []
    if "sb" in itypes:
        frame_contacts += compute_salt_bridges(molid, frame, index_to_atom, sele1, sele2, geom_criteria, ligand_anions, ligand_cations)
    if "pc" in itypes:
        frame_contacts += compute_pi_cation(molid, frame, index_to_atom, sele1, sele2, geom_criteria)
    if "ps" in itypes:
        frame_contacts += compute_pi_stacking(molid, frame, index_to_atom, sele1, sele2,
                                              sele1_atoms, sele2_atoms, geom_criteria)
    if "ts" in itypes:
        frame_contacts += compute_t_stacking(molid, frame, index_to_atom, sele1, sele2,
                                             sele1_atoms, sele2_atoms, geom_criteria)
    if "vdw" in itypes:
        frame_contacts += compute_vanderwaals(molid, frame, index_to_atom, sele1, sele2, geom_criteria)
    if "hb" in itypes:
        frame_contacts += compute_hydrogen_bonds(molid, frame, index_to_atom,
                                                 sele1, sele2, sele1_atoms, sele2_atoms, geom_criteria)
    if "hp" in itypes:
        frame_contacts += compute_hydrophobics(molid, frame, index_to_atom, sele1, sele2, geom_criteria)

    return frame_contacts


def compute_fragment_contacts(frag_idx, beg_frame, end_frame, top, traj, itypes, geom_criterion_values, stride, distout,
                              sele1, sele2, sele1_atoms, sele2_atoms, index_to_atom, ligand_anions, ligand_cations):
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
    distout: bool
        Whether to write out distances or not
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
    molid = load_traj(top, traj, beg_frame, end_frame, stride)
    fragment_contacts = []

    # Compute contacts for each frame
    num_frag_frames = molecule.numframes(molid)
    for frame_idx in range(num_frag_frames):
        # if frame_idx > 1: break
        fragment_contacts += compute_frame_contacts(molid, frame_idx, itypes, geom_criterion_values,
                                                    sele1, sele2, sele1_atoms, sele2_atoms, index_to_atom, ligand_anions, ligand_cations)

    if distout:
        for contact in fragment_contacts:
            frame = contact[0]
            a1_id = contact[2].split(":")[-1]
            a2_id = contact[3].split(":")[-1]
            distance = compute_distance(molid, frame, a1_id, a2_id)
            contact.append("%.3f" % distance)

    # Delete trajectory fragment to clear memory
    molecule.delete(molid)

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
                     beg, end, stride, distout, ligand_sele, solv_sele, lipid_sele, sele1, sele2):
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
    stride: int
        The number of frames to increment after each read frame
    distout: bool
        Indicates whether to write out distances
    ligand_sele: string
        Ligand VMD selection
    solv_sele: string
        Solvent VMD selection
    lipid_sele: string
        Lipid VMD selection
    sele1: string, default = None
        Compute contacts on subset of atom selection based on VMD query
    sele2: string, default = None
        If second VMD query is specified, then compute contacts between atom selection 1 and 2 
    """

    index_to_atom = gen_index_to_atom(top, traj)
    sim_length = simulation_length(top, traj)
    configure_solv(top, traj, solv_sele)
    configure_lipid(top, traj, lipid_sele)
    configure_ligand(top, traj, ligand_sele, sele1, sele2)
    ligand_anions, ligand_cations = extract_ligand_features(top, traj, index_to_atom)

    trajid = load_traj(top, traj, 0, 1, 1)
    sele1_atoms = get_selection_indices(trajid, 0, sele1)
    sele2_atoms = get_selection_indices(trajid, 0, sele2)
    molecule.delete(trajid)

    if len(sele1_atoms) == 0:
        print("The --sele selection is empty!")
        sys.exit(1)
    if len(sele2_atoms) == 0:
        print("The --sele2 selection is empty!")
        sys.exit(1)

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
                        stride, distout, sele1, sele2, sele1_atoms, sele2_atoms, index_to_atom, ligand_anions, ligand_cations))

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
        contact_consumer(resultsqueue, output_fd, itypes, beg, end, stride, distout)
        output_fd.close()

    else:
        workers = [Process(target=contact_worker, args=(inputqueue, resultsqueue)) for _ in range(num_workers)]
        for w in workers:
            w.start()

        # Set up and start consumer process which takes contact results and saves them to output
        output_fd = open(output, "w")
        consumer = Process(target=contact_consumer, args=(resultsqueue, output_fd, itypes, beg, end, stride, distout))
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
        if args == "DONE":
            return
        contacts = compute_fragment_contacts(*args)
        frag_idx = args[0]
        resultsqueue.put((frag_idx, contacts))


def contact_consumer(resultsqueue, output_fd, itypes, beg, end, stride, distout):
    import heapq

    total_frames = math.ceil((end - beg + 1) / stride)
    output_fd.write("# total_frames:%d beg:%d end:%d stride:%d interaction_types:%s\n" %
                    (total_frames, beg, end, stride, ",".join(itypes)))

    if distout:
        output_fd.write("# Columns: frame, interaction_type, atom_1, atom_2[, atom_3[, atom_4]], distance_1-2\n")
    else:
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
                for a in range(2, len(interaction) - 1 if distout else len(interaction)):
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
