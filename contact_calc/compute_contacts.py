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
from multiprocessing import *
from vmd import *  # Loads the static `molecule` object

from .contact_utils import *
from .aromatics import *
from .hbonds import *
from .salt_bridges import *
from .pi_cation import *
from .vanderwaals import *

##############################################################################
# Global Variables
##############################################################################
TRAJ_FRAG_SIZE = 100
full_name_dirs = {'hbbb': 'hydrogen_bonds/backbone_backbone_hydrogen_bonds',
                  'hbsb': 'hydrogen_bonds/sidechain_backbone_hydrogen_bonds',
                  'hbss': 'hydrogen_bonds/sidechain_sidechain_hydrogen_bonds',
                  'vdw': 'van_der_Waals',
                  'sb': 'salt_bridges',
                  'ts': 't_stacking',
                  'ps': 'pi_stacking',
                  'pc': 'pi_cation',
                  'wb': 'hydrogen_bonds/water_mediated_hydrogen_bonds',
                  'wb2': 'hydrogen_bonds/extended_water_mediated_hydrogen_bonds',
                  'lhb': 'ligand_hydrogen_bonds/hydrogen_bonds',
                  'hlb': 'ligand_hydrogen_bonds/backbone_hydrogen_bonds',
                  'hls': 'ligand_hydrogen_bonds/sidechain_hydrogen_bonds',
                  'lwb': 'ligand_hydrogen_bonds/water_mediated_hydrogen_bonds',
                  'lwb2': 'ligand_hydrogen_bonds/extended_water_mediated_hydrogen_bonds',
                  }

##############################################################################
# Functions
##############################################################################

# Note: Trying to write directly to output and doing postprocessing after
# May save a lot of memory instead of having massive arrays, but there is
# also the IO time (which would happen anyways). Figure out best output
# format and most efficient way to write to disk.


def compute_frame_contacts(traj_frag_molid, frag_idx, frame_idx, ITYPES, geom_criterion_values, solvent_resn, sele_id,
                           ligand, index_to_label):
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

    # Extract geometric criterion
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
    if "sb" in ITYPES:
        frame_contacts += compute_salt_bridges(traj_frag_molid, frame_idx, sele_id, SALT_BRIDGE_CUTOFF_DISTANCE)
    if "pc" in ITYPES:
        frame_contacts += compute_pi_cation(traj_frag_molid, frame_idx, index_to_label, sele_id, PI_CATION_CUTOFF_DISTANCE, PI_CATION_CUTOFF_ANGLE)
    if "ps" in ITYPES:
        frame_contacts += compute_pi_stacking(traj_frag_molid, frame_idx, index_to_label, sele_id, PI_STACK_CUTOFF_DISTANCE, PI_STACK_CUTOFF_ANGLE, PI_STACK_PSI_ANGLE)
    if "ts" in ITYPES:
        frame_contacts += compute_t_stacking(traj_frag_molid, frame_idx, index_to_label, sele_id, T_STACK_CUTOFF_DISTANCE, T_STACK_CUTOFF_ANGLE, T_STACK_PSI_ANGLE)
    if "vdw" in ITYPES:
        frame_contacts += compute_vanderwaals(traj_frag_molid, frame_idx, index_to_label, sele_id, VDW_EPSILON)
    if "hb" in ITYPES:
        frame_contacts += compute_hydrogen_bonds(traj_frag_molid, frame_idx, index_to_label, solvent_resn, sele_id, None, HBOND_CUTOFF_DISTANCE, HBOND_CUTOFF_ANGLE)
    if "lhb" in ITYPES:
        frame_contacts += compute_hydrogen_bonds(traj_frag_molid, frame_idx, index_to_label, solvent_resn, sele_id, ligand, HBOND_CUTOFF_DISTANCE, HBOND_CUTOFF_ANGLE)

    toc = datetime.datetime.now()
    print("Finished computing contacts for frag %s frame %s in %s s" % (frag_idx, frame_idx, (toc-tic).total_seconds()))
    return frame_contacts


def compute_fragment_contacts(frag_idx, beg_frame, end_frame, top, traj, output, itypes, geom_criterion_values, stride, solvent_resn, sele_id, ligand, index_to_label):
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
    top: str
        Topology in .pdb or .mae format
    traj: str
        Trajectory in .nc or .dcd format
    output: str
        Path to file where results should be written
    itypes: list
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
    traj_frag_molid = load_traj(top, traj, beg_frame, end_frame, stride)
    fragment_contacts = []

    # Compute contacts for each frame
    num_frag_frames = molecule.numframes(traj_frag_molid)
    for frame_idx in range(num_frag_frames):
        # if frame_idx > 1: break
        fragment_contacts += compute_frame_contacts(traj_frag_molid, frag_idx, frame_idx, itypes, geom_criterion_values, solvent_resn, sele_id, ligand, index_to_label)

    # Delete trajectory fragment to clear memory
    molecule.delete(traj_frag_molid)

    # Update frame-number so it's not relative to beg_frame
    for fc in fragment_contacts:
        fc[0] = beg_frame + fc[0]

    # Write directly out to temporary output
    # print("Writing output to seperate files, one for each itype ...")
    #
    # fd_map = {itype: open(OUTPUT + "_" + itype + "_frag_" + str(frag_idx) + ".txt", 'w') for itype in contact_types}
    # for contact in fragment_contacts:
    #     itype_key = contact[-1]
    #     output_string = str(frag_idx) + "\t" + "\t".join(map(str, contact)) + "\n"
    #     fd_map[itype_key].write(output_string)
    #
    # for itype in fd_map:
    #     fd_map[itype].close()
    #
    # return frag_idx, num_frag_frames - 1

    return fragment_contacts


def compute_fragment_contacts_helper(args):
    return compute_fragment_contacts(*args)


# def stitch_fragment_contacts(itype, OUTPUT_DIR, frag_contact_files, frag_idx_to_length):
#     """
#     Stitch together multiple fragments of non-covalent contacts into
#     single file and delete individual fragment files.
#
#     Parameters
#     ----------
#     itype: Type of non-covalent contact
#     OUTPUT_DIR: string
#         Absolute path to output directory
#     frag_contact_files: list of strings
#         List of paths to fragment contact files
#     frag_idx_to_length: dict from int to int
#         Map the fragment index to length of fragment
#     """
#     print("Stitching %s ..." % (itype))
#     # stitched_filename = OUTPUT_DIR + "/" + itype + ".txt"
#     # stitched_filename = OUTPUT_DIR + full_name_dirs[contact_type] + '/' + 'raw_frame_output.txt'
#     stitched_filename = OUTPUT_DIR + full_name_dirs[itype] + '/' + itype + '.txt'
#     fo = open(stitched_filename, 'w')
#
#     num_frames = 0
#     for frag_contact_file in frag_contact_files:
#         frag_idx = int(frag_contact_file.split("/")[-1].strip(".txt").split("_")[2])
#         ffrag = open(frag_contact_file, 'r')
#         for line in ffrag:
#             linfo = line.split("\t")
#             frame_idx = int(linfo[1])
#             new_frame_idx = num_frames + frame_idx - 1
#             output_string = str(new_frame_idx) + "\t" + "\t".join(linfo[2:])
#             fo.write(output_string)
#         num_frames += frag_idx_to_length[frag_idx]
#         ffrag.close()
#         os.remove(frag_contact_file)
#
#     fo.close()
#
#     return stitched_filename


def compute_contacts(top, traj, output, itypes, geom_criterion_values, cores, stride, solvent_resn, sele_id, ligand):
    """ Computes non-covalent contacts across the entire trajectory and writes them to `output`.

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
    cores: int, default = 6
        Number of CPU cores to parallelize over
    stride: int, default = 1
        Frequency to skip frames in trajectory
    solvent_resn: string, default = TIP3
        Denotes the resname of solvent in simulation
    sele_id: string, default = None
        Compute contacts on subset of atom selection based on VMD query
    ligand: string, default = None
        Include ligand resname if computing contacts between ligand and binding pocket residues
    """

    contact_types = []
    for itype in itypes:
        if itype == "hb":
            contact_types += ["hbbb", "hbsb", "hbss", "wb", "wb2"]
        elif itype == "lhb":
            contact_types += ["hls", "hlb", "lwb", "lwb2"]
        else:
            contact_types += [itype]

    index_to_label = gen_index_to_atom_label(top, traj)
    sim_length = simulation_length(top, traj)
    input_args = []

    # Generate input arguments for each trajectory piece
    print("Processing %s with %s total frames and stride %s" % (traj, str(sim_length), str(stride)))
    for frag_idx, beg_frame in enumerate(range(0, sim_length, TRAJ_FRAG_SIZE)):
        # if frag_idx > 0: break
        end_frame = beg_frame + TRAJ_FRAG_SIZE - 1
        print("Preparing fragment %s, beg_frame:%s end_frame:%s" % (frag_idx, beg_frame, end_frame))
        input_args.append((frag_idx, beg_frame, end_frame, top, traj, output, itypes, geom_criterion_values,
                           stride, solvent_resn, sele_id, ligand, index_to_label))

    # Parallel computation
    pool = Pool(processes=cores)
    contacts = pool.map(compute_fragment_contacts_helper, input_args)
    pool.close()
    pool.join()
    contacts = [x for y in contacts for x in y]  # Flatten

    # Serial computation: Use this mode to debug since multiprocessing module doesn't trace back to bugs. 
    # contacts = compute_fragment_contacts_helper(input_args[0])

    # Sort and write to output-file
    contacts.sort(key=lambda i: i[0])
    with open(output, "w") as output_fd:
        output_fd.write("# total_frames:%d interaction_types:%s\n" % (sim_length, ",".join(itypes)))
        output_fd.write("# Columns: frame, interaction_type, atom_1, atom_2[, atom_3[, atom_4]]\n")
        for interaction in contacts:
            # Strip vmd ID from atom strings
            for a in range(2, len(interaction)):
                atom_str = interaction[a]
                interaction[a] = atom_str[0:atom_str.rfind(":")]

            # Write to file
            output_fd.write("\t".join(map(str, interaction)))
            output_fd.write("\n")

