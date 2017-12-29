from __future__ import division
from collections import defaultdict
import json


def get_residue_from_atom(atom_name):
    split_atom_name = atom_name.split(':')
    return "%s:%s" % (split_atom_name[1], split_atom_name[2])


def create_atom_file(itype, output_dir, stitched_filename, sim_length, atompair_to_frames):
    # Assemble and write out .json
    edges = []
    for atom_pair in atompair_to_frames:  # Loop over all respairs
        atom1, atom2 = atom_pair.split()
        edge = {}
        edge["name1"] = atom1
        edge["name2"] = atom2
        edge["frames"] = list(atompair_to_frames[atom_pair])
        edges.append(edge)
    with open(output_dir + '/' + itype + "_atom.json", 'w+') as fq:
        fq.write(json.dumps({"edges":edges}))


def create_frequencies_file(itype, output_dir, stitched_filename, sim_length, respair_to_frames):
    # Write out .txt file
    with open(output_dir + '/' + itype + "_frequencies.txt", 'w+') as fq:
        fq.write('NumFrames:\t%d\n' % (sim_length))
        for interaction in respair_to_frames:
            interaction_frequency = len(respair_to_frames[interaction])/sim_length
            fq.write("%s\t%.4f\t%d\n" % (interaction, interaction_frequency, len(respair_to_frames[interaction])))

    # Assemble and write out .json
    edges = []
    for res_pair in respair_to_frames:  # Loop over all respairs
        res1, res2 = res_pair.split()
        edge = {}
        edge["name1"] = res1
        edge["name2"] = res2
        edge["frames"] = list(respair_to_frames[res_pair])
        edges.append(edge)
    with open(output_dir + '/' + itype + "_res.json", 'w+') as fq:
        fq.write(json.dumps({"edges":edges}))


def get_atompair_set(itype, stitched_lines):
    atompair_to_frames = defaultdict(set)
    for stitched_line in stitched_lines:
        frame = int(stitched_line[0])
        atom1 = stitched_line[1]
        atom2 = stitched_line[-2]
        atom_pair = '\t'.join(sorted([atom1, atom2]))
        atompair_to_frames[atom_pair].add(frame)
    return atompair_to_frames


def get_respair_set(itype, stitched_lines):
    """Creates a dict from residue_pair to all of the frames in which that residue_pair interaction appears"""
    respair_to_frames = defaultdict(set)
    for stitched_line in stitched_lines:
        frame = int(stitched_line[0])
        atom1 = stitched_line[1]
        atom2 = stitched_line[-2]

        res_pair = '\t'.join(sorted([get_residue_from_atom(atom1), get_residue_from_atom(atom2)]))
        respair_to_frames[res_pair].add(frame)
    return respair_to_frames


def make_frequencies_file(itype, output_dir, stitched_filename, sim_length):
    with open(stitched_filename) as stitched_open:
        stitched_lines = [line.strip().split() for line in stitched_open.readlines()]
    respair_to_frames = get_respair_set(itype, stitched_lines)
    atompair_to_frames = get_atompair_set(itype, stitched_lines)
    create_frequencies_file(itype, output_dir, stitched_filename, sim_length, respair_to_frames)
    create_atom_file(itype, output_dir, stitched_filename, sim_length, atompair_to_frames)

