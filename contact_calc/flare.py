"""
Functions to transform and write flare objects. A "flare" is a json object / dictionary with the following properties
 * It must have a "edges" key where the value is a list of edge-objects. Each edge-object has the keys
   - "name1": The identifier of an end-vertex
   - "name2": The identifier of the other end-vertex
   - "frames": a list of non-negative integers indicating which time-points or conditions the edge is active
   - "color": The edge color (not mandatory, default: gray)
   - "weight": The edge weight (not mandatory, default: 1)
 * It can have a "trees" key where the value is a list of tree-objects, each with the following keys:
   - "treeLabel": The name of the tree. If 'default' it will be the default tree when shown on the flareplot page
   - "treeProperties": A list of objects with "path", and "key" keys used to identify the tree-path of vertices
 * It can have a "tracks" key where the value is a list of tree-objects, each with the following keys:
   - "trackLabel": The name of the track. If 'default' it will be the default track when shown on the flareplot page
   - "trackProperties": A list of objects with "nodeName", "color", and "size" keys indicating the track-objects for
     each node.
 * It can have a "frameDict" key (TODO: expand)
 * It can have a "defaults" key (TODO: expand)

There are three subtypes of flare objects:
 * Single-flare: has no "frameDict" key and all "frames"-lists have exactly one entry
 * Time-flare: has no "frameDict" key
 * Compare-flare: has a "frameDict" key and the frames all index an entry in it
(TODO: Expand a bit with examples)
"""

import json
import re
import sys
from contact_calc.transformations import res_contacts


def write_json(flare, fstream):
    """
    Serialize the flare object as json with two space indentation and write to fname. The output is expected to have a
    number of integer-lists (frames) that will be compressed into single lines.

    Parameters
    ----------
    flare: dict of (str, list)
        Flare object to write

    fstream: str | TextIOWrapper
        Filename or stream to write
    """
    pretty_json = json.dumps(flare, indent=2)
    # Put list-of-numbers on a single line
    pretty_json = re.sub(r"(?<=\d,)\n *|(?<=\[)\n *(?=\d)|(?<=\d)\n *(?=\])", "", pretty_json, flags=re.MULTILINE)

    if type(fstream) == str:
        with open(fstream, "w") as f:
            f.write(pretty_json)
        return True
    elif fstream is not None:
        fstream.write(pretty_json)
        return True
    else:
        print(pretty_json)
        return False


def is_single_flare(flare):
    if "frameDict" in flare:
        return False
    if any([len(e["frames"] != 1 or e["frames"][0] != 0 for e in flare["edges"])]):
        return False
    return True


def is_time_flare(flare):
    if "frameDict" in flare:
        return False
    return True


def is_compare_flare(flare):
    if "frameDict" not in flare:
        return False

    # Check that all frames correspond to a condition in the frameDict
    conditions = range(len(flare["frameDict"]))
    for edge in flare["edges"]:
        if any([f not in conditions for f in edge["frames"]]):
            return False

    return True


def create_flare(contacts, resi_labels):
    """
    Creates a flare from a contact-list and residue labels. If `resi_labels` isn't `None` then the "trees" and "tracks"
    will be generated as well.

    Parameters
    ----------
    contacts : list of list
        Each entry specifies a frame-number, an interaction type, and 2 to 4 atom-tuples depending on the interaction
        type. Water mediated and water-water mediated interactions will have waters in the third and fourth tuples.
    #contacts : list of tuples of (str, str, tuple, tuple [[, tuple], tuple])
    #    Each entry specifies a frame-number, an interaction type, and 2 to 4 atom-tuples depending on the interaction
    #    type. Water mediated and water-water mediated interactions will have waters in the third and fourth tuples.

    resi_labels : dict of (str : dict of (str : str))
        Each key is a residue identifier and the associated value is a dictionary with the label, tree-path, and color
        that are consistent with the format of flareplots.

    Returns
    -------
    dict of str : list
        The types of the list contents varies depending on the key, but the format corresponds to the specifications of
        jsons used as input for flareplot. For example:
            {
              "edges": [
                {"name1": "ARG1", "name2": "ALA3", "frames": [0,4,10]},
                {"name1": "ALA3", "name2": "THR2", "frames": [1,2]}
              ],
              "trees": [
                {"treeName": "DefaultTree", "treePaths": ["Group1.ARG1", "Group1.THR2", "Group2.ALA3", "Group2.CYS4"]}
              ],
              "tracks": [
                {"trackName": "DefaultTrack", "trackProperties": [
                  {"nodeName": "ARG1", "size": 1.0, "color": "red"},
                  {"nodeName": "THR2", "size": 1.0, "color": "red"},
                  {"nodeName": "ALA3", "size": 1.0, "color": "blue"},
                  {"nodeName": "CYS4", "size": 1.0, "color": "blue"}
                ]}
              ]
            }
    """
    ret = {
        "edges": []
    }

    # Strip atom3, atom4, and atom names
    # unique_chains = set([c[2][0] for c in contacts] + [c[3][0] for c in contacts])
    # contacts = [(c[0], c[1], c[2][0:3], c[3][0:3]) for c in contacts]
    rcontacts = res_contacts(contacts)

    resi_edges = {}
    for contact in rcontacts:
        # Compose a key for atom1 and atom2 that ignores the order of residues
        a1_key = contact[1]
        a2_key = contact[2]
        # a1_key = ":".join(contact[2][0:3])
        # a2_key = ":".join(contact[3][0:3])
        # if a1_key == a2_key:
        #     continue
        # if a1_key > a2_key:
        #     a1_key, a2_key = a2_key, a1_key
        contact_key = a1_key + a2_key

        # Look up labels
        if resi_labels:
            if a1_key not in resi_labels or a2_key not in resi_labels:
                print("create_flare: Omitting contact "+str(contact)+" as it doesn't appear in flare-label file")
                continue
            a1_label = resi_labels[a1_key]["label"]
            a2_label = resi_labels[a2_key]["label"]
        else:
            a1_label = a1_key
            a2_label = a2_key

        # Create contact_key if it doesn't exist
        if contact_key not in resi_edges:
            edge = {"name1": a1_label, "name2": a2_label, "frames": []}
            resi_edges[contact_key] = edge
            ret["edges"].append(edge)

        resi_edges[contact_key]["frames"].append(int(contact[0]))

    # Sort edge frames and ensure that there are no duplicates
    for e in ret["edges"]:
        e["frames"] = sorted(set(e["frames"]))

    # Create "trees" and "tracks" sections if resi_labels specified
    if resi_labels is not None:
        tree = {"treeLabel": "DefaultTree", "treePaths": []}
        ret["trees"] = [tree]

        track = {"trackLabel": "DefaultTrack", "trackProperties": []}
        ret["tracks"] = [track]

        for rlabel in resi_labels.values():
            tree["treePaths"].append(rlabel["treepath"])
            track["trackProperties"].append({
                "nodeName": rlabel["label"],
                "color": rlabel["color"],
                "size": 1.0
            })

    return ret


def compose_frequencytable(freq_table, column_headers, freq_cutoff):
    """
    asdf
    Parameters
    ----------
    freq_table: dict of tuple: array
        Each key is a tuple of residue ids and the value is a numpy array indicating frequencies for each column
        Example:
            {
                ('A:ARG:293', 'A:GLN:297'): array([ 1.,  1.]),
                ('A:ASN:280', 'A:TRP:246'): array([ 1.,  1.]),
                ('A:ADN:400', 'A:ASN:181'): array([ 1.,  0.]),
                ('A:ARG:293', 'A:ASN:280'): array([ 0.,  1.]),
                ('A:ASN:181', 'A:TRP:246'): array([ 0.,  1.]),
                ('A:ADN:400', 'A:ARG:293'): array([ 0.,  1.]),
                ('A:GLN:297', 'A:GLN:297'): array([ 0.,  1.])
            }

    column_headers: list of str
        Name of each condition / column

    freq_cutoff: float
        Frequency values above this will be considered an edge
    """
    ret = {
        "edges": [],
        "frameDict": {str(i): col_head for i, col_head in enumerate(column_headers)}
    }

    for (res1, res2), freq_arr in freq_table.items():
        edge = {"name1": res1, "name2": res2, "frames": []}
        for i in range(freq_arr.shape[0]):
            if freq_arr[i] >= freq_cutoff:
                edge["frames"].append(i)

        ret["edges"].append(edge)

    return ret


def compose_flares(singleflares, names):
    """
    Take a list of single-flare files and convert them to a multi-flare file.
    Briefly, a single-flare is a flareplot input file where the value of the `frames` attribute of each edge is
    `[0]`. This corresponds e.g. to the contact network of a single protein
    structure. A multi-frame has edges with multiple values in the `frames` list,
    but each index has a label specified in a `frameDict` entry. Multi-frames are
    used to generate finger-print comparison flareplots.

    Parameters
    ----------
    singleflares : list of dict
        Single-flare objects representing graph in a certain condition

    names: list of str
        Names to assign each flare

    Returns
    -------
    dict
        A multiflare object where the entries of the 'frameDict' correspond to the indicated names
    """
    assert all([is_single_flare(f) for f in singleflares])

    ret = {
        "edges": [],
        "frameDict": {}
    }

    # Compose edges.
    def findedge(edge):
        node_names = set([edge["name1"], edge["name2"]])
        for e in ret["edges"]:
            if e["name1"] in node_names and e["name2"] in node_names:
                return e
        return None

    for flareidx, flare in enumerate(singleflares):
        for edge in flare["edges"]:
            existing_edge = findedge(edge)
            if existing_edge is None:
                existing_edge = {"name1": edge["name1"],
                                 "name2": edge["name2"],
                                 "frames": [],
                                 "colors": [],
                                 "widths": []}
                ret["edges"].append(existing_edge)

            existing_edge["frames"].append(flareidx)
            if "color" in edge:
                existing_edge["colors"].append(edge["color"])
            if "width" in edge:
                existing_edge["widths"].append(float(edge["width"]))

    for edge in ret["edges"]:
        if edge["widths"]:
            width_sum = sum(edge["widths"]) + (len(singleflares) - len(edge["widths"])) * 1.0
            edge["width"] = width_sum / len(singleflares)
        del edge["widths"]

        if edge["colors"]:
            edge["color"] = edge["colors"][0]  # TODO: Average colors
        del edge["colors"]

    # Compose trees - currently it just takes the trees of the first flare
    if any(map(lambda f: "trees" in f and len(f["trees"]) > 1, singleflares)):
        print("Error, can't compose flares that have multiple trees")
        # TODO: Print which flares have multiple trees
        sys.exit(1)

    if any(map(lambda f: "trees" in f and len(f["trees"]) == 1, singleflares)):
        ret["trees"] = [{"treeLabel": "DefaultTree", "treePaths": []}]

        def findpath(p):
            p_label = p[p.rfind(".")+1:]
            for tp in ret["trees"][0]["treePaths"]:
                tp_label = tp[tp.rfind(".")+1:]
                if tp_label == p_label:
                    return tp
            return None

        for flare in singleflares:
            if "trees" in flare and len(flare["trees"]) > 0:
                for treepath in flare["trees"][0]["treePaths"]:
                    existing_path = findpath(treepath)
                    if existing_path is None:
                        ret["trees"][0]["treePaths"].append(treepath)
                    elif existing_path != treepath:
                        fidx = singleflares.index(flare)
                        print("Can't compose conflicting tree-paths:")
                        print("> "+existing_path)
                        print("> "+treepath)
                        print("> "+names[fidx])
                        exit(1)

    # Compose tracks - currently adds all tracks
    # TODO: Check for identical tracks and only add one
    if any(map(lambda f: "tracks" in f and len(f["tracks"]) > 0, singleflares)):
        ret["tracks"] = []
        for flare, flarename in zip(singleflares, names):
            if "tracks" in flare and len(flare["tracks"]) > 0:
                for track in flare["tracks"]:
                    ret["tracks"].append({"trackLabel": track["trackLabel"]+"_"+flarename,
                                          "trackProperties": track["trackProperties"]})

    # Compose frameDict
    for frameidx, name in enumerate(names):
        ret["frameDict"][str(frameidx)] = name

    return ret


__author__ = 'Rasmus Fonseca <fonseca.rasmus@gmail.com>'
__license__ = "Apache License 2.0"
