__author__ = 'Rasmus Fonseca <fonseca.rasmus@gmail.com>'
__license__ = "APACHE2"

__all__ = ['parse_contacts', 'parse_residuelabels', 'res_contacts']


def parse_contacts(input_lines, itypes=None):
    """
    Read a contact-file (tab-separated file with columns: frame, i-type, atomid1, atomid2[, atomid3[, atomid4]] and
    return it as a list of lists with frames converted to ints. The total number of frames is also returned.

    Example
    -------
        parse_contacts([
            "# total_frames:2\n",
            "0  hbbb    A:ALA:1:N   A:THR:10:O\n",
            "0  vdw     A:ALA:1:CB  B:CYS:3:H\n",
            "1  vdw     A:ALA:1:N   A:THR:10:C\n"
        ])
        # returns:
        # ([
        #        [0, "hbbb", "A:ALA:1:N", "A:THR:10:O"],
        #        [0, "vdw", "A:ALA:1:CB", "B:CYS:3:H"],
        #        [1, "vdw", "A:ALA:1:N", "A:THR:10:C"]
        #  ], 2)

    Parameters
    ----------
    input_lines: iterable
        Iterator of over a set of strings. Can be a file-handle

    itypes: set of str | None
        Interactions to include in the output

    Returns
    -------
    (list of list, int)
        The list of interactions and the total number of frames
    """
    ret = []
    total_frames = 0
    for line in input_lines:
        line = line.strip()
        if "total_frames" in line:
            tokens = line.split(" ")
            total_frames = int(tokens[1][tokens[1].find(":")+1:])

        if len(line) == 0 or line[0] == "#":
            continue

        tokens = line.split("\t")
        tokens[0] = int(tokens[0])

        if itypes is None or tokens[1] in itypes:
            ret.append(tokens)

    return ret, total_frames


def parse_residuelabels(label_file):
    """
    Parses a residue-label file and generates a dictionary mapping residue identifiers (e.g. A:ARG:123) to a
    user-specified label, trees that can be parsed by flareplots, and a color indicator for vertices.

    Parameters
    ----------
    label_file : file
        A flare-label file where each line contains 2-3 columns formatted as
          - CHAIN:RESN:RESI (e.g. A:ARG:123)
          - [[TOPLEVEL.]MIDLEVEL.]LABEL (e.g. Receptor.Helix2.2x44)
          - COLOR (e.g. #FF0000 or white)

    Returns
    -------
    dict of str : (dict of str : str)
        Keys are all residue identifiers and values are dicts that hold both the LABEL by itself (key "label", the full
        tree-path (key "treepath") and a CSS-compatible color string (key "color").

    Raises
    ------
    AssertionError
        if a residue identifier (CHAIN:RESN:RESI) is specified twice in the file, or if a LABEL appears twice.
    """
    if label_file is None:
        return None

    ret = {}
    flarelabels = set()  # Only used to check for duplicates
    for line in label_file:
        line = line.strip()
        if not line:
            continue  # Ignore empty lines

        columns = line.split("\t")
        residentifier = columns[0]
        flaretreepath = columns[1] if len(columns) > 1 else columns[0]
        flarelabel = flaretreepath.split(".")[-1]
        flarecolor = columns[2] if len(columns) > 2 else "white"
        if residentifier in ret:
            raise AssertionError("Residue identifier '"+residentifier+"' appears twice in "+label_file.name)
        if flarelabel in flarelabels:
            raise AssertionError("Flare label '"+flarelabel+"' used twice in "+label_file.name)

        ret[residentifier] = {"label": flarelabel, "treepath": flaretreepath, "color": flarecolor}
        flarelabels.add(flarelabel)

    return ret


def res_contacts(contacts):
    """
    Convert atomic contacts into unique residue contacts.

    Example
    -------
        res_frequencies([
            [0, 'hbbb', 'A:ASN:108:O', 'A:ARG:110:N'],
            [0, 'vdw', 'A:ARG:110:N', 'A:ASN:108:CB']
        ])
        # Returns: [[0, 'A:ASN:108', 'A:ARG:110']]

    Parameters
    ----------
    contacts: List of list
        List of contacts, where each contact is a list of frame-num, i-type, and atom names

    Returns
    -------
    List of list
        Each entry is a list with a frame and two residue identifiers
    """
    from collections import defaultdict
    # Associates a frame-number with a set of contacts
    frame_dict = defaultdict(set)

    for atom_contact in contacts:
        frame = atom_contact[0]
        resi1 = ":".join(atom_contact[2].split(":")[0:3])
        resi2 = ":".join(atom_contact[3].split(":")[0:3])
        if resi1 < resi2:
            resi1, resi2 = resi2, resi1
        frame_dict[frame].add((resi1, resi2))

    ret = []
    for frame in sorted(frame_dict):
        for resi1, resi2 in frame_dict[frame]:
            ret.append([frame, resi1, resi2])

    return ret





