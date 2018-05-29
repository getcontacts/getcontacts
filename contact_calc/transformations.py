

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
    Parses a residue-label file and generates two dictionaries mapping residue identifiers (e.g. A:ARG:123) to a
    user-specified label, trees that can be parsed by flareplots, and a color indicator for vertices respectively.

    Example
    -------
        parse_residuelabels([
            "A:ALA:4    A4  white\n",
            "A:CYS:5    C5  yellow\n",
            "A:TRP:6    W6\n"
        ])
        # => ({"A:ALA:4": "A4", "A:CYS:5": "C5", "A:TRP:6": "W6" }, { "A:ALA:4": "white", "A:CYS:5": "yellow"})

    Parameters
    ----------
    label_file : file
        A flare-label file where each line contains 2-3 columns formatted as
          - CHAIN:RESN:RESI (e.g. A:ARG:123)
          - [[TOPLEVEL.]MIDLEVEL.]LABEL (e.g. Receptor.Helix2.2x44)
          - COLOR (e.g. #FF0000 or white)

    Returns
    -------
    tuple[dict[str:str]]
        Both tuple-entries have residue identifiers as keys. The values are strings that hold the label and the
        CSS-compatible color string respectively.

    Raises
    ------
    AssertionError
        if a residue identifier (CHAIN:RESN:RESI) is specified twice in the file, or if a LABEL appears twice.
    """
    if label_file is None:
        return None

    ret = ({}, {})
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
    Convert atomic contacts into unique residue contacts. The interaction type is removed as well as any third or
    fourth atoms that are part of the interaction (e.g. water-bridges). Finally, the order of the residues within an
    interaction is such that the first is lexicographically smaller than the second ('A:ARG' comes before 'A:CYS').

    Example
    -------
        res_frequencies([
            [0, 'hbbb', 'A:ASN:108:O', 'A:ARG:110:N'],
            [0, 'vdw', 'A:ARG:110:N', 'A:ASN:108:CB']
        ])
        # => [[0, 'A:ARG:110', 'A:ASN:108']]

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
        if resi2 < resi1:
            resi1, resi2 = resi2, resi1
        frame_dict[frame].add((resi1, resi2))

    ret = []
    for frame in sorted(frame_dict):
        for resi1, resi2 in frame_dict[frame]:
            ret.append([frame, resi1, resi2])

    return ret


def parse_frequencyfiles(freq_files, freq_cutoff):
    import numpy as np

    columns = len(freq_files)
    ret = {}
    for fidx, freq_file in enumerate(freq_files):
        for line in freq_file:
            line = line.strip()
            if len(line) == 0 or line[0] == "#":
                continue

            tokens = line.split("\t")
            res1 = tokens[0]
            res2 = tokens[1]
            freq = float(tokens[2])

            if not (res1, res2) in ret:
                ret[(res1, res2)] = np.zeros(columns)

            ret[(res1, res2)][fidx] = freq

    # Remove entries where no frequency exceeds 0.6
    ret = {key: val for key, val in ret.items() if np.amax(val) > freq_cutoff}
    return ret


def gen_frequencies(count_list):
    """
    Take a list of residue contact counts (see output of `gen_counts`) and compute total counts and frequencies.

    Example:
        clist = [
            (4, {("A1", "R4"): 4, ("A1", "C5"): 3}),  # First simulation has 4 frames and two contacts
            (3, {("A1", "R4"): 2})                    # Second simulation has 3 frames and one contact
        ]
        gen_frequencies(clist)
        # Returns: (7, {("A1", "R4"): (6, 0.857), ("A1", "C5"): (3, 0.429)})

    Parameters
    ----------
    count_list: list of (int, dict of (str, str): int)
        List with individual frame counts and dictionaries mapping residue pairs to frame-counts

    Return
    ------
    (int, dict of (str, str): (int, float))
        Total framecount and mapping of residue ID pairs to the number of frames in which they contact and the frequency
    """
    from collections import defaultdict
    rescontact_count = defaultdict(int)
    total_frames = 0
    for frames, rescount_dict in count_list:
        total_frames += frames

        for (res1, res2), count in rescount_dict.items():
            rescontact_count[(res1, res2)] += count

    respair_freqs = {respair: (count, float(count) / total_frames) for respair, count in rescontact_count.items()}
    return total_frames, respair_freqs


def relabel(residue_contacts, residuelabels):
    """
    Relabel and filter residue contacts according to the residuelabel dictionary

    For example:
        contacts = [
            [0, 'A:ARG:110', 'A:ASN:108'],
            [0, 'A:FOL:200', 'A:ASN:108'],
            [0, 'A:NDP:201', 'A:ASN:108'],
            [1, 'A:ARG:110', 'A:NDP:201']
        ]
        labels = {'A:ARG:110': 'R110', 'A:ASN:108': 'N108', 'A:NDP:201': 'Ligand', 'A:FOL:200': 'Ligand'}
        relabel(contacts, labels)
        # => [(0, 'N108', 'R110'), (0, 'Ligand', 'N108'), (1, 'Ligand', 'R110')]

    Parameters
    ----------
    residue_contacts: list of list
        List of residue contacts specified as a framenumber and two residue identifiers
    residuelabels: dict of (str: str)
        Remaps and filters residuelabels, e.g. {"A:ARG:4": "R4"}

    Returns
    -------
    List of tuple
        Remapped and filtered residue-contacts specified as tuple of (int, string, string). Frames are guaranteed
        to be sorted
    """
    print('relabel')
    ret = set()
    for frame, res1, res2 in residue_contacts:
        if res1 in residuelabels and res2 in residuelabels:
            res1 = residuelabels[res1]
            res2 = residuelabels[res2]
            if res1 == res2:
                continue
            if res1 > res2:
                res1, res2 = res2, res1
            ret.add((frame, res1, res2))
    return sorted(list(ret), key=lambda c: c[0])


def gen_counts(residue_contacts):
    """
    Computer interaction-counts for each residue pair.

    For example:
        gen_counts([(0, 'Ligand', 'R110'), (0, 'Ligand', 'N108'), (1, 'Ligand', 'R110')])
        #  => { ("Ligand", "R110"): 2, ("Ligand", "N108"): 1 }

    Parameters
    ----------
    residue_contacts: list[tuple[int, str, str]]
        Residue interactions specified by tuple of frame and two residue ids

    Returns
    -------
    dict of (str, str): int
        Mapping of residue-residue interactions to frame-count
    """
    from collections import defaultdict
    rescontact_frames = defaultdict(set)
    for frame, res1, res2 in residue_contacts:
        rescontact_frames[(res1, res2)].add(frame)

    rescontact_counts = {(res1, res2): len(frames) for (res1, res2), frames in rescontact_frames.items()}
    return rescontact_counts


def gen_counts_old(input_lines, interaction_types=None, residuelabels=None):
    """
    Parse each line in `input_lines` as a line from get_*_contacts and return interaction-counts for each residue pair.
    If `residuelabels` is defined it is used to modify residue identifiers and to filter out residues not indicated.

    For example:
        inputs = [
            "# total_frames: 3",
            "\t".join(["0", "hbbb", "A:ALA:1:N", "A:ARG:4:O"]),
            "\t".join(["0", "vdw", "A:ALA:1:CB", "A:ARG:4:CA"]),
            "\t".join(["1", "vdw", "A:ALA:1:N", "A:CYS:5:CA"]),
            "\t".join(["2", "hbbb", "A:THR:2:N", "A:CYS:5:O"]),
            "\t".join(["2", "hbss", "A:ALA:1:N", "A:CYS:5:O"])
        ]
        labels = {"A:ALA:1": "A1", "A:ARG:4": "R4", "A:CYS:5": "C5"}

        # Only consider hbbb and vdw, filter away THR, and map to single-letter labels
        gen_counts(inputs, ["hbbb", "vdw"], labels)
        #  => { ("A1", "R4"): 1, ("A1", "C5"): 1 }

    Parameters
    ----------
    input_lines: Iterable[str]
        Interactions formatted as get_*_contacts output, e.g. ["0\thbbb\tA:ALA:1:N\tA:ARG:4:H", ...]
    interaction_types: list of str
        Which interaction types to consider
    residuelabels: dict of (str: str)
        Remaps and filters residuelabels, e.g. {"A:ARG:4": "R4"}

    Returns
    -------
    (int, dict of (str, str): int)
        Total frame-count and mapping of residue-residue interactions to frame-count
    """
    from collections import defaultdict

    def atomid_to_resid(atom):
        return atom[0:atom.rfind(":")]
        # return ':'.join(atom.split(':')[1:3])

    # Maps residue pairs to set of frames in which they're present
    rescontact_frames = defaultdict(set)
    total_frames = 0

    for line in input_lines:
        line = line.strip()
        if "total_frames" in line:
            tokens = line.split(" ")
            total_frames = int(tokens[1][tokens[1].find(":")+1:])

        if len(line) == 0 or line[0] == "#":
            continue

        tokens = line.split("\t")

        # Check that the interaction type is specified
        if interaction_types and tokens[1] not in interaction_types:
            continue

        frame = int(tokens[0])
        if frame + 1 > total_frames:
            total_frames = frame + 1

        res1 = atomid_to_resid(tokens[2])
        res2 = atomid_to_resid(tokens[3])

        # Change residue id according to `residuelabels` or skip if any of the residues are not present
        if residuelabels is not None:
            if res1 not in residuelabels or res2 not in residuelabels:
                continue
            res1 = residuelabels[res1]
            res2 = residuelabels[res2]

        # Ensure lexicographical order of residue names
        if res2 < res1:
            res1, res2 = res2, res1

        rescontact_frames[(res1, res2)].add(frame)

    # Insted of returning list of frames for each interaction, only return number of frames
    rescontact_counts = {(res1, res2): len(frames) for (res1, res2), frames in rescontact_frames.items()}
    return total_frames, rescontact_counts


__author__ = 'Rasmus Fonseca <fonseca.rasmus@gmail.com>'
__license__ = "APACHE2"
