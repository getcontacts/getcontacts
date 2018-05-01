__author__ = 'Rasmus Fonseca <fonseca.rasmus@gmail.com>'
__license__ = "APACHE2"


def parse_contacts(input_file):
    """
    Read a contact-file (tab-separated file with columns: frame, i-type, atomid1, atomid2[, atomid3[, atomid4]] and
    return it as a list of lists with frames converted to ints. The total number of frames is also returned

    Parameters
    ----------
    input_file: str
        Valid and readable input-file

    Returns
    -------
    (list of list, int)
        The list of interactions and the total number of frames
    """
    with open(input_file) as f:
        ret = []
        total_frames = 0
        for line in f:
            line = line.strip()
            if "total_frames" in line:
                tokens = line.split(" ")
                total_frames = int(tokens[1][tokens[1].find(":")+1:])

            if len(line) == 0 or line[0] == "#":
                continue

            tokens = line.split("\t")
            tokens[0] = int(tokens[0])

            ret.append(tokens)

        return ret, total_frames


def res_contacts(contacts):
    """
    Given a list of atomic contacts, compute the residue contacts

    Example
    -------
        res_frequencies([[0, 'hbbb', 'A:ASN:108:O', 'A:ARG:110:N'], [0, 'vdw', 'A:ARG:110:N', 'A:ASN:108:CB']])
        # ==> [[0, 'A:ASN:108', 'A:ARG:110']]

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

    # todo convert to list


