#!/bin/sh
# Shebang-hack for launching pymol
''':'
exec pymol -q "$0" -- "$@"
'''

"""
PyMol visualization of weighted contact network and communication pathways. 

Inputs:
- topology file
- contact frequency file 
- edge list for communication pathway or visualizing portions of full network 

"""

def fix_amino_acid_names(key):
    key = key.replace("HSD", "HIS")
    key = key.replace("HSE", "HIS")
    key = key.replace("HSP", "HIS")
    key = key.replace("HIE", "HIS")
    key = key.replace("HIP", "HIS")
    key = key.replace("HID", "HIS")
    key = key.replace("GLH", "GLU")
    key = key.replace("ASH", "ASP")
    key = key.replace("CYP", "CYS")
    key = key.replace("CYX", "CYS")
    return key


def get_edge_to_weight_map(edge_weights):
    """
    Map the edge between a pair of residues to weight

    Parameters
    ----------
    edge_weights: string
        Path to table associated edge between a pair of amino acids to an edge weight

    Returns
    -------
    edge_weight_map: dictionary
        {("B:SER:60", "B:TYR:99"): 0.64032, ...}
    """

    edge_weight_map = {}
    f = open(edge_weights, 'r')

    for line in f:
        linfo = line.strip().split("\t")
        aa1, aa2, freq = linfo[0], linfo[1], float(linfo[2])
        aa1 = fix_amino_acid_names(aa1)
        aa2 = fix_amino_acid_names(aa2)
        edge_weight_map[(aa1, aa2)] = freq

    return edge_weight_map


def gen_atom_selection_str(top_name, node):
    """
    Generate atom selection string

    Parameters
    ----------
    top_name: string
        Name of the protein structure
    node: string
        ":" delimited atom or residue string

    Return
    ------
    selection_str: string
        Selection string
    """

    node_info = node.split(":")

    # Residue level interactions
    if len(node_info) == 3:
        chain, resn, resi = node_info
        selection_str = "%s and chain %s and resn %s and resi %s and name CA" % (top_name, chain, resn, resi)

    # Atomic level interactions
    elif len(node_info) == 4:
        chain, resn, resi, name = node_info
        selection_str = "%s and chain %s and resn %s and resi %s and name %s" % (top_name, chain, resn, resi, name)

    return selection_str


def get_line_thickness(weight):
    """
    Determine radius of cylinder based on weight of contact edge

    Parameters
    ----------
    weight: float
        Edge weight

    Return
    ------
    radius: float
        Radius of sphere
    """
    min_width = 0.01
    max_width = 0.5
    radius = float(weight)*max_width
    if radius < min_width:
        radius = min_width
    return radius


def draw_edge(top_name, node1, node2, weight):
    """
    Draw an edge between two nodes with specified weight

    Parameters
    ----------
    top_name: string
        Name of the protein structure
    node1: string
        Represents first residue or atom
    node2: string
        Represents second residue or atom
    weight: float
        Edge between between node1 and node2
    """

    label = "%s_contact_network" % top_name
    sel1 = gen_atom_selection_str(top_name, node1)
    sel2 = gen_atom_selection_str(top_name, node2)

    # Draw cylinder
    xyz1 = cmd.get_coords(sel1, 1)[0]
    xyz2 = cmd.get_coords(sel2, 1)[0]
    x1, y1, z1 = xyz1
    x2, y2, z2 = xyz2
    r1, g1, b1 = 0.5, 0.5, 0.5  # color (red)
    r2, g2, b2 = 0.5, 0.5, 0.5  # color (red)
    radius = get_line_thickness(weight)
    cmd.load_cgo([9.0, x1, y1, z1, x2, y2, z2, radius, r1, g1, b1, r2, g2, b2], "cylinder_%s_%s" % (node1, node2))

    cmd.do("show spheres, %s" % sel1)
    cmd.do("show spheres, %s" % sel2)


def visualize_protein_network(structure, edge_weights, sub_network=[], cutoff=0.0):
    """
    PyMol rendering of protein structure network with edge weights derived from MD simulation

    Parameters
    ----------
    structure: string
        Path to structure to superimpose network on
    edge_weights: string
        Path to table associated edge between a pair of amino acids to an edge weight
    sub_network: list of tuples
        Default is empty array. User can use this argument to draw a subnetwork or shortest path
        between a pair of nodes
    cutoff: float
        Default 0.0, cutoff value for frequency of contact
    """

    # Load in background protein structure
    top_name = top.split("/")[-1].split(".")[0]
    cmd.bg_color("white")
    cmd.load(top)
    cmd.hide()
    cmd.show("cartoon")
    cmd.cartoon("loop")
    cmd.color("white", top_name)
    cmd.do("set cartoon_transparency, 0.7")
    cmd.do("set sphere_scale, 0.35")
    cmd.do("set sphere_transparency, 0.5")
    cmd.do("set sphere_color, deepteal")

    # Draw protein network with edge thickness representing the frequency of a contact in simulation
    edge_weight_map = get_edge_to_weight_map(edge_weights)

    if not sub_network:  # If sub_network argument empty then draw entire network
        for edge in edge_weight_map:
            aa1, aa2 = edge
            weight = edge_weight_map[edge]
            if aa1 == aa2:
                continue
            if weight < cutoff:
                continue
            draw_edge(top_name, aa1, aa2, weight)
    else:
        for aa1, aa2 in sub_network:
            aa1 = fix_amino_acid_names(aa1)
            aa2 = fix_amino_acid_names(aa2)
            k1, k2 = (aa1, aa2), (aa2, aa1)
            if k1 in edge_weight_map:
                weight = edge_weight_map[k1]
            elif k2 in edge_weight_map:
                weight = edge_weight_map[k2]
            if aa1 == aa2:
                continue
            if weight < cutoff:
                continue
            draw_edge(top_name, aa1, aa2, weight)


# Example 
top = "../example/5xnd_topology.pdb"
contact_freq = "../example/5xnd_contact_freq.tsv"
sub_network = [('A:VAL:44', 'A:SER:40'), ('A:PHE:71', 'A:VAL:12'), ('A:ASP:26', 'A:SER:24'), ('A:VAL:12', 'A:ALA:10'), ('A:PHE:67', 'A:GLN:69'), ('A:LEU:64', 'A:GLU:61'), ('A:PHE:31', 'A:LYS:28'), ('A:ASN:70', 'A:MET:1'), ('A:ASP:54', 'A:LYS:55'), ('A:PHE:71', 'A:LEU:16'), ('A:PHE:67', 'A:GLU:63'), ('A:LYS:72', 'A:GLY:74'), ('A:PHE:86', 'A:LYS:88'), ('A:ILE:51', 'A:GLN:53'), ('A:LEU:64', 'A:GLU:62'), ('A:VAL:44', 'A:ASP:42'), ('A:GLY:23', 'A:ALA:21'), ('A:ILE:50', 'A:VAL:6'), ('A:ILE:107', 'A:ALA:104'), ('A:PHE:67', 'A:PHE:31'), ('A:PHE:86', 'A:LYS:84'), ('A:PHE:31', 'A:SER:37'), ('A:PHE:71', 'A:LEU:36'), ('A:PHE:67', 'A:LEU:78'), ('A:ASN:70', 'A:SER:0'), ('A:PHE:31', 'A:ALA:33'), ('A:LEU:87', 'A:GLY:96'), ('A:PHE:103', 'A:ILE:100'), ('A:LYS:28', 'A:GLY:109'), ('A:PHE:67', 'A:ILE:51'), ('A:PHE:86', 'A:ALA:89'), ('A:GLU:63', 'A:ASP:52'), ('A:ILE:59', 'A:LYS:97'), ('A:GLY:96', 'A:ASP:95'), ('A:ILE:51', 'A:PHE:103'), ('A:LEU:36', 'A:GLU:43'), ('A:ASN:70', 'A:ALA:2'), ('A:LEU:16', 'A:LYS:20'), ('A:ASN:70', 'A:PHE:3'), ('A:GLU:82', 'A:ALA:22'), ('A:PHE:25', 'A:GLY:23'), ('A:ALA:15', 'A:ASP:17'), ('A:PHE:86', 'A:THR:83'), ('A:PHE:86', 'A:GLY:90'), ('A:ILE:51', 'A:ALA:47'), ('A:LEU:78', 'A:ASP:80'), ('A:ASN:70', 'A:ILE:50'), ('A:GLN:69', 'A:ALA:75'), ('A:LEU:87', 'A:ASP:91'), ('A:PHE:67', 'A:LEU:64'), ('A:PHE:86', 'A:ILE:98'), ('A:VAL:12', 'A:LYS:8'), ('A:ALA:15', 'A:THR:13'), ('A:ILE:51', 'A:ALA:49'), ('A:LEU:16', 'A:ALA:14'), ('A:ASN:70', 'A:LEU:66'), ('A:PHE:67', 'A:ASN:70'), ('A:ASN:70', 'A:LYS:72'), ('A:ILE:51', 'A:ILE:59'), ('A:LEU:87', 'A:ALA:85'), ('A:GLU:102', 'A:ASP:93'), ('A:PHE:30', 'A:ALA:15'), ('A:LEU:78', 'A:CYS:19'), ('A:ILE:50', 'A:LYS:46'), ('A:PHE:31', 'A:ILE:107'), ('A:GLY:96', 'A:GLY:94'), ('A:LEU:64', 'A:LEU:87'), ('A:LYS:65', 'A:SER:79'), ('A:MET:106', 'A:LYS:108'), ('A:PHE:86', 'A:HIS:27'), ('A:PHE:71', 'A:CYS:34'), ('A:GLY:35', 'A:LYS:39'), ('A:LEU:78', 'A:ARG:76'), ('A:PHE:30', 'A:LYS:32'), ('A:VAL:12', 'A:ASP:9'), ('A:PHE:48', 'A:LYS:45'), ('A:GLU:63', 'A:PHE:58'), ('A:PHE:103', 'A:ASP:101'), ('A:PHE:48', 'A:GLY:57'), ('A:LEU:64', 'A:GLU:60'), ('A:GLY:90', 'A:SER:92'), ('A:PHE:31', 'A:LYS:29'), ('A:PHE:67', 'A:LYS:65'), ('A:PHE:30', 'A:ASP:26'), ('A:ALA:2', 'A:SER:5'), ('A:ILE:51', 'A:PHE:48'), ('A:GLU:63', 'A:ASP:54'), ('A:ASN:70', 'A:LEU:68'), ('A:ALA:75', 'A:ALA:73'), ('A:ALA:2', 'A:ALA:4'), ('A:PHE:3', 'A:LEU:7'), ('A:ALA:33', 'A:GLU:11'), ('A:PHE:31', 'A:GLY:35'), ('A:PHE:86', 'A:GLU:82'), ('A:GLY:90', 'A:GLU:102'), ('A:VAL:44', 'A:THR:41'), ('A:ALA:15', 'A:GLY:18'), ('A:ILE:107', 'A:ALA:105'), ('A:PHE:67', 'A:PHE:30'), ('A:LYS:84', 'A:ALA:81'), ('A:PHE:67', 'A:PHE:86'), ('A:PHE:67', 'A:PHE:71'), ('A:PHE:58', 'A:GLY:99'), ('A:GLU:63', 'A:SER:56'), ('A:PHE:86', 'A:MET:106'), ('A:PHE:67', 'A:PHE:25'), ('A:LYS:65', 'A:ALA:77'), ('A:ILE:107', 'A:VAL:44'), ('A:GLY:35', 'A:GLY:38')]
visualize_protein_network(top, contact_freq, sub_network)
