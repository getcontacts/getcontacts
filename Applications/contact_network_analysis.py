#!/usr/bin/env python3

"""
Suite of network analyses tools for identifying high centrality protein residues and 
communication pathways. 

Example:

from contact_network_analysis import *
contact_freq="../example/5xnd_contact_freq.tsv"

G = create_graph(contact_freq)
betweenness_centrality_dist(G)
degree_centrality_dist(G, True)

communication_pathway(G, ["A:PHE:67"])
sp_edges = communication_pathway(G, ["A:PHE:67"], ["A:CYS:19"])

See contact_network_pymol_viz.py for visualization
"""

import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns 


# Utils
def get_edge_weight(graph, node1, node2):
    """
    Returns edge weight between two nodes

    Parameters
    ----------
    graph: networkx object
        protein interaction graph
    node1: string
        residue 1
    node2: string
        residue 2

    Return
    ------
    output: float
        interaction frequency between two residues
    """
    try:
        freq = graph[node1][node2]
        return freq
    except KeyError:
        print("Key Error: Edge between %s -- %s doesn't exist" % (node1, node2))


def create_graph(contact_frequency):
    """
    Create networkx graph from edge list with contact frequencies

    Parameters
    ----------
    contact_frequency: string
        Path to contact_frequency file containing weighted edge list in format of <res1> <res2> <frequency>

    Return
    ------
    graph: networkx object
        Residue interaction graph object

    """

    # Parse the contact_frequency graph file
    f = open(contact_frequency, 'r')
    nodes, edges = set(), set()
    for line in f:
        linfo = line.strip().split("\t")
        res1 = linfo[0]
        res2 = linfo[1]
        freq = float(linfo[2])

        nodes.add(res1)
        nodes.add(res2)
        edges.add((res1, res2, freq))

    # Construct networkx graph
    graph = nx.Graph()
    for res in nodes:
        graph.add_node(res)

    for res1, res2, freq in edges:
        graph.add_edge(res1, res2, weight=freq)

    return graph


# Network Analysis
def betweenness_centrality_dist(graph, plot=False):
    """
    Calculate and plot distribution for betweenness centrality for all nodes in protein network

    Parameters
    ----------
    graph: networkx object
        Residue interaction graph object
    plot: bool
        Default = False. Display distribution
    Returns
    -------
    bc: dictionary
        Mapping between node to betweenness centrality values
    """

    bc = nx.betweenness_centrality(graph)
    centrality_values = []
    for key, value in reversed(sorted(bc.items(), key=lambda item: (item[1], item[0]))):
        print("%s: %s" % (key, value))
        centrality_values += [value]
    if plot:
        sns.distplot(centrality_values)
        plt.xlabel("Betweenness Centrality")
        plt.ylabel("Frequency")
        plt.title("Betweenness Centrality Distribution")
        plt.show()
    return bc


def degree_centrality_dist(graph, plot=False):
    """
    Calculate and plot distribution for degree centrality for all nodes in protein network

    Parameters
    ----------
    graph: networkx object
        Residue interaction graph object
    plot: bool
        Default = False. Display distribution
    Returns
    -------
    dc: dictionary
        Mapping between node to degree centrality values
    """

    dc = nx.degree_centrality(graph)
    centrality_values = []
    for key, value in reversed(sorted(dc.items(), key=lambda item: (item[1], item[0]))):
        print("%s: %s" % (key, value))
        centrality_values += [value]
    if plot:
        sns.distplot(centrality_values)
        plt.xlabel("Degree Centrality")
        plt.ylabel("Frequency")
        plt.title("Degree Centrality Distribution")
        plt.show()
    return dc


def communication_pathway(graph, src_nodes, dest_nodes=[], draw=False):
    """
    Calculate shortest pathways between specified source and destination nodes in the weighted protein network.

    Parameters
    ----------
    graph: networkx object
        Residue interaction graph object
    src_nodes: list of strings
        Source nodes
    dest_nodes: list of strings
        Default = []. Destination nodes
    draw: bool
        Default = False. Interface with PyMol to launch session/generate .pse to superimpose graph on the protein.
        (Have yet to connect module to contact_network_pymol_viz.py)
    Returns
    -------
    sp_edges: list of tuples
        Generate list of edges representing the subnetwork corresponding to shortest communication pathways.
    """

    sp_edges = set()
    if not dest_nodes:  # Edge case: if no destination nodes specified, then show all shortest pathways
        for src in src_nodes:
            pathways = nx.shortest_path(graph, src)
            for dest in pathways:
                sp = pathways[dest]
                for idx in range(len(sp) - 1):
                    sp_edges.add((sp[idx], sp[idx + 1]))
    else:
        for src in src_nodes:
            for dest in dest_nodes:
                sp = nx.shortest_path(graph, src, dest)
                for idx in range(len(sp) - 1):
                    sp_edges.add((sp[idx], sp[idx + 1]))
    sp_edges = list(sp_edges)
    return sp_edges
