from vis_tools import *
import numpy as np
import copy
import random
from bioservices import UniProt

def generate_graph(path_nodes, edge_weights):
    """
    Generate a directed graph based on a list of paths
    Parameters:
    - path_nodes (list): A list of paths, where each path is represented as a list of nodes
    - edge_weights (list): A list of lists containing edge weights corresponding to each path
    Returns:
    - G (DiGraph): A directed graph representing the paths
    """
    G = nx.DiGraph()
    for path in path_nodes:
        for (source, target) in zip(path[:-1], path[1:]):
            weight = edge_weights[(source, target)]
            G.add_edge(source, target, weight=weight)
    return G


def calculate_path_score(graph, path):
    """
    Calculate the score of a path in a graph based on edge weights.
    Parameters:
        graph (NetworkX graph): The graph containing edge weights.
        path (list): The list of nodes representing the path.
    Returns:
        tuple: A tuple containing the total score of the path and a list of individual edge weights.
    """
    score_list = list()
    score = 1
    for i in range(len(path) - 1):
        u, v = path[i], path[i + 1]
        if graph.has_edge(u, v):
            weight = graph.get_edge_data(u, v)['weight']
            score *= 10**(-weight)
            score_list.append(round(10**(-weight), 3))
    return (round(score, 3), score_list)


def is_acyclic(graph, path):
    # Check if adding the next edge creates a cycle
    graph_copy = copy.deepcopy(graph)
    for i in range(len(path) - 1):  # iterating on a window of 2 nodes
        u, v = path[i], path[i + 1]
        if graph_copy.has_edge(u, v):
            graph_copy.remove_edge(u, v)
            if nx.has_path(graph, v, u):  # Check if there is a path back to the source
                graph_copy.add_edge(u, v)
                print(
                    f"there is a feedback loop between: {(u , v)}")
                return False
            graph_copy.add_edge(u, v)
    return True


def generate_path(graph):
    no_nodes = graph.number_of_nodes()
    list_of_nodes = list(graph.nodes())
    while True:  # iterating indefinitely util finding a testing path
        # select two nodes randomly and uniformly
        candidate_nodes_indices = random.sample(range(no_nodes), 2)
        v, u = list_of_nodes[candidate_nodes_indices[0]
                             ], list_of_nodes[candidate_nodes_indices[1]]
        # if the two nodes has no edge and there is path between them, to get a large path.
        if not graph.has_edge(v, u) and nx.has_path(graph, v, u):
            return v, u


def get_node_predecessors(graph, protein_of_interest):
    k = graph.in_degree(protein_of_interest)  # degree of such node
    predecessors = list(graph.predecessors(protein_of_interest))
    weight_list = list()
    edge_weights = nx.get_edge_attributes(graph, "weight")
    for protein in predecessors:
        weight_list.append(edge_weights[(protein, protein_of_interest)])
    # in_degree, in_conneced_proteins
    return k, list(zip(predecessors, weight_list))


def get_node_successors(graph, protein_of_interest):
    k = graph.out_degree(protein_of_interest)  # degree of such node
    successors = list(graph.successors(protein_of_interest))
    weight_list = list()
    edge_weights = nx.get_edge_attributes(graph, "weight")
    for protein in successors:
        weight_list.append(edge_weights[(protein_of_interest, protein)])
    # out_degree, out_conneced_proteins
    return k, list(zip(successors, weight_list))


def get_conversion_map(uniprot_ids):
    """
    Retrieve gene names associated with UniProt IDs using the UniProt service.
    Parameters:
    - uniprot_ids (list): A list of UniProt IDs for which gene names are to be fetched.
    Returns:
    - dict: A dictionary mapping each UniProt ID to its associated gene names.
    """
    lookup = dict()
    u = UniProt()
    for protein in uniprot_ids:
        gene_namse = u.get_df(protein)["Gene Names"]
        if isinstance(gene_namse.values[0], str):
            lookup[protein] = gene_namse.str.split(" ").values[0]
        else:
            lookup[protein] = gene_namse.values[0]
    return lookup


def basic_net_stats(G):
    print(f"number of nodes {G.number_of_nodes()}")
    print(f"number of edges {G.number_of_edges()}")
    degree_sequence = [d for n, d in G.degree()]
    print(f"average degree {(np.mean(degree_sequence)):.2f}")


def get_graph_hubs(G, k=10):
    degree_sequence = dict(G.degree())
    sorted_degrees = dict(sorted(degree_sequence.items(),
                          key=lambda item: item[1], reverse=True)[:k])
    return sorted_degrees
