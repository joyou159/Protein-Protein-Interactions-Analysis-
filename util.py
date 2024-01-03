from matplotlib.lines import Line2D
from matplotlib.colors import to_rgba
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import copy
from math import log10
from collections import Counter
import random
from bioservices import UniProt
import warnings
warnings.filterwarnings("ignore")


def visualize_graph(df):
    PPI_graph = nx.DiGraph()
    for i in range(len(df)):
        tail = df.loc[i, :][0]
        head = df.loc[i, :][1]
        weight = df.loc[i, :][2]
        PPI_graph.add_edge(tail, head, weight=weight)
    betCent = nx.betweenness_centrality(
        PPI_graph, normalized=True, endpoints=True)
    node_color = [20000.0 * PPI_graph.degree(v) for v in PPI_graph]
    node_size = [v * 1000000 for v in betCent.values()]
    plt.figure(figsize=(40, 40))
    pos = nx.spring_layout(PPI_graph)
    nx.draw(PPI_graph, node_color=node_color, edge_color="gray", pos=pos,
            font_color='white', font_weight='bold', node_size=node_size)
    plt.show()
    return PPI_graph


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


def plot_graph(G, pos, ax, path_nodes, total_scores):
    """
    Plot the directed graph along with path information.
    Parameters:
    - G (DiGraph): A directed graph
    - pos (dict): A dictionary with node positions for the graph
    - ax (matplotlib.axes.Axes): The axes on which to draw the graph
    - path_nodes (list): A list of paths, where each path is represented as a list of nodes
    - total_scores (list): A list of total scores corresponding to each path
    """
    colormap = plt.cm.get_cmap('viridis', len(path_nodes))
    table_data = [["Total Score"]]
    custom_legend = []
    for i, path in enumerate(path_nodes):
        path_edges = list(zip(path[:-1], path[1:]))
        color = to_rgba(colormap(i))
        nx.draw_networkx_nodes(G, pos, node_size=200,
                               node_color='lightblue', ax=ax)
        nx.draw_networkx_edges(G, pos, edgelist=path_edges,
                               edge_color=[color], ax=ax)
        nx.draw_networkx_edge_labels(
            G, pos, edge_labels=nx.get_edge_attributes(G, 'weight'))
        nx.draw_networkx_labels(G, pos, labels={
                                node: node for node in G.nodes()}, font_size=8, font_color='black', ax=ax)
        table_data.append([f"{total_scores[i]:.4f}"])
        custom_legend.append(
            Line2D([0], [0], color=color, lw=2, label=f"Path {i + 1}"))
    ax.legend(handles=custom_legend, loc='upper left')
    table = ax.table(cellText=table_data, colLabels=None,
                     cellLoc='center', loc='bottom', bbox=[0, -0.1, 1, 0.1])
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    ax.set_title("DAC with K-Shortest-Paths")


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


def Export_PDF(graph, path_nodes, total_scores, pdf_filename):
    """
    Export a PDF file containing a plot of the directed graph with path information.
    Parameters:
    - graph (DiGraph): A directed graph representing the network.
    - path_nodes (list): A list of paths, where each path is represented as a list of nodes
    - total_scores (list): A list of total scores corresponding to each path
    - pdf_filename (str): The name of the PDF file to be created
    """

    if not np.any(path_nodes):
        print("terminate")
        return
    edge_weights = nx.get_edge_attributes(graph, 'weight')
    G = generate_graph(path_nodes, edge_weights)
    pos = nx.fruchterman_reingold_layout(G)
    with PdfPages(pdf_filename) as pdf:
        fig, ax = plt.subplots(figsize=(5, 5))
        plot_graph(G, pos, ax, path_nodes, total_scores)
        pdf.savefig(fig)
        plt.close(fig)


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


def find_acyclic_K_shortest_paths(graph, source, target, pdf_filename, k=100):
    """
    Find the top K acyclic shortest paths in a directed graph and export a PDF with visualization.
    Parameters:
    - graph (DiGraph): A directed graph.
    - source (node): The starting node for paths.
    - target (node): The target node for paths.
    - pdf_filename (str): The name of the PDF file to be created.
    - k (int): The number of top paths to find (default is 100).
    """
    graph_copy = graph.copy()
    for u, v, data in graph_copy.edges(data=True):
        probability = data.get("weight")
        if probability == 0:
            data['weight'] = np.Inf
        else:
            data['weight'] = round(-log10(probability), 3)

    K_shortest_paths = list(nx.shortest_simple_paths(
        graph_copy, source, target, weight="weight"))  # using Yen's algorithm
    print(f"All paths before cyclic filtration:{K_shortest_paths}")
    K_acyclic_shortest_paths = np.array(
        [path for path in K_shortest_paths if is_acyclic(graph_copy, path)], dtype=object)[:k]
    print(f"All paths after cyclic filtration:{K_acyclic_shortest_paths}")
    paths_score = np.array([calculate_path_score(graph_copy, path)[
        0] for path in K_acyclic_shortest_paths])
    non_zero_paths_score = np.where(paths_score != 0)[0]
    K_acyclic_shortest_paths = K_acyclic_shortest_paths[non_zero_paths_score]
    paths_score = paths_score[non_zero_paths_score]
    Export_PDF(graph_copy, K_acyclic_shortest_paths, paths_score, pdf_filename)


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


def export_protein_data(graph, protein_of_interest, file_name):
    """
    Export protein data, including in-degree and out-degree information, to a text file.
    Parameters:
        graph (NetworkX graph): The graph containing protein interaction data.
        protein_of_interest (str): The name of the target protein.
        file_name (str): The name of the output text file.
    Returns:
        None: The function writes the protein data to the specified file.
    Notes:
        The function calculates the in-degree and out-degree of the target protein,
        retrieves the connections with their weights, and writes the information
        to the specified text file. If there are no in-degree or out-degree connections,
        appropriate messages are written to the file.
    """
    in_degree_data = get_node_predecessors(graph, protein_of_interest)
    out_degree_data = get_node_successors(graph, protein_of_interest)
    with open(file_name, 'w') as file:
        file.write(f"target protein name:{protein_of_interest}\n")
        file.write(f"# in degree: {in_degree_data[0]}\n")
        file.write(f"# out degree: {out_degree_data[0]}\n")
        file.write("\nin_degree_connections::\n")
        if in_degree_data[0] > 0:  # Check if in_degree is greater than 0
            for connection in in_degree_data[1]:
                file.write(f"{connection[0]} || {connection[1]}\n")
        else:
            file.write("There's No In Degree Connection\n")
        file.write("\nout_degree_connections::\n")
        if out_degree_data[0] > 0:  # Check if out_degree is greater than 0
            for connection in out_degree_data[1]:
                file.write(f"{connection[0]} || {connection[1]}\n")
        else:
            file.write("There's No Out Degree Connection\n")


def draw_degree_histogram(graph, protein_set, output_pdf, flag="all"):
    """
    Draw a histogram of protein degrees and save the sorted degree data to a text file.
    Parameters:
        graph (NetworkX graph): The graph containing protein interaction data.
        protein_set (list): List of protein names.
        output_pdf (str): The name of the output PDF file for the histogram and degree data.
        flag (str): Specifies the type of degree to consider ("all", "in", or "out").
    Returns:
        None: The function displays the histogram and saves the degree data to a text file.
    Notes:
        The function calculates the specified type of degree for each protein in the input set,
        draws a histogram, and saves the sorted degree data to a text file. The histogram is displayed,
        and the degree data includes the protein names and their corresponding degrees.
    """
    degree_function = {"all": graph.degree,
                       "in": graph.in_degree,
                       "out": graph.out_degree}
    protein_degrees = np.array(
        [degree_function[flag](protein) for protein in protein_set])
    sorted_degrees_indices = np.argsort(protein_degrees)[::-1]
    sorted_degrees = protein_degrees[sorted_degrees_indices]
    protein_set = np.array(protein_set)[sorted_degrees_indices]
    _, ax = plt.subplots(figsize=(5, 5))
    ax.hist(sorted_degrees, bins=30, color='skyblue',
            edgecolor='black',  density=True)
    ax.set_title('Degree Of Protein Histogram')
    ax.set_xlabel('Total Degree')
    ax.set_ylabel('Frequency')
    plt.show()
    with open(output_pdf, 'w') as degree_file:
        degree_file.write(f"protein_name  {flag}_degree\n")
        degree_file.write("-----------------------\n")
        for protein, degree in zip(protein_set, sorted_degrees):
            degree_file.write(f"{protein} \t\t  {degree}\n")


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


def saving_as_unweighed_graph_as_adjMatrix(graph, file_name):
    unweighted_graph = graph.copy()
    for _, _, data in unweighted_graph.edges(data=True):
        data["weight"] = 1
    adjacency_matrix = nx.adjacency_matrix(unweighted_graph)
    df = pd.DataFrame.sparse.from_spmatrix(
        adjacency_matrix, index=list(graph.nodes()), columns=list(graph.nodes()))
    df.to_csv(file_name)


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
