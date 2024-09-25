import numpy as np
from vis_tools import * 
from math import log10
from util import (generate_graph,
                    calculate_path_score,
                    is_acyclic, get_node_predecessors,
                    get_node_successors) 

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


def saving_as_unweighed_graph_as_adjMatrix(graph, file_name):
    unweighted_graph = graph.copy()
    for _, _, data in unweighted_graph.edges(data=True):
        data["weight"] = 1
    adjacency_matrix = nx.adjacency_matrix(unweighted_graph)
    df = pd.DataFrame.sparse.from_spmatrix(
        adjacency_matrix, index=list(graph.nodes()), columns=list(graph.nodes()))
    df.to_csv(file_name)
