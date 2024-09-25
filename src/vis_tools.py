import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgba
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import numpy as np
import networkx as nx

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
