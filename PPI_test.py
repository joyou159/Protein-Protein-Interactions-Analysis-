from src.vis_tools import * 
from src.report import *
from src.util import (basic_net_stats, 
                  get_graph_hubs,
                  generate_path, 
                  get_conversion_map) 
import copy
from collections import Counter
from bioservices import UniProt



file_path = "Dataset\PathLinker_2018_human-ppi-weighted-cap0_75.txt"
df = pd.read_csv(file_path, sep='\t')
df = df[:5000]  # taking a portion of the graph

print("------------------------------Data Inspection and Visualization---------------------------------------")
# checking for any NaN values in each row and column
print("Check for any missing data: ", df.isna().any().any())  # return false

# checking the range of probabilities in edge_weight column
print(f"The minimum interaction probability {min(df['edge_weight'])}")  # 0
print(f"The maximum interaction probability {max(df['edge_weight'])}")  # 0.75

# first required point
PPI_graph = visualize_graph(df)  # custom function for visualization

basic_net_stats(PPI_graph)  # For the whole graph
""" number of nodes 17168
number of edges 612516
 average degree 71.36 """

# return the most important 20 hubs in the PPIN
hubs = get_graph_hubs(PPI_graph, 20)
print(f"The highest proteins in degree (hubs) {hubs}")

# Strongly connected component in the PPIN for Diameter computation.
components_PPI = np.array([PPI_graph.subgraph(
    c) for c in nx.strongly_connected_components(PPI_graph)], dtype=object)
sorted_indices = np.argsort([len(comp) for comp in components_PPI])[::-1]
max_comp_graph = components_PPI[sorted_indices]
print(f"The 1st largest connected Component: {max_comp_graph[0]}")
print(f"The 2nd largest connected Component: {max_comp_graph[1]}")


# for the whole graph it's equal to 8
graph_diameter = nx.diameter(max_comp_graph[0])
print(f"The maximum shortest path in the PPIN of {graph_diameter} Steps")


# transitivity of the interactome as directed as it is.
# transitivity_value = nx.transitivity(PPI_graph) # computational exhaustive, equal to 0.06446 for whole graph

# Transitivity of the interactome as undirected compared to random graph (Transitivity is a measure of modularity as a global measure)
PPI_graph_undirected = copy.deepcopy(PPI_graph).to_undirected()
original_transitivity = nx.transitivity(PPI_graph_undirected)
print("Original PPIN transitivity Coefficient (no direction considered):",
      original_transitivity)  # 0.0679 for whole graph

random_graph = nx.gnm_random_graph(
    PPI_graph.number_of_nodes(), PPI_graph.number_of_edges())
random_transitivity = nx.transitivity(random_graph)
print("Random Graph transitivity Coefficient (no direction considered):",
      random_transitivity)  # 0.00414 for whole graph


# Clustering coefficient (local measure)
original_clustering_coefficient = nx.average_clustering(PPI_graph_undirected)
random_clustering_coefficient = nx.average_clustering(random_graph)
print("Original PPIN Clustering Coefficient (no direction considered):",
      original_clustering_coefficient)  # 0.135 # for whole graph
print("Random Graph Clustering Coefficient (no direction considered):",
      random_clustering_coefficient)  # 0.0021 for whole graph

print("------------------------------k-Cyclic Shortest Paths---------------------------------------")

# testing example for acyclic shortest path computation (taken from lecture's slides)
G = nx.DiGraph()
G.add_edges_from([(1, 2, {'weight': 0.1}), (2, 3, {'weight': 0.2}), (3, 4, {
                 'weight': 0.3}), (1, 5, {'weight': 0.3}), (5, 4, {'weight': 0.5})])
# ,(2, 3, {'weight': 1})
# Specify the layout for better visualization
pos = nx.circular_layout(G)

# Draw the graph
# nx.draw(G, pos, with_labels=True, node_size=700, node_color='skyblue', font_size=10, font_color='black', font_weight='bold')

# Draw edge labels
edge_labels = nx.get_edge_attributes(G, 'weight')
# nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels);

file_path = "Exported_files/acyclic_shortest_paths_test1.pdf"
find_acyclic_K_shortest_paths(G, 1, 4, file_path)


# generate a path and get the shortest path
source, target = generate_path(PPI_graph)
print("\nThe source and target proteins generated randomly in order", source, target)
file_path = "Exported_files/acyclic_shortest_paths_test2.pdf"
find_acyclic_K_shortest_paths(PPI_graph, source, target, file_path, 2)


# these source and target proteins have been deduced from generate_path function.
source = "Q00535"
target = "P14373"
file_path = "Exported_files/acyclic_shortest_paths_test3.pdf"
find_acyclic_K_shortest_paths(PPI_graph, source, target, file_path, 2)
# the actual path: ['Q00535', 'Q9NXR1', 'P14373']

# no feedback loop between the first two nodes
print("checking if there is a feedback loop between two nodes in a path")
print(
    f'check if there is an edge between ("Q00535", "Q9NXR1") {PPI_graph.has_edge("Q00535", "Q9NXR1")}')
print(
    f'check if there a path back from ("Q9NXR1", "Q00535") {nx.has_path(PPI_graph, "Q9NXR1", "Q00535")}')

# list paths' (source & target) that are acyclic  all of them of length 3 (small world effect)
path_list = [
    ("Q9NXR8", "Q92793"), ("Q00535", "P14373"), ("Q00535",
                                                 "Q9BRK4"), ("Q92831", "A8MW92"), ("Q92831", "Q9UBN7"),
    ("Q9NXR8", "Q9Y618"), ("Q9NXR8", "Q9UKV0"), ("Q9NXR8",
                                                 "O14744"), ("Q92831", "Q9BY41"),
    ("Q92831", "P28749"), ("Q9NXR8", "P02751"), ("Q92831",
                                                 "Q15573"), ("Q9NXR8", "O75376"),
    ("Q9NXR8", "Q96NS5"), ("Q9NXR8", "P00533"), ("Q00535",
                                                 "Q08379"), ("Q92831", "Q93034"),
    ("Q9NXR8", "P10070"), ("Q9NXR8", "Q9BW71"), ("Q9NXR8",
                                                 "Q04206"), ("Q92831", "Q8IZL8"),
    ("Q9NXR8", "Q8WTS6"), ("Q92831", "P11509"), ("Q9NXR8",
                                                 "P33993"), ("Q9NXR8", "Q14997"),
    ("Q9NXR8", "Q9UBU9"), ("Q9NXR8", "Q14566"), ("Q92831",
                                                 "Q6PL18"), ("Q9NXR8", "O14929"),
    ("Q92831", "Q14134"), ("Q92831", "A8MW92")]

# candidates for hubs ("Q9NXR8", "Q00535", "Q92831"), most repeated in acyclic paths, they act as source
# targets are unique
# a question is raised about these candidates. lists report about its degree
betCent = nx.betweenness_centrality(PPI_graph, normalized=True, endpoints=True)

print("\nsource nodes")
print(f"Degree of Q9NXR8: {PPI_graph.degree('Q9NXR8')}",
      f"Betweenness Centrality of Q92831: {betCent['Q92831']}")
print(f"Degree of Q00535: {PPI_graph.degree('Q00535')}",
      f"Betweenness Centrality of Q00535: {betCent['Q00535']}")
print(f"Degree of Q92831: {PPI_graph.degree('Q92831')}",
      f"Betweenness Centrality of Q92831: {betCent['Q92831']}")
print("target nodes")
print(f"Degree of Q9BY41: {PPI_graph.degree('Q9BY41')}",
      f"Betweenness Centrality of Q9BY41: {betCent['Q9BY41']}")
print(f"Degree of Q8IZL8: {PPI_graph.degree('Q8IZL8')}",
      f"Betweenness Centrality of Q8IZL8: {betCent['Q8IZL8']}")
print(f"Degree of P00533: {PPI_graph.degree('P00533')}",
      f"Betweenness Centrality of P00533: {betCent['P00533']}")


# let's examine the paths that are most rejected, and see the nature of the source nodes. The majority are 4 in length.
most_reject_paths = [("P05388", "O60216"), ("Q15572", "Q9Y383"), ("Q92833", "P31327"),
                     ("P20933", "P40692"), ("P02765", "Q11203"), ("Q9H840",
                                                                  "O14979"), ("P02766", "P49619"),
                     ("P02766", "P00734"), ("P00734", "P62753"), ("P40429",
                                                                  "Q14103"), ("P02768", "Q01850"),
                     ("P02766", "Q9Y3A3"), ("P05387", "Q8NHV4"), ("Q9UBM7",
                                                                  "P36406"), ("Q9H840", "P49757"),
                     ("Q9UBM7", "Q9H1K0"), ("P05387", "Q9NPI1"), ("P05387",
                                                                  "P56134"), ("P02765", "O75496"),
                     ("P05388", "Q9NZI8"), ("P20933", "O60925"), ("P39023",
                                                                  "Q8IXZ2"), ("P02760", "P23396"),
                     ("Q9UBM7", "Q9UKR5"), ("Q15572", "Q9UHJ3")]
# often the longer paths always have a feedback at its beginning


print("\nsource nodes")
print(f"Degree of Q9UBM7: {PPI_graph.degree('Q9UBM7')}",
      f"Betweenness Centrality of Q9UBM7: {betCent['Q9UBM7']}")
print(f"Degree of P02766: {PPI_graph.degree('P02766')}",
      f"Betweenness Centrality of P02766: {betCent['P02766']}")
print(f"Degree of P05387: {PPI_graph.degree('P05387')}",
      f"Betweenness Centrality of P05387: {betCent['P05387']}")
print("target nodes")
print(f"Degree of P62753: {PPI_graph.degree('P62753')}",
      f"Betweenness Centrality of P62753: {betCent['P62753']}")
print(f"Degree of O75496: {PPI_graph.degree('O75496')}",
      f"Betweenness Centrality of O75496: {betCent['O75496']}")
print(f"Degree of Q9UHJ3: {PPI_graph.degree('Q9UHJ3')}",
      f"Betweenness Centrality of Q9UHJ3: {betCent['Q9UHJ3']}")

print("------------------------------Degree Analysis---------------------------------------")


# third point
tail_counter = Counter(df["#tail"])
head_counter = Counter(df["head"])

most_tail_repeated = tail_counter.most_common(1)[0]
most_head_repeated = head_counter.most_common(1)[0]
print(most_tail_repeated, "the most repeated head")
print(most_head_repeated, "the most repeated tail")

export_protein_data(
    PPI_graph, most_head_repeated[0], "Exported_files/protein_data.txt")

print("------------------------------Degree Distribution (free-scale property)---------------------------------------")

# fourth point
protein_set = list(PPI_graph.nodes())
draw_degree_histogram(PPI_graph, protein_set,
                      "Exported_files/degree_report.txt", flag="all")

print("------------------------------Conversion Map---------------------------------------")

# fifth point
uniprot_ids = list(PPI_graph.nodes())[:5]
conversion_map = get_conversion_map(uniprot_ids)
print(f"Sample of conversion map{conversion_map}")


print("debugging efforts")
print(f"Number of UniProt IDs from graph: {len(uniprot_ids)}")
print(f"Number of UniProt IDs in conversion map: {len(conversion_map)}")
missing_data = list()
for protein in uniprot_ids:
    if not protein in conversion_map:
        missing_data.append(protein)
print(f"Missing data list: {missing_data}")

# this protein does not have a name.
print("the genes name of protein 'Q9Y2S0'")
curr_protein = 'Q9Y2S0'
u = UniProt()
entries = u.get_df(curr_protein)
print(list(entries["Gene Names"]))  # the gene name is not known

# sample of Uniport output
print("the genes name of protein 'Q8NI27'")
curr_protein = 'Q8NI27'
u = UniProt()
entries = u.get_df(curr_protein)
print(list(entries["Gene Names"]))  # the gene name is known


# sixth point
file_name = "Exported_files/unweighted_graph.csv"
unweighted_graph = saving_as_unweighed_graph_as_adjMatrix(PPI_graph, file_name)
print("adjacency matrix has been generated")
