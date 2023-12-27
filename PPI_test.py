from util import *

file_path = "PathLinker_2018_human-ppi-weighted-cap0_75.txt"
df = pd.read_csv(file_path, sep='\t')
df = df[:1000]

# checking for any NaN values in each row and column
print("Check for any missing data: ", df.isna().any().any())

# first point
# plotting just 1000 edge for simplicity and clarity
PPI_graph = visualize_graph(df)


# second point
# these source and target proteins have been deduced from generate_path function.
source = "Q00535"
target = "P14373"
file_path = "acyclic_shortest_paths_test2.pdf"
find_acyclic_K_shortest_paths(PPI_graph, source, target, file_path, 2)

# third point
tail_counter = Counter(df["#tail"])
head_counter = Counter(df["head"])

most_tail_repeated = tail_counter.most_common(1)[0]
most_head_repeated = head_counter.most_common(1)[0]

export_protein_data(PPI_graph, most_head_repeated[0], "protein_data.txt")


# fourth point
protein_set = list(PPI_graph.nodes())
draw_degree_histogram(PPI_graph, protein_set, "degree_report.txt", flag="all")

# fifth point
uniprot_ids = list(PPI_graph.nodes())[:5]
conversion_map = get_conversion_map(uniprot_ids)
print(conversion_map)


# sixth point
file_name = "unweighted_graph.csv"
unweighted_graph = saving_unweighed_graph_as_adjMatrix(PPI_graph, file_name)
