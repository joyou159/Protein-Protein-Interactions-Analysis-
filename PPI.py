import util
file_path = "PathLinker_2018_human-ppi-weighted-cap0_75.txt"
df = pd.read_csv(file_path, sep='\t')
df = df[:3000]
df.head()
