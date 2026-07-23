import networkx as nx
import matplotlib.pyplot as plt
from Bio import SeqIO

# Set input/output files (replace input.fasta, edges.txt and clusters.tsv with custom filenames)
FASTA_FILE = "input.fasta"
EDGE_FILE = "edges.txt"
CLUSTER_FILE = "clusters.tsv"

# Build graph
G = nx.Graph()

# Add nodes from FASTA_FILE (contigs)
all_contigs = [record.id for record in SeqIO.parse(FASTA_FILE, "fasta")]
for contig in all_contigs:
    G.add_node(contig)

print(f"Added {len(all_contigs)} contigs as nodes.")

# Add edges from EDGE_FILE
with open(EDGE_FILE, "r") as f:
    for line in f:
        query, target = line.strip().split("\t")
        G.add_edge(query, target)

print(f"Added edges from '{EDGE_FILE}'.")

# Extract clusters (connected components)
clusters = list(nx.connected_components(G))

# Write clusters
with open(CLUSTER_FILE, "w") as out_map:
    out_map.write("fragment_id\tcluster_id\n")
    for i, cluster in enumerate(clusters, 1):
        cluster_id = f"Cluster_{i}"
        for frag in cluster:
            out_map.write(f"{frag}\t{cluster_id}\n")

print(f"Found {len(clusters)} clusters. Clusters saved to '{CLUSTER_FILE}'.")