import pandas as pd

# Set input/output files (replace BLAST_results.tsv and BLAST_deduplicated.tsv with custom filenames)
INPUT_FILE = "BLAST_results.tsv"
OUTPUT_FILE = "BLAST_deduplicated.tsv"

# Load the BLAST output and add column names
df = pd.read_csv(INPUT_FILE, sep="\t", header=None)
df.columns = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send"]

# Filter out self-hits (query == target)
df = df[df["qseqid"] != df["sseqid"]]

# Deduplicate reciprocal hits
df["pair_key"] = df.apply(lambda row: tuple(sorted([row["qseqid"], row["sseqid"]])), axis=1)
df = df.drop_duplicates(subset="pair_key").drop(columns=["pair_key"])

# Save results
df.to_csv(OUTPUT_FILE, sep="\t", index=False)

print(f"Removed self-hits and deduplicated pairs. Output saved to {OUTPUT_FILE}")
