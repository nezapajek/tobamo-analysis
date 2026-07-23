import pandas as pd

# Set input/output files
INPUT_FILE = "BLAST_deduplicated.tsv"
OUTPUT_EDGES = "edges.txt"
OUTPUT_TABLE = "BLAST_filtered.tsv"

# Filtering parameters
MIN_IDENTITY = 90.0
MIN_LENGTH = 50

# Load BLAST table
df = pd.read_csv(INPUT_FILE, sep="\t")

# Apply filtering
df = df[(df["pident"] >= MIN_IDENTITY) & (df["length"] > MIN_LENGTH)]

# Save edge list
with open(OUTPUT_EDGES, "w") as out:
    for _, row in df.iterrows():
        out.write(f"{row['qseqid']}\t{row['sseqid']}\n")

# Save full filtered table
df.to_csv(OUTPUT_TABLE, sep="\t", index=False)

# Print summary
print(f"Edge list complete: {len(df)} edges saved to '{OUTPUT_EDGES}'. "
      f"Full filtered BLAST table saved to '{OUTPUT_TABLE}'.")
