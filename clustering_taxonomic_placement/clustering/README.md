# Contig Clustering Pipeline

This folder contains scripts for clustering contigs obtained from the Snakemake workflow for downstream phylogenetic analysis.

Sequences are clustered by first performing an all-vs-all blastn search to asses pairwise similarity between contigs. Hits that meet the specified criteria for identity and alignment length are used to build a network graph, with clusters defined as connected components within the graph.


## Overview

The clustering pipeline groups similar contigs based on sequence similarity using BLAST-based network analysis. This enables:
- Identification of closely related viral sequences
- Selection of representative sequences for phylogenetic analysis
- Reduction of redundancy in large contig datasets
- Preparation of clustered groups for comparative genomics studies

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)
- [Workflow](#workflow)
- [Output](#output)
- [Parameters](#parameters)

## Installation

### Prerequisites

- BLAST+ (v.2.15.0 or later)
- Python 3.11+

### Option 1: Manual Installation

Install required software and Python packages:

```bash
# Install BLAST+ (system-dependent)
# Ubuntu/Debian:
sudo apt-get install ncbi-blast+

# macOS:
brew install blast

# Install Python packages
pip install pandas networkx matplotlib biopython
```

### Option 2: Conda Installation (Recommended)

```bash
# Create and activate conda environment
conda env create -f clustering_env.yaml
conda activate clustering_env
```

## Usage

### Quick Start

```bash
# 1. Configure input/output paths in scripts
# 2. Run the complete pipeline
./blast.sh
python filter_self_hits.py
python filter_blast_results.py
python clustering.py
```

## Workflow

The clustering pipeline consists of 4 main steps:

### 1. All-vs-All BLASTN Search

**Script:** `blast.sh`

Performs comprehensive sequence similarity search between all contigs.

**Configuration:**
Edit the script to define input/output files:
```bash
FASTA_FILE="path/to/contigs.fasta"     # Input contigs
DB_NAME="contig_db"                    # BLAST database name
OUTPUT_FILE="blast_results.tsv"        # BLASTN output
```

**Usage:**
```bash
./blast.sh
```

**Input:**
- FASTA file with contigs to cluster (`.fasta`)

**Output:**
- BLAST database files (multiple `.nhr`, `.nin`, `.nsq` files)
- Complete BLASTN results (`.tsv`)

### 2. Remove Self-hits and Duplicates

**Script:** `filter_self_hits.py`

Cleans BLASTN results by removing self-alignments and duplicate hits.

**Configuration:**
Edit the script to define:
```python
INPUT_FILE = "blast_results.tsv"        # Raw BLASTN results
OUTPUT_FILE = "blast_cleaned.tsv"       # Cleaned results
```

**Usage:**
```bash
python filter_self_hits.py
```

**Input:**
- Raw BLASTN results (`.tsv`)

**Output:**
- Cleaned BLASTN results (`.tsv`)

### 3. Filter by Similarity Parameters

**Script:** `filter_blast_results.py`

Filters BLASTN hits based on identity and alignment length thresholds.

**Configuration:**
Edit the script to define:
```python
INPUT_FILE = "blast_cleaned.tsv"        # Cleaned BLASTN results
OUTPUT_EDGES = "network_edges.txt"      # Network edges file
OUTPUT_TABLE = "filtered_results.txt"   # Filtered results table

# Filtering parameters
MIN_IDENTITY = 90.0    # Minimum % identity
MIN_LENGTH = 0.8       # Minimum aligned fraction
```

**Usage:**
```bash
python filter_blast_results.py
```

**Input:**
- Cleaned BLASTN results (`.tsv`)

**Output:**
- Network edges file (`.txt`)
- Filtered BLASTN results table (`.txt`)

### 4. Network Clustering

**Script:** `clustering.py`

Creates network graph from filtered hits and extracts connected components as clusters.

**Configuration:**
Edit the script to define:
```python
FASTA_FILE = "contigs.fasta"           # Original contigs
OUTPUT_EDGES = "network_edges.txt"     # Network edges
CLUSTER_FILE = "clusters.tsv"          # Output clusters
```

**Usage:**
```bash
python clustering.py
```

**Input:**
- Original FASTA file with contigs (`.fasta`)
- Network edges file (`.txt`)

**Output:**
- Cluster assignments file (`.tsv`)

## Output

### Main Output Files

- **`clusters.tsv`** - Final cluster assignments with contig IDs and cluster IDs
- **`network_edges.txt`** - Network edges representing sequence similarities
- **`filtered_results.txt`** - Filtered BLASTN results meeting similarity criteria
- **`blast_cleaned.tsv`** - Cleaned BLASTN results (self-hits removed)

### Output Structure

```
clustering/
├── blast_results.tsv          # Raw BLASTN all-vs-all results
├── blast_cleaned.tsv          # Cleaned results (no self-hits)
├── filtered_results.txt       # Results passing similarity filters
├── network_edges.txt          # Network edges for clustering
├── clusters.tsv               # Final cluster assignments
└── contig_db.*                # BLAST database files
```

## Parameters

### Similarity Filtering Parameters

Configure in `filter_blast_results.py`:

- **`MIN_IDENTITY`** - Minimum percentage identity (default: 90.0%)
  - Higher values = stricter clustering
  - Range: 70.0-100.0%

- **`MIN_LENGTH`** - Minimum aligned fraction (default: 0.8)
  - Fraction of contig length that must align
  - Range: 0.5-1.0

**Example configurations:**
```python
# Strict clustering (highly similar sequences)
MIN_IDENTITY = 95.0
MIN_LENGTH = 0.9

# Relaxed clustering (distantly related sequences)
MIN_IDENTITY = 80.0
MIN_LENGTH = 0.6
```

### Performance Considerations

- **Memory usage:** Scales with number of contigs squared
- **Runtime:** BLASTN step is most time-consuming
- **Disk space:** BLAST results can be large for many contigs

**Recommendations:**
- For >10,000 contigs: Consider pre-filtering by length or quality
- Use multiple CPU cores for BLAST step (`-num_threads` parameter)
- Monitor disk space during BLAST database creation

## Contact

For questions and support, contact lana.vogrinec@nib.si