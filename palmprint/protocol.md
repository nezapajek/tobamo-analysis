# Palmprint Analysis Protocol

This protocol describes the process for identifying viral palmprints (RNA-dependent RNA polymerase conserved domains) in contig sequences using ORF prediction and PalmScan analysis.

## Prerequisites

- Docker installed and running
- Input file: `../data/contigs/contigs_non_cellular_filtered.fasta`
- PalmScan Docker image available

## Setup

First, navigate to the palmprint analysis directory and create the necessary output directories:

```bash
cd tobamo-analysis/palmprint
mkdir -p results/palmscan
mkdir -p results/getorf
mkdir -p results/orfipy
```

## Step 1: ORF (Open Reading Frame) Prediction

We use two different approaches to predict ORFs from the input sequences:

### Option A: GETORF (EMBOSS)

**Method 1: Find first ORF only (-find 1)**
```bash
getorf -sequence ../data/contigs/contigs_non_cellular_filtered.fasta \
       -outseq results/getorf/getorf_output_find1.fasta \
       -table 0 \
       -find 1 \
       -minsize 270
```

**Method 2: Find all ORFs between stops (-find 0)**
```bash
getorf -sequence ../data/contigs/contigs_non_cellular_filtered.fasta \
       -outseq results/getorf/getorf_output_find0.fasta \
       -table 0 \
       -find 0 \
       -minsize 270
```

### Option B: ORFIPY

```bash
orfipy ../data/contigs/contigs_non_cellular_filtered.fasta \
       --min 270 \
       --pep PEP \
       --bed BED \
       --between-stops \
       --outdir results/orfipy
```

> **Note:** `getorf -find 0` is equivalent to `orfipy --between-stops`

## Step 2: PalmScan Analysis

Run PalmScan on the predicted ORFs to identify viral palmprints:

### Analysis 1: PalmScan with GETORF (-find 1)

```bash
docker run -it --rm \
    -v "$PWD":/usr/src/palmscan/data \
    -w /usr/src/palmscan/data \
    palmscan \
    -search_pssms results/getorf/getorf_output_find1.fasta \
    -tsv results/palmscan/palmscan_hits_find1.tsv \
    -report_pssms results/palmscan/palmscan_report_find1 \
    -fasta results/palmscan/palmscan_pp_find1.fa
```

### Analysis 2: PalmScan with GETORF (-find 0) or ORFIPY

```bash
docker run -it --rm \
    -v "$PWD":/usr/src/palmscan/data \
    -w /usr/src/palmscan/data \
    palmscan \
    -search_pssms results/getorf/getorf_output_find0.fasta \
    -tsv results/palmscan/palmscan_hits_find0.tsv \
    -report_pssms results/palmscan/palmscan_report_find0 \
    -fasta results/palmscan/palmscan_pp_find0.fa
```

## Output Files

After running the protocol, you will have:

- **TSV files**: Tab-separated values with palmprint hit information
- **Report files**: Detailed PSSM search reports
- **FASTA files**: Sequences containing identified palmprints

## Parameters Explained

- `--min 270` / `-minsize 270`: Minimum ORF size of 270 nucleotides (90 amino acids)
- `-table 0`: Use standard genetic code
- `-find 1`: Find only the first ORF in each sequence
- `-find 0`: Find all ORFs between stop codons
- `--between-stops`: Equivalent to `-find 0`, finds ORFs between stop codons

## Troubleshooting

- Ensure Docker is running before executing PalmScan commands
- Check that input files exist in the specified paths
- Verify that output directories have write permissions