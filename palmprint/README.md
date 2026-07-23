# Palmprint Analysis

This directory contains tools and workflows for identifying viral palmprints (conserved RNA-dependent RNA polymerase domains) in sequence data using PalmScan.

## Overview

Palmprints are highly conserved structural domains found in viral RNA-dependent RNA polymerases (RdRp). This analysis pipeline:

1. Predicts open reading frames (ORFs) from nucleotide sequences
2. Searches for palmprint domains using position-specific scoring matrices (PSSMs)
3. Extracts and reports sequences containing viral palmprints

## Prerequisites

### Software Requirements

- **Conda/Miniconda**: For installing EMBOSS
- **Docker**: For running PalmScan
- **Git**: For cloning repositories (if needed)

### Installation

1. **Install EMBOSS (for getorf)**:
   ```bash
   conda install bioconda::emboss
   ```

2. **Build PalmScan Docker image**:
   ```bash
   docker build -t palmscan .
   ```

## Usage

### Quick Start

1. **Get help and available options**:
   ```bash
   docker run -it --rm \
     -v "$PWD":/usr/src/palmscan/data \
     -w /usr/src/palmscan/data \
     palmscan -help
   ```

2. **Run PalmScan analysis**:
   ```bash
   docker run -it --rm \
     -v "$PWD":/usr/src/palmscan/data \
     -w /usr/src/palmscan/data \
     palmscan \
     -search_pssms seqs.fasta \
     -tsv hits.tsv \
     -report_pssms report \
     -fasta pp.fa
   ```

### Input Requirements

- **Input file** (`seqs.fasta`): Amino acid sequences in FASTA format
- **File location**: Must be in the current working directory or mounted volume

### Output Files

- **`hits.tsv`**: Tab-separated file containing palmprint hit information
- **`report`**: Detailed report files (generates seqA, seqB, seqC files)
- **`pp.fa`**: FASTA file containing complete palmprint sequences

### Docker Command Breakdown

```bash
docker run -it --rm \
  -v "$PWD":/usr/src/palmscan/data \    # Mount current directory
  -w /usr/src/palmscan/data \           # Set working directory
  palmscan \                            # Docker image name
  -search_pssms seqs.fasta \           # Input amino acid sequences
  -tsv hits.tsv \                      # Output TSV file
  -report_pssms report \               # Output report prefix
  -fasta pp.fa                         # Output palmprint sequences
```

## Complete Workflow

For a complete analysis pipeline including ORF prediction, see the [protocol.md](protocol.md) file.

## File Structure

```
palmprint/
├── README.md           # This file
├── protocol.md         # Complete analysis protocol
├── Dockerfile          # Docker build instructions
├── data/              # Input data directory
├── results/           # Output results directory
│   ├── getorf/        # GETORF output
│   ├── orfipy/        # ORFIPY output
│   └── palmscan/      # PalmScan results
└── scripts/           # Analysis scripts
```

## Troubleshooting

### Common Issues

1. **Docker build fails**:
   - Ensure Docker is running
   - Check internet connection for downloading dependencies
   - Verify Dockerfile syntax

2. **Permission denied errors**:
   - Check file permissions in mounted directories
   - Ensure Docker has access to the working directory

3. **Input file not found**:
   - Verify file paths are correct
   - Ensure files are in the mounted directory (`$PWD`)

4. **Empty output files**:
   - Check input sequence format (should be amino acids, not nucleotides)
   - Verify sequences meet minimum length requirements
   - Check for valid FASTA format

### Getting Help

- Run `palmscan -help` for command-line options
- Check the [protocol.md](protocol.md) for step-by-step instructions
- Review input file formats and requirements

## Notes

- PalmScan requires **amino acid sequences** as input, not nucleotide sequences
- Use GETORF or ORFIPY to translate nucleotide sequences to amino acids first
- The analysis identifies conserved palmprint domains characteristic of viral RdRp
- Results can be used for viral classification and phylogenetic analysis