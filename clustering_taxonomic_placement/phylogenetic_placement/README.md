# Phylogenetic Placement Pipeline

Phylogenetic placement of viral contigs onto reference trees using EPA-ng (Evolutionary Placement Algorithm).

## Overview

This pipeline places query sequences (viral contigs) onto a pre-computed reference phylogenetic tree. It consists of three main steps:

1. **Prepare reference data**: Align reference sequences and build phylogenetic tree
2. **Phylogenetic placement**: Align contigs to reference and place onto tree using EPA-ng
3. **Visualization**: Convert placement results to annotated trees using gappa

## Quick Start

```bash
# Activate environment
conda env create -f environment.yml
conda activate phylo_placement

# 1. Prepare reference alignment and tree
python scripts/prep_refs.py \
  --ref-fasta data/reference.fasta \
  --output-dir results/untrimmed \
  --threads 32

# 2. Place contigs onto reference tree
python scripts/phylo_placement.py \
  --from-preprocess results/untrimmed/preprocess_summary.json \
  --contig-fasta data/contigs.fasta \
  --output-dir results/untrimmed/placement

# 3. Generate visualization
gappa examine graft \
  --jplace-path results/untrimmed/placement/epa_output/epa_result.jplace \
  --out-dir results/untrimmed/gappa
```

## Scripts

- `prep_refs.py`: Align reference sequences with MAFFT, optionally trim with trimAl, and build tree with IQ-TREE
- `phylo_placement.py`: Align contigs to reference and perform phylogenetic placement with EPA-ng
- `add_rc.py`: Add reverse complement sequences to FASTA file (for testing both orientations)

## Important Notes

### FASTA Header Format

FASTA files must have clean headers without special characters. Use simple identifiers without delimiters or whitespace. Periods are generally safe.

### Contig Orientation

Contigs from the Snakemake workflow have unknown orientation (forward or reverse strand). Before phylogenetic placement, ensure all contigs are in the correct forward orientation. Mixed orientations will produce incorrect alignments and unreliable placements.

To test both orientations:

```bash
python scripts/add_rc.py data/contigs.fasta data/contigs_rc.fasta
```

This adds reverse complement sequences with `rc_` prefix, doubling the input size.

## Documentation

- `protocol.md`: Detailed step-by-step protocol with all options and examples
- `environment.yml`: Conda environment specification with all required tools

