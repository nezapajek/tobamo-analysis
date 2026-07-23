# Phylogenetic Placement Protocol

This protocol describes the phylogenetic placement pipeline using EPA-ng. The workflow consists of three main steps:

1. **Optional**: Add reverse complements to contigs
2. Prepare reference alignment and tree
3. Perform phylogenetic placement of contigs

---

## Setup

```bash
# Activate conda environment
conda activate phylo_placement
```

---

## Workflow

### Step 0: Add Reverse Complements (Optional)

If you want to test both orientations of your contigs:

```bash
python scripts/add_rc.py data/contigs.fasta data/contigs_rc.fasta
```

This doubles the number of sequences, adding `rc_` prefix to reverse complement sequences.

---

### Step 1: Prepare Reference Alignment and Tree

Build a multiple sequence alignment and phylogenetic tree from reference sequences.

#### Option A: Without Trimming

```bash
python scripts/prep_refs.py \
  --ref-fasta data/reference.fasta \
  --output-dir results/untrimmed \
  --threads 32
```

#### Option B: With Trimming (recommended for divergent sequences)

```bash
python scripts/prep_refs.py \
  --ref-fasta data/reference.fasta \
  --output-dir results/trimmed \
  --trimal \
  --threads 32
```

**Outputs:**
- `reference.aln` - MAFFT alignment
- `reference_trimmed.aln` - trimAl output (if `--trimal` used)
- `reference.treefile` - IQ-TREE phylogenetic tree
- `reference.iqtree` - IQ-TREE model information
- `preprocess_summary.json` - Summary file for next step

---

### Step 2: Phylogenetic Placement

Place contig sequences onto the reference tree using EPA-ng.

#### Using the summary file (recommended)

```bash
python scripts/phylo_placement.py \
  --from-preprocess results/untrimmed/preprocess_summary.json \
  --contig-fasta data/contigs.fasta \
  --output-dir results/untrimmed/placement
```

#### Or specify paths manually

```bash
python scripts/phylo_placement.py \
  --contig-fasta data/contigs.fasta \
  --ref-alignment results/untrimmed/reference.aln \
  --tree results/untrimmed/reference.treefile \
  --model results/untrimmed/reference.iqtree \
  --output-dir results/untrimmed/placement
```

#### Debug mode (test with first 2 contigs only)

```bash
python scripts/phylo_placement.py \
  --from-preprocess results/untrimmed/preprocess_summary.json \
  --contig-fasta data/contigs.fasta \
  --output-dir results/untrimmed/placement_test \
  --debug
```

**Outputs:**
- `combined_aligned.fasta` - Reference + contig alignment
- `contigs_aligned.fasta` - Aligned contigs only
- `epa_output/epa_result.jplace` - EPA-ng placement results

---

### Step 3: Analyze Placement Results with gappa

Visualize and analyze the phylogenetic placements.

#### View placement information

```bash
gappa examine info --jplace-path results/untrimmed/placement/epa_output/epa_result.jplace
```

#### Generate visualization (grafted tree)

```bash
gappa examine graft \
  --jplace-path results/untrimmed/placement/epa_output/epa_result.jplace \
  --out-dir results/untrimmed/gappa
```

This creates a Newick tree file with placed sequences that can be visualized in tools like FigTree, iTOL, or ggtree.

---

## Complete Example Workflows

### Workflow A: Untrimmed alignment

```bash
# 1. Prepare references
python scripts/prep_refs.py \
  --ref-fasta data/reference.fasta \
  --output-dir results/untrimmed \
  --threads 32

# 2. Place contigs
python scripts/phylo_placement.py \
  --from-preprocess results/untrimmed/preprocess_summary.json \
  --contig-fasta data/contigs.fasta \
  --output-dir results/untrimmed/placement

# 3. Visualize
gappa examine graft \
  --jplace-path results/untrimmed/placement/epa_output/epa_result.jplace \
  --out-dir results/untrimmed/gappa
```

### Workflow B: Trimmed alignment with reverse complements

```bash
# 0. Add reverse complements
python scripts/add_rc.py data/contigs.fasta data/contigs_rc.fasta

# 1. Prepare references with trimming
python scripts/prep_refs.py \
  --ref-fasta data/reference.fasta \
  --output-dir results/trimmed \
  --trimal \
  --threads 32

# 2. Place contigs (including reverse complements)
python scripts/phylo_placement.py \
  --from-preprocess results/trimmed/preprocess_summary.json \
  --contig-fasta data/contigs_rc.fasta \
  --output-dir results/trimmed/placement

# 3. Visualize
gappa examine graft \
  --jplace-path results/trimmed/placement/epa_output/epa_result.jplace \
  --out-dir results/trimmed/gappa
```

---

## Tips

1. **Reference sequences**: Include at least one outgroup sequence for rooting the tree
2. **Threads**: Use `--threads -1` (default) to auto-detect available cores, or specify a number
3. **Trimming**: Use `--trimal` if references are divergent or have many gaps
4. **Debug mode**: Test with `--debug` flag to process only 2 contigs before running full analysis
5. **Redo**: Use `--redo` flag to recompute existing output files

---

## Manual Step-by-Step Commands (Advanced)

If you prefer to run individual tools directly:

### Without Trimming

```bash
# 1. Align references
mafft --auto --thread -1 data/reference.fasta > results/reference.aln

# 2. Build reference tree
iqtree3 -s results/reference.aln -m MFP -bb 1000 -T -1 --prefix results/reference

# 3. Align contigs to reference (add fragments)
mafft --addfragments data/contigs.fasta --keeplength results/reference.aln > results/combined.aln

# 4. Extract only contigs from alignment
# (Manual filtering or use phylo_placement.py)

# 5. Run EPA-ng placement
epa-ng --ref-msa results/reference.aln \
  --tree results/reference.treefile \
  --query results/contigs_aligned.fasta \
  --model results/reference.iqtree \
  -w results/epa_output
```

### With Trimming

```bash
# 1. Align references
mafft --auto --thread -1 data/reference.fasta > results/reference.aln

# 2. Trim alignment
trimal -in results/reference.aln -out results/reference_trimmed.aln -automated1

# 3. Build reference tree
iqtree3 -s results/reference_trimmed.aln -m MFP -bb 1000 -T -1 --prefix results/reference

# 4. Align contigs to trimmed reference
mafft --addfragments data/contigs.fasta --keeplength results/reference_trimmed.aln > results/combined.aln

# 5. Extract only contigs from alignment
# (Manual filtering or use phylo_placement.py)

# 6. Run EPA-ng placement
epa-ng --ref-msa results/reference_trimmed.aln \
  --tree results/reference.treefile \
  --query results/contigs_aligned.fasta \
  --model results/reference.iqtree \
  -w results/epa_output
```
