# Phylogenetic Placement — Rerun on Orientation-Corrected Contigs

Command log for rerunning phylogenetic placement after contigs found to be
in the wrong orientation were manually corrected to the forward strand.
Otherwise the same pipeline as `protocol.md`, run for both the untrimmed
and trimmed reference trees.

# make a reference multiple alignment and make a tree
python scripts/prep_refs.py \
  --ref-fasta data/reference.fasta \
  --threads 32

# make a reference multiple alignment, trim, and make a tree
python scripts/prep_refs.py \
  --ref-fasta data/reference.fasta \
  --trimal \
  --threads 32

# create query alignments and make phylogenetic placements
# untrimmed ref tree
python scripts/phylo_placement.py \
  --from-preprocess results/untrimmed/preprocess_summary.json \
  --contig-fasta data/contigs.fasta \
  --output-dir results/untrimmed/placement

## trimmed ref tree
python scripts/phylo_placement.py \
  --from-preprocess results/trimmed/preprocess_summary.json \
  --contig-fasta data/contigs.fasta \
  --output-dir results/trimmed/placement

# make a tree with each of the query sequences represented as a pendant edge
# untrimmed ref tree
gappa examine graft \
  --jplace-path results/untrimmed/placement/epa_output/epa_result.jplace \
  --out-dir results/untrimmed/gappa

# trimmed ref tree
gappa examine graft \
  --jplace-path results/trimmed/placement/epa_output/epa_result.jplace \
  --out-dir results/trimmed/gappa
