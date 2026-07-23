# Clustering and Taxonomic Placement

**Purpose:** Group the curated candidate contigs (from `manual_curation_annotation/`)
into clusters of closely related sequences, then place representative
sequences from each cluster onto a reference tobamovirus phylogeny — the two
halves of the article's "Clustering and taxonomic placement" module (Fig. 1).

## `clustering/`

All-vs-all BLASTN across the curated contigs, filtered to high-identity hits
(≥90% identity, ≥50nt aligned), built into a similarity network with
NetworkX; connected components define clusters. One representative sequence
per cluster feeds phylogenetic analysis. See
[`clustering/README.md`](clustering/README.md) for the full workflow and
parameters.

## `phylogenetic_placement/`

Two related but distinct analyses, both built on the cluster representatives:

1. **Phylogenetic inference** — four separate per-ORF trees (ORF1–ORF4), done
   externally by Lana Vogrinec — **no code for this step lives in this
   repo**. Per ORF: cluster members were aligned with MUSCLE in MEGA11
   (v11.0.13), and one representative per cluster was selected (the sequence
   spanning the longest region aligned across most contigs). Cluster
   representatives were queried against NCBI nr via BLASTX (e-value 0.05,
   word size 5, gap open 11, gap extend 1, matrix BLOSUM62), and the top 10
   unique hits (by E-value) per representative were aligned together with the
   cluster representatives, ICTV VMR reference tobamoviral sequences (MSL40,
   v1.20250307), and an outgroup — Oat golden stripe virus (*Furovirus*) for
   the ORF1/ORF2/ORF3 trees, pepper ringspot virus (*Tobravirus*) for ORF4.
   Sequences that didn't align to the region shared by most representatives
   were discarded, ends were manually trimmed to the longest aligned segment,
   then further automated trimming was applied with trimAl (v1.13,
   `automated1`) — these four trimmed alignments are
   `tobamo-supp-data/supplementary_data/orf{1,2,3,4}_tree_alignment.fasta`
   (Supplementary Files S1–S4). Maximum-likelihood trees were inferred from
   the trimmed alignments in IQ-TREE (v3.0.1), using automatic
   substitution-model selection and bootstrap support from the bootstrap
   splits method (default parameters).
2. **Phylogenetic placement** (the code that *is* here) — MAFFT-aligns query
   contigs against a fixed reference alignment, places them onto a
   pre-computed IQ-TREE reference tree with EPA-ng, then visualizes/annotates
   with gappa + iTOL. See [`phylogenetic_placement/README.md`](phylogenetic_placement/README.md)
   and [`protocol.md`](phylogenetic_placement/protocol.md) for the full
   workflow.

Note: `phylogenetic_placement/scripts/check_contig_orietntation.ipynb` (an
earlier orientation-checking exploration) was not migrated — its inputs
(`data/representative_contigs*.fasta`) were superseded/archived versions no
longer present. `add_rc.py` is the current, documented way to test contig
orientation (see `phylogenetic_placement/README.md`).
