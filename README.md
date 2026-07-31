# tobamo-analysis

Dataset discovery, manual curation, clustering/phylogenetic placement, and
machine-learning classification code behind:

> Pajek Arambašič N, Bačnik K, Kranjc L, Vogrinec L, Curk T, Kutnjak D.
> **Data mining of global sequence datasets markedly expands known diversity
> within *Tobamovirus* genus.** Manuscript in preparation.

This is the second of three repos this project's code/data was split into:

1. [tobamo-snakemake](https://github.com/nezapajek/tobamo-snakemake) —
   automated assembly and candidate contig discovery (quality control, *de
   novo* assembly, similarity search, preliminary taxonomic assignment).
2. **tobamo-analysis** (this repo) — everything downstream of the Snakemake
   pipeline's output.
3. tobamo-supp-data — supplementary tables and data for the article.

<p align="center">
  <img src="images/pipeline.png" alt=" Flowchart illustrating a five-stage viral discovery pipeline. The process begins with dataset discovery via viral RNA-dependent RNA polymerase mining, followed by automated assembly and contig discovery using Snakemake. Candidate contigs then undergo parallel tracks for manual curation, clustering for phylogenetic analysis, and automated classification using a supervised machine learning model" width="500">
</p>
Schematics were created with draw.io

## Structure

Mirrors the five-module pipeline in the paper's Figure 1 (the first module,
automated assembly, is `tobamo-snakemake`; the other four live here):

```
tobamo-analysis/
├── dataset_discovery/              # Serratus PalmID hits -> 253 selected SRRs
├── manual_curation_annotation/     # Snakemake output -> 510 curated candidate contigs
├── palmprint/                      # viral palmprint (RdRp) domain identification
├── clustering_taxonomic_placement/ # clustering + phylogenetic placement of novel contigs
├── machine_learning/               # two-stage ORF/contig classifier
└── data/                           # shared curated datasets and reference collections
```

- **`dataset_discovery/`** — Serratus PalmID search with 35 ICTV-recognized
  tobamovirus palmprints against all SRA short-read datasets (6,870 hits),
  filtered to `50 < pident < 95` → 253 candidate SRRs.
- **`manual_curation_annotation/`** — narrows the Snakemake pipeline's raw
  output (2,567 contigs) down to 510 curated candidate contigs (control-SRR
  removal, dedup, cellular-organism removal, one anomalous-dataset
  exclusion), then supports domain-expert manual labeling into six
  categories (known/novel/divergent tobamoviral, related/unrelated other
  viruses, misassembled).
- **`palmprint/`** — palmprint (RdRp) domain identification via ORF
  prediction + PalmScan; not a named module in the paper's figure, but
  informative for the "why" behind curation and the ML classifier's feature
  design.
- **`clustering_taxonomic_placement/`** — BLAST+NetworkX clustering of novel
  tobamoviral contigs, then MAFFT+EPA-ng phylogenetic placement onto a fixed
  reference tree (external MEGA11/MUSCLE/trimAl/IQ-TREE tree-building step
  not included — done by another researcher).
- **`machine_learning/`** — two-stage supervised classifier (Random Forest
  ORF classifier → Logistic Regression contig-level aggregation via stacked
  generalization) that automates curation of new candidate contigs; ~0.974
  cross-validation accuracy, ~0.899 AUC on this study's discovered contigs.
- **`data/`** — shared curated datasets (contigs, reference sequences, domain
  scientist-provided ground truth/annotations) consumed by the modules above.

Each subfolder has its own `README.md` with setup/usage details.

## License

MIT — see `LICENSE`.
