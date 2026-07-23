#!/usr/bin/env python3
"""
Reference Preparation Pipeline
Aligns reference sequences, optionally trims, and builds phylogenetic tree
"""

import argparse
import subprocess
import sys
from pathlib import Path
import json


def run_command(cmd, description):
    """Run a shell command and handle errors."""
    print(f"\n{'='*60}")
    print(f"{description}")
    print(f"{'='*60}")
    print(f"Command: {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"ERROR: {description} failed!")
        print(f"STDERR: {result.stderr}")
        sys.exit(1)

    if result.stdout:
        print(result.stdout)

    return result


def run_mafft(input_fasta, output_fasta, threads=-1, redo=False):
    """Run MAFFT multiple sequence alignment."""
    if not redo and output_fasta.exists():
        print(f"\n[SKIP] MAFFT alignment already exists: {output_fasta}")
        return output_fasta

    cmd = ["mafft", "--auto", "--thread", str(threads), str(input_fasta)]

    result = run_command(cmd, "Running MAFFT alignment")

    # Write output to file
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(output_fasta, "w") as f:
        f.write(result.stdout)

    print(f"Alignment written to: {output_fasta}")
    return output_fasta


def run_trimal(input_fasta, output_fasta, redo=False):
    """Run trimAl to trim alignment."""
    if not redo and output_fasta.exists():
        print(f"\n[SKIP] trimAl output already exists: {output_fasta}")
        return output_fasta

    cmd = ["trimal", "-in", str(input_fasta), "-out", str(output_fasta), "-automated1"]

    run_command(cmd, "Running trimAl")
    return output_fasta


def run_iqtree(input_fasta, output_prefix, threads=-1, redo=False):
    """Run IQ-TREE to build reference tree."""
    treefile = Path(str(output_prefix) + ".treefile")
    iqtree_file = Path(str(output_prefix) + ".iqtree")

    if not redo and treefile.exists() and iqtree_file.exists():
        print(f"\n[SKIP] IQ-TREE output already exists: {treefile}")
        return output_prefix

    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "iqtree3",
        "-s",
        str(input_fasta),
        "-m",
        "MFP",
        "-bb",
        "1000",
        "-T",
        str(threads),
        "--prefix",
        str(output_prefix),
    ]

    run_command(cmd, "Running IQ-TREE")
    return output_prefix


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare reference alignment and phylogenetic tree",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Without trimming (auto output to results/untrimmed)
  python prep_refs.py \\
    --ref-fasta data/reference.fasta \\
    --threads 32

  # With trimming (auto output to results/trimmed)
  python prep_refs.py \\
    --ref-fasta data/reference.fasta \\
    --trimal \\
    --threads 32

  # Custom output directory
  python prep_refs.py \\
    --ref-fasta data/reference.fasta \\
    --output-dir results/custom \\
    --threads 32
        """,
    )

    # Required arguments
    parser.add_argument(
        "--ref-fasta", required=True, help="Path to reference FASTA file (include one outgroup for rooting)"
    )

    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory for all output files (default: results/untrimmed or results/trimmed based on --trimal flag)",
    )

    # Optional flags
    parser.add_argument("--trimal", action="store_true", help="Run trimAl on the alignment")

    parser.add_argument(
        "--threads", type=int, default=-1, help="Number of threads for MAFFT and IQ-TREE (default: -1 = auto-detect)"
    )

    parser.add_argument("--redo", action="store_true", help="Redo all steps even if output files already exist")

    return parser.parse_args()


def main():
    args = parse_args()

    # Setup paths
    # Auto-determine output directory based on --trimal flag if not specified
    if args.output_dir is None:
        output_dir = Path("results/trimmed" if args.trimal else "results/untrimmed")
    else:
        output_dir = Path(args.output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    ref_fasta = Path(args.ref_fasta)

    print("\n" + "=" * 60)
    print("REFERENCE PREPARATION PIPELINE")
    print("=" * 60)
    print(f"Reference FASTA: {ref_fasta}")
    print(f"Output directory: {output_dir}")
    print(f"Run trimAl: {args.trimal}")
    print(f"Threads: {args.threads}")
    print(f"Redo existing: {args.redo}")

    # Step 1: Multiple alignment of references
    alignment_output = output_dir / "reference.aln"
    run_mafft(ref_fasta, alignment_output, threads=args.threads, redo=args.redo)

    # Step 2: trimAl (optional)
    if args.trimal:
        trimmed_output = output_dir / "reference_trimmed.aln"
        run_trimal(alignment_output, trimmed_output, redo=args.redo)
        tree_input = trimmed_output
    else:
        print("\n[SKIP] trimAl disabled")
        tree_input = alignment_output

    # Step 3: Build reference tree with IQ-TREE
    iqtree_dir = output_dir / "iqtree"
    iqtree_prefix = iqtree_dir / "reference"
    run_iqtree(tree_input, iqtree_prefix, threads=args.threads, redo=args.redo)

    # Create a summary file with paths for phylo_placement.py
    summary = {
        "ref_alignment": str(tree_input),
        "tree": str(iqtree_prefix) + ".treefile",
        "model": str(iqtree_prefix) + ".iqtree",
        "trimmed": args.trimal,
    }

    summary_file = output_dir / "preprocess_summary.json"
    with open(summary_file, "w") as f:
        json.dump(summary, f, indent=2)

    # Summary
    print("\n" + "=" * 60)
    print("PIPELINE COMPLETED SUCCESSFULLY!")
    print("=" * 60)
    print(f"\nOutput files:")
    print(f"  - Reference alignment: {alignment_output}")
    if args.trimal:
        print(f"  - Trimmed alignment: {trimmed_output}")
    print(f"  - Tree file: {iqtree_prefix}.treefile")
    print(f"  - Model file: {iqtree_prefix}.iqtree")
    print(f"  - Summary file: {summary_file}")

    print("\nNext steps:")
    print("  - Review the reference tree")
    print("  - Run phylogenetic placement:")
    print(f"    python scripts/phylo_placement.py \\")
    print(f"      --from-preprocess {summary_file} \\")
    print(f"      --contig-fasta data/contigs.fasta \\")
    print(f"      --output-dir {output_dir}/placement")


if __name__ == "__main__":
    main()
