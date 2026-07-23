#!/usr/bin/env python3
"""
Phylogenetic Placement Script using EPA-ng
Iterates over contigs, performs MAFFT --addfragments, filters alignments,
runs EPA-ng placement, and combines results into a single .jplace file
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path
from Bio import SeqIO


def run_command(cmd, description, capture_stdout=False):
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

    if result.stdout and not capture_stdout:
        print(result.stdout)

    return result


def mafft_add_fragments(contig_record, ref_alignment, contig_name):
    """Add contig fragment to reference alignment using MAFFT and return aligned records in memory."""
    import tempfile

    # Create temporary file for single contig
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp:
        SeqIO.write([contig_record], tmp, "fasta")
        tmp_contig = tmp.name

    try:
        cmd = ["mafft", "--addfragments", tmp_contig, "--keeplength", str(ref_alignment)]
        result = run_command(cmd, f"Adding fragments from {contig_name}", capture_stdout=True)

        # Parse output and return records in memory
        import io

        records = list(SeqIO.parse(io.StringIO(result.stdout), "fasta"))
        return records
    finally:
        # Clean up temporary file
        Path(tmp_contig).unlink(missing_ok=True)


def filter_contig_only(aligned_records, contig_id):
    """Filter alignment records to keep only the contig sequence (remove references)."""
    # Filter to keep only the contig
    contig_records = [rec for rec in aligned_records if rec.id == contig_id]

    if len(contig_records) == 0:
        print(f"WARNING: Contig {contig_id} not found in alignment!")
        return None

    return contig_records[0]


def run_epa_ng(ref_alignment, tree_file, query_fasta, model_file, output_dir, redo=False):
    """Run EPA-ng for phylogenetic placement."""
    jplace_file = output_dir / "epa_result.jplace"

    if not redo and jplace_file.exists():
        print(f"\n[SKIP] EPA-ng result already exists: {jplace_file}")
        return jplace_file

    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "epa-ng",
        "--ref-msa",
        str(ref_alignment),
        "--tree",
        str(tree_file),
        "--query",
        str(query_fasta),
        "--model",
        str(model_file),
        "-w",
        str(output_dir),
        "--redo",
    ]

    run_command(cmd, "Running EPA-ng phylogenetic placement")

    return jplace_file


def parse_args():
    parser = argparse.ArgumentParser(
        description="Phylogenetic placement pipeline using EPA-ng",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Using preprocess summary file
  python phylo_placement.py \\
    --from-preprocess results/preprocess_summary.json \\
    --contig-fasta data/contigs.fasta \\
    --output-dir results/placement

  # Or specify paths manually
  python phylo_placement.py \\
    --contig-fasta data/contigs.fasta \\
    --ref-alignment results/reference.aln \\
    --tree results/reference.treefile \\
    --model results/reference.iqtree \\
    --output-dir results/placement
        """,
    )

    # Required argument
    parser.add_argument(
        "--contig-fasta", required=True, help="Path to contig FASTA file (can contain multiple contigs)"
    )

    # Option 1: Use preprocess summary file
    parser.add_argument(
        "--from-preprocess", help="Path to preprocess_summary.json from prep_refs.py (auto-fills ref paths)"
    )

    # Option 2: Manual path specification (required if not using --from-preprocess)
    parser.add_argument("--ref-alignment", help="Path to reference alignment (from prep_refs.py)")

    parser.add_argument("--tree", help="Path to reference tree file (.treefile from prep_refs.py)")

    parser.add_argument("--model", help="Path to IQ-TREE model file (.iqtree from prep_refs.py)")

    parser.add_argument("--output-dir", required=True, help="Directory for all output files")

    # Optional arguments
    parser.add_argument("--redo", action="store_true", help="Redo all steps even if output files already exist")

    parser.add_argument(
        "--debug", action="store_true", help="Debug mode: process only the first 2 contigs to test the pipeline"
    )

    args = parser.parse_args()

    # Load from preprocess summary if provided
    if args.from_preprocess:
        summary_file = Path(args.from_preprocess)
        if not summary_file.exists():
            print(f"ERROR: Summary file not found: {summary_file}")
            sys.exit(1)

        with open(summary_file, "r") as f:
            summary = json.load(f)

        # Override with summary values if not manually specified
        if not args.ref_alignment:
            args.ref_alignment = summary.get("ref_alignment")
        if not args.tree:
            args.tree = summary.get("tree")
        if not args.model:
            args.model = summary.get("model")

        print(f"Loaded configuration from: {summary_file}")

    # Validate required arguments
    if not args.ref_alignment or not args.tree or not args.model:
        parser.error("Either --from-preprocess or all of (--ref-alignment, --tree, --model) are required")

    return args


def main():
    args = parse_args()

    # Setup paths
    contig_fasta = Path(args.contig_fasta)
    ref_alignment = Path(args.ref_alignment)
    tree_file = Path(args.tree)
    model_file = Path(args.model)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 60)
    print("PHYLOGENETIC PLACEMENT PIPELINE (EPA-ng)")
    print("=" * 60)
    print(f"Contig FASTA: {contig_fasta}")
    print(f"Reference alignment: {ref_alignment}")
    print(f"Tree file: {tree_file}")
    print(f"Model file: {model_file}")
    print(f"Output directory: {output_dir}")
    print(f"Redo existing: {args.redo}")
    print(f"Debug mode: {args.debug}")

    # Read contigs
    print("\nReading contigs...")
    contigs = list(SeqIO.parse(contig_fasta, "fasta"))
    print(f"Found {len(contigs)} contig sequences")

    if len(contigs) == 0:
        print("ERROR: No contigs found in input file!")
        sys.exit(1)

    # Debug mode: limit to first 2 contigs
    if args.debug:
        original_count = len(contigs)
        contigs = contigs[:2]
        print(f"\n[DEBUG MODE] Processing only first {len(contigs)} of {original_count} contigs")
        print(f"Debug contigs: {[c.id for c in contigs]}")

    # Step 1: Process all contigs and collect aligned records
    print(f"\n{'='*60}")
    print(f"Step 1: Aligning contigs with reference")
    print(f"{'='*60}")

    ref_records = None  # Will store reference sequences once
    filtered_contig_records = []

    # Process each contig individually for alignment and filtering
    for idx, contig in enumerate(contigs, start=1):
        contig_id = contig.id
        print(f"\n[{idx}/{len(contigs)}] Processing contig: {contig_id}")

        # MAFFT add fragments (returns records in memory)
        aligned_records = mafft_add_fragments(contig, ref_alignment, contig_id)

        # Store reference sequences only from the first contig
        if ref_records is None:
            # Separate references from contig in the first iteration
            ref_records = [rec for rec in aligned_records if rec.id != contig_id]
            print(f"  Collected {len(ref_records)} reference sequences")

        # Filter to keep only contig sequence for output #2
        filtered_contig = filter_contig_only(aligned_records, contig_id)
        if filtered_contig:
            filtered_contig_records.append(filtered_contig)

    if len(filtered_contig_records) == 0:
        print("\nERROR: No filtered contigs to process!")
        sys.exit(1)

    # Step 2: Write output files
    print(f"\n{'='*60}")
    print(f"Step 2: Writing output files")
    print(f"{'='*60}")

    # Output #1: Combined reference and contig alignment (refs once + all contigs)
    combined_ref_contig_fasta = output_dir / "combined_aligned.fasta"
    combined_records = ref_records + filtered_contig_records
    SeqIO.write(combined_records, combined_ref_contig_fasta, "fasta")
    print(
        f"✓ Output 1: Combined reference+contig alignment ({len(ref_records)} refs + {len(filtered_contig_records)} contigs = {len(combined_records)} total)"
    )
    print(f"  → {combined_ref_contig_fasta}")

    # Output #2: Filtered contigs only (for EPA-ng)
    filtered_contigs_fasta = output_dir / "contigs_aligned.fasta"
    SeqIO.write(filtered_contig_records, filtered_contigs_fasta, "fasta")
    print(f"✓ Output 2: Aligned contigs only ({len(filtered_contig_records)} sequences)")
    print(f"  → {filtered_contigs_fasta}")

    # Step 3: Run EPA-ng once on the combined filtered contigs
    print(f"\n{'='*60}")
    print(f"Step 3: Running EPA-ng phylogenetic placement")
    print(f"{'='*60}")

    epa_output_dir = output_dir / "epa_output"
    final_output = run_epa_ng(
        ref_alignment, tree_file, filtered_contigs_fasta, model_file, epa_output_dir, redo=args.redo
    )
    print(f"✓ Output 3: EPA-ng placement results")
    print(f"  → {final_output}")

    # Summary
    print("\n" + "=" * 60)
    print("PHYLOGENETIC PLACEMENT COMPLETED!")
    print("=" * 60)
    print(f"\nTotal contigs processed: {len(contigs)}")
    print(f"\nOutput files:")
    print(f"  1. Combined ref+contig alignment: {combined_ref_contig_fasta}")
    print(f"  2. Aligned contigs only: {filtered_contigs_fasta}")
    print(f"  3. EPA-ng results: {final_output}")

    print("\nNext steps:")
    print("  - Visualize placements using a .jplace viewer")
    print("  - Analyze placement results with gappa or other tools")


if __name__ == "__main__":
    main()
