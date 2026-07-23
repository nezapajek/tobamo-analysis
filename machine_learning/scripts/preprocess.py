#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
from utils import *


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Preprocess data for Tobamovirus classification")

    # Common arguments
    parser.add_argument("refs_path", help="Path to reference Excel file")
    parser.add_argument("orf_fasta_path", help="Path to ORF FASTA file")
    parser.add_argument("pairwise_path", help="Path to pairwise alignment CSV")
    parser.add_argument("output_dir", help="Directory name where output files will be saved in results/")
    parser.add_argument(
        "--contig_fasta_path", help="Path to contig FASTA file (required for non-snakemake contigs)", default=None
    )

    # Data type flags (mutually exclusive)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--train", action="store_true", help="Process training data")
    group.add_argument("--test", action="store_true", help="Process test data")

    # Training-specific arguments
    parser.add_argument("--contigs", help="Path to contigs FASTA file (required for training)", default=None)

    # Test-specific arguments
    parser.add_argument(
        "--snakemake",
        action="store_true",
        help="Use snakemake mode for contigs from snakemake, relevant for parsing names",
    )

    return parser.parse_args()


def preprocess_training(refs_path, contig_fasta_path, orf_fasta_path, pairwise_path, output_dir):
    """
    Process training data
    """
    print("Starting training data preprocessing...")

    # Check if files exist and are in the correct format
    check_file_exists(refs_path, "Excel")
    check_file_exists(contig_fasta_path, "FASTA")
    check_file_exists(pairwise_path, "CSV")
    check_file_exists(orf_fasta_path, "FASTA")

    # Create output directory
    os.makedirs(f"results/{output_dir}", exist_ok=True)

    # Load data
    print("Loading reference data and pairwise alignments")
    refs = pd.read_excel(refs_path)
    pairwise = pd.read_csv(pairwise_path)

    # Filter out parent amino acid refs
    print("Filtering out parent amino acid refs")
    pairwise_minus, removed_refs = filter_pairwise(pairwise, refs)

    # Make mappers from refs
    print("Creating reference mappers")
    aa_id2type, aa_id2prot, inv_type_acc_dict_path, inv_nt_virusname2id_path = make_mappers(refs)

    # Aggregate and pivot
    print("Aggregating and pivoting alignment data")
    agg = aggregate_df(pairwise_minus, aa_id2type, aa_id2prot)
    pvtd = pivot_df(agg)

    # Add basic info
    print("Adding basic sequence information")
    basic = add_info_basic(pvtd, orf_fasta_path, snakemake=False)

    # Add training-specific info
    print("Adding training-specific information")
    input_df = add_info_to_training_input_df(basic, orf_fasta_path, inv_nt_virusname2id_path, inv_type_acc_dict_path)

    # Save results
    print("Saving processed training data")
    input_df.to_csv(f"results/{output_dir}/training_input.csv", index=False)
    removed_refs.to_csv(f"results/{output_dir}/removed_refs.csv", index=False)

    print(f"Training preprocessing complete. Output saved to results/{output_dir}/")
    return input_df


def preprocess_testing(refs_path, orf_fasta_path, pairwise_path, output_dir, snakemake=True, contig_fasta_path=None):
    """
    Process test data
    """
    print("Starting test data preprocessing...")

    # Check if files exist and are in the correct format
    check_file_exists(refs_path, "Excel")
    check_file_exists(pairwise_path, "CSV")
    check_file_exists(orf_fasta_path, "FASTA")

    # Create output directory
    os.makedirs(f"results/{output_dir}", exist_ok=True)

    # Load data
    print("Loading reference data and pairwise alignments")
    refs = pd.read_excel(refs_path)
    pairwise = pd.read_csv(pairwise_path)

    # Make mappers from refs
    print("Creating reference mappers")
    aa_id2type, aa_id2prot, _, _ = make_mappers(refs)

    # Aggregate and pivot
    print("Aggregating and pivoting alignment data")
    agg = aggregate_df(pairwise, aa_id2type, aa_id2prot)
    pvtd = pivot_df(agg)

    # Add basic info
    print("Adding basic sequence information")

    if snakemake:
        input_df = add_info_basic(pvtd, orf_fasta_path, snakemake=True)
    else:
        if contig_fasta_path is None:
            raise ValueError("Contig FASTA path must be provided for non-snakemake mode")
        input_df = add_info_basic(pvtd, orf_fasta_path, snakemake=False, contig_fasta_path=contig_fasta_path)

    # Save results
    print("Saving processed test data")
    input_df.to_csv(f"results/{output_dir}/testing_input.csv", index=False)

    print(f"Test preprocessing complete. Output saved to results/{output_dir}/")
    return input_df


def main():
    """Main function"""
    args = parse_args()

    # Validate arguments
    if args.train and args.contigs is None:
        print("Error: --contigs is required when using --train")
        sys.exit(1)

    # Process data based on specified mode
    if args.train:
        preprocess_training(args.refs_path, args.contigs, args.orf_fasta_path, args.pairwise_path, args.output_dir)
    elif args.test:
        preprocess_testing(
            args.refs_path,
            args.orf_fasta_path,
            args.pairwise_path,
            args.output_dir,
            args.snakemake,
            args.contig_fasta_path,
        )

    print("Preprocessing completed successfully!")


if __name__ == "__main__":
    main()
