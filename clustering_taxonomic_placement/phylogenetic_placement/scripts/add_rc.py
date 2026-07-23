#!/usr/bin/env python3
"""
Add reverse complement sequences to a FASTA file.
Each reverse complement sequence is added with 'rc_' prefix in the ID.
"""

from Bio import SeqIO
from Bio.Seq import Seq
import sys
import argparse


def add_reverse_complements(input_fasta, output_fasta):
    """
    Read a FASTA file and write a new one with both original sequences
    and their reverse complements.

    Args:
        input_fasta: Path to input FASTA file
        output_fasta: Path to output FASTA file
    """
    records_with_rc = []

    for record in SeqIO.parse(input_fasta, "fasta"):
        # Add original sequence
        records_with_rc.append(record)

        # Create reverse complement
        rc_record = record.reverse_complement(id=f"rc_{record.id}", description="")
        records_with_rc.append(rc_record)

    # Write all records to output
    SeqIO.write(records_with_rc, output_fasta, "fasta")

    print(f"Processed {len(records_with_rc) // 2} sequences")
    print(f"Output contains {len(records_with_rc)} sequences (original + reverse complement)")
    print(f"Saved to: {output_fasta}")


def main():
    parser = argparse.ArgumentParser(description="Add reverse complement sequences to a FASTA file with 'rc_' prefix")
    parser.add_argument("input", help="Input FASTA file")
    parser.add_argument("output", help="Output FASTA file")

    args = parser.parse_args()

    add_reverse_complements(args.input, args.output)


if __name__ == "__main__":
    main()
