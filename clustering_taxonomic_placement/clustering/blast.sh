#!/usr/bin/env bash

# Set input/output paths (replace input.fasta, blast_db and output_BLAST.tsv with custom filenames)
FASTA_FILE="input.fasta"
DB_NAME="blast_db"
OUTPUT_FILE="BLAST_results.tsv"

# Make BLAST database
makeblastdb -in "$FASTA_FILE" -dbtype nucl -out "$DB_NAME"

# Run BLASTN (all-vs-all)
blastn -query "$FASTA_FILE" \
       -db "$DB_NAME" \
       -out "$OUTPUT_FILE" \
       -outfmt "6 qseqid sseqid pident length qstart qend sstart send"

echo "BLAST search completed. Results saved to $OUTPUT_FILE"
