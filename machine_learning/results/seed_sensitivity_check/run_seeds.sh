#!/bin/bash
set -e
cd /home/tobamo/analize/tobamo-analysis/machine_learning
for seed in 43 44; do
  echo "=== seed $seed start: $(date) ==="
  conda run -n tobamo-model python -u scripts/train_model_pipeline.py \
    results/training/training_input.csv \
    ../data/tobamo/reference_database.xlsx \
    results/training/sampling/2025-07-11_sampled_contigs_30.fasta \
    --stage evaluate \
    --iterations 30 \
    --sample_depth 30 \
    --seed "$seed" \
    --n-jobs -1 \
    --outdir "seed_sensitivity_check/seed_$seed"
  echo "=== seed $seed end: $(date) ==="
done
