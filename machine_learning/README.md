# *Tobamovirus Contig Classification with Random Forest and Logistic Regression*

This repository provides Python scripts for training and using a **Random Forest Classifier** and **Logistic Regression** to predict **Tobamoviral sequences**. The workflow includes simulating training data, training the model, and using it to classify query contigs. The process ensures a realistic approach to handling sequencing data and contig assembly.

---

## **Table of Contents**
1. [**Installation**](#1-installation)
2. [**Workflow Overview**](#2-workflow-overview)
3. [**Training the Model**](#3-training-the-model)
   - 3.1 [Fitting the curve on Snakemake output data](#31-fitting-the-curve-on-snakemake-output-data-for-weighted-random-sampling)
   - 3.2 [Simulating Sequencing and Assembly for Training Data](#32-simulating-sequencing-and-assembly-for-training-data)
   - 3.3 [Finding ORFs and Pairwise Alignment](#33-finding-orfs-and-pairwise-alignment)
   - 3.4 [Data Processing and Training Input Generation](#34-data-processing-and-training-input-generation)
   - 3.5 [Model Training and Evaluation](#35-model-training-and-evaluation)
4. [**Using the Model for Classification**](#4-using-the-model-for-classification)
   - 4.1 [Preprocessing Query Contigs](#41-preprocessing-query-contigs)
   - 4.2 [Predicting Tobamovirus Contigs](#42-predicting-tobamovirus-contigs)
5. [**Example Usage**](#5-example-usage)
   - 5.1 [Complete Training Workflow Example](#51-complete-training-workflow-example)
   - 5.2 [External Validation](#52-external-validation)
     - 5.2.1 [This Study's Candidate Contigs (Snakemake Output)](#521-this-studys-candidate-contigs-snakemake-output)
     - 5.2.2 [Viral Sequences Outside Virgaviridae](#522-viral-sequences-outside-virgaviridae)
   - 5.3 [Expected Directory Structure After Training](#53-expected-directory-structure-after-training)
6. [**Notes**](#6-notes)
7. [**Contact**](#7-contact)

---

# **1. Installation**
Before running the scripts, you need to set up the environment and install the required tools. We recommend using **Conda** to manage dependencies.  

Run the following commands to set up the environment:  

```bash
# Create and activate the environment
conda env create -f analysis_environment.yml
conda activate tobamo-model
```
Note: You can rename environment in analysis_environment.yml file.

---

# **2. Workflow Overview**

The workflow is divided into two main parts:
  - **Training the Model** (Section 3)
  - **Using the Model for Classification** (Section 4)

---

# **3. Training the Model**

This process involves several steps. Each step corresponds to a specific script that performs part of the pipeline.

## **3.1 Fitting the curve on Snakemake output data for weighted random sampling**

First we take a look inside the Snakemake pipeline output. This part demands some manual checkup and needs to be tailored for each study. In our case, we removed contigs that had hits on *cellular organisms* and contigs from *SRR6846476*, after consulting domain scientists. We then fitted a curve on contig length distribution of the selected contigs, which we'll later use for random weighted sampling of reference genomes, to generate training data.

Example implementation is available in [`notebooks/fit_distribution_curve.ipynb`](notebooks/fit_distribution_curve.ipynb).

## **3.2 Simulating Sequencing and Assembly for Training Data**

This step fragments reference genomes to generate contigs with realistic lengths, ensuring that the training data resembles actual sequencing data.

**Command**
```bash
python scripts/sample_refs.py <path/to/reference.fasta> <out_dir_name> <sampling_num> <subsampling_num> <path/to/lens_freq.json>
```
**Arguments**
| Argument                | Description                                              |
|-------------------------|----------------------------------------------------------|
| `<path/to/reference.fasta>` | Path to the reference FASTA file.                       |
| `<out_dir_name>`        | Directory where the output files will be saved.           |
| `<sampling_num>`        | Total number of samples to generate.                      |
| `<subsampling_num>`     | Number of subsamples per reference sequence.              |
| `<path/to/lens_freq.json>`| Path to the lens_freq.json (Step 1 output)              |

## **3.3 Finding ORFs and Pairwise Alignment**

This step identifies **Open Reading Frames (ORFs)** in contigs and performs **pairwise alignment** with reference proteins. It uses **Orfipy**  and **biopython Bio.Seq.Seq.translate** method to detect ORFs and **MAFFT** to perform pairwise alignments against known reference sequences, such as RdRp ORF1, RdRp ORF2, and Coat Protein from species within the family *Virgaviridae*. 

**Command**
```bash
python getorfs_pairwise_aln.py <path/to/contig.fasta> <out_dir_name> <contig_orientation>
```
**Arguments**
| Argument                  | Description                                                     |
|---------------------------|-----------------------------------------------------------------|
| `<path/to/contig.fasta>`   | Path to the contig FASTA file to be processed.                  |
| `<out_dir_name>`           | Directory where the output files will be saved.                 |
| `<contig_orientation>`     | Orientation of the contigs (e.g., forward or unknown).|

note: if you want to use different reference proteins, change reference database and reference fasta files accordingly (see ../data/tobamo/*).


## **3.4 Data Processing and Training Input Generation**

This step processes reference data and pairwise alignment results to create a training input dataset. It performs data filtering, aggregation, and enrichment with additional sequence information to prepare the final training dataset.

**Command**
```bash
python scripts/preprocess.py <path/to/reference_database.xlsx> <path/to/orf.fasta> <path/to/pairwise_aln.csv> <output_dir> --train --contigs <path/to/sampled_contigs.fasta>
```

**Arguments**
| Argument                           | Description                                                     |
|------------------------------------|-----------------------------------------------------------------|
| `<path/to/reference_database.xlsx>`| Path to the Excel file containing reference protein data.       |
| `<path/to/orf.fasta>`              | Path to the ORF FASTA file from previous step.                 |
| `<path/to/pairwise_aln.csv>`       | Path to the CSV file containing pairwise alignment results.     |
| `<output_dir>`                     | Directory name where output files will be saved in results/.    |
| `--train`                          | Flag to specify processing for training data.                  |
| `--contigs`                        | Path to the contigs FASTA file (required for training).        |

**Processing Steps**
1. **File Validation**: Checks that all input files exist and are in correct formats (Excel, CSV, FASTA)
2. **Data Loading**: Loads reference data from Excel and pairwise alignment results from CSV
3. **Filtering**: Removes parent amino acid references from pairwise data (training mode only)
4. **Mapper Creation**: Creates mappings between amino acid IDs, protein types, and virus names
5. **Data Aggregation**: Aggregates pairwise alignment data by protein type and ID
6. **Data Pivoting**: Restructures data for model input format
7. **Information Enrichment**: Adds basic sequence information and metadata
8. **Output Generation**: Saves processed input data

**Output Files**
- `results/<output_dir>/training_input.csv` - Final training dataset (when using --train)


## **3.5 Model Training and Evaluation**

Our machine learning pipeline employs a three-stage approach to ensure optimal classification performance.

### **3.5.1 Model Selection**

First, we conduct comprehensive algorithm selection through systematic evaluation of multiple classification models:

- **Algorithms Tested**: Logistic Regression, Random Forest, SVM, Decision Tree, Naive Bayes, and KNN
- **Hyperparameter Optimization**: Grid search across all relevant parameters for each algorithm
- **Evaluation Method**: 5-fold cross-validation with stratified sampling
- **Results**: Random Forest demonstrated superior performance across evaluation metrics

### **3.5.2 Cross-Validation Evaluation**

We then compare two different strategies for contig-level prediction:

- **ORF Prediction**: All models use Random Forest for Open Reading Frame (ORF) classification (best performing model from previous step)
- **Contig Prediction Methods**:
  - **Extreme Method**: Uses the most confident ORF prediction score
    - **Binned Prediction Method (stacked generalization)**: Bins ORF predictions and uses Logistic Regression
- **Validation Process**: 5 fold Cross-Validation repeated 30 times
- **Performance Assessment**: Comprehensive metrics including accuracy, F1 score, precision, and recall
- **Winner**: Binned prediction approach (bins=10) achieved superior performance

### **3.5.3 Final Model Training**

The production model combines the best components from our evaluation:

- **ORF Classifier (Morf)**: Random Forest trained on all available training data
- **Contig Classifier (Mc)**: Logistic Regression using binned ORF prediction probabilities
- **Feature Importance**: Analysis reveals most informative sequence characteristics
- **Serialization**: Models saved as joblib files for deployment in production pipeline

This multi-stage approach ensures robust performance across diverse viral sequence data while maintaining interpretability

**Command**
```bash
python scripts/train_model_pipeline.py <path/to/training_input.csv> <path/to/references.xlsx> <path/to/contigs.fasta> --stage <stage> [options]
```

**Arguments**
| Argument                           | Description                                                     |
|------------------------------------|-----------------------------------------------------------------|
| `<path/to/training_input.csv>`     | Path to the training input CSV from Step 4.                     |
| `<path/to/references.xlsx>`        | Path to the Excel file containing reference protein data.       |
| `<path/to/contigs.fasta>`          | Path to the contigs FASTA file.                                |
| `--stage`                          | Pipeline stage to run (`select`, `evaluate`, or `final`).       |
| `--outdir`                         | Output directory name (default: "default").                     |
| `--iterations`                     | Number of iterations for cross-validation (default: 30).        |
| `--sample_depth`                   | Number of contigs to sample per species (default: 30).          |
| `--seed`                           | Random seed for reproducibility (default: 42).                  |
| `--n-jobs`                         | Number of parallel workers for model fitting (`-1` uses all cores; default: 2). |
| `--bin-num`                        | Number of bins for the final model (`--stage final`; default: 10). |
| `--threshold`                      | Custom classification threshold (used in evaluate/final; overrides CV threshold in evaluate when combined with `--use_fixed_threshold`; default behavior uses CV-optimized thresholds in evaluate and 0.5 in final). |
| `--use_fixed_threshold`            | In `--stage evaluate`, use a fixed threshold (`--threshold` or 0.5) instead of CV-optimized thresholds. |
| `--bins`                           | Bin counts to evaluate in `--stage evaluate` (space-separated list, e.g. `--bins 5 10 15`; default: 10). |

**Pipeline Stages**

1. **ORF Model Selection** (`--stage select`)
   - Performs grid search to find the best-performing model type for ORF classification
   - Tests multiple models (RandomForest, LogisticRegression, SVM, etc.)
   - Evaluates models using 5-fold cross-validation
   - Saves performance metrics for each model
   
   ```bash
   python scripts/train_model_pipeline.py training_input.csv references.xlsx contigs.fasta --stage select --outdir model_selection
   ```

2. **Cross-Validation Evaluation** (`--stage evaluate`)
   - Performs extensive evaluation using multiple iterations of 5-fold cross-validation
     - Compares two prediction methods for contig classification:
         - Extreme prediction (most extreme ORF probability)
         - Binned prediction (logistic regression on binned ORF probabilities, bins=10)
     - Threshold behavior for binned predictions:
         - Default: CV-optimized threshold (tuned per fold)
         - Optional fixed threshold: add `--use_fixed_threshold --threshold 0.5`
   - Generates comprehensive performance metrics
   
   ```bash
   python scripts/train_model_pipeline.py training_input.csv references.xlsx contigs.fasta --stage evaluate --iterations 30 --sample_depth 30 --outdir evaluation_results
    ```

    Example with fixed threshold:

    ```bash
    python scripts/train_model_pipeline.py training_input.csv references.xlsx contigs.fasta --stage evaluate --iterations 30 --sample_depth 30 --use_fixed_threshold --threshold 0.5 --outdir evaluation_results
    ```

3. **Final Model Training** (`--stage final`)
   - Trains the final production model on all training data
    - Uses the best-performing approach (Random Forest for ORF classification + Logistic Regression binned prediction for contig-level classification)
   - Saves all necessary model files for deployment
   
   ```bash
   python scripts/train_model_pipeline.py training_input.csv references.xlsx contigs.fasta --stage final --outdir final_model
   ```

**Output Files**
- **Model Selection**:
    - `results/<outdir>/all_performance_metrics.csv` - Combined performance metrics across all iterations/folds
    - `results/<outdir>/model_selection_summary.txt` - Summary of model rankings and best-performing model

- **Cross-Validation Evaluation**:
    - `results/<outdir>/extreme_predictions_results.csv` - Results using extreme prediction
    - `results/<outdir>/binned_predictions_results.csv` - Binned prediction output (default)
    - `results/<outdir>/binned_predictions_results_tuned.csv` - Binned prediction output with CV-optimized thresholds
    - `results/<outdir>/binned_predictions_results_fixed_<threshold>.csv` - Binned prediction output with fixed threshold (e.g. `fixed_0p5`)
    - `results/<outdir>/method_comparison_stats.csv` - Detailed performance comparison between methods
    - `results/<outdir>/method_comparison_simplified.csv` - Mean ± std summary comparison between methods
    - `results/<outdir>/best_method.txt` - Information about the best-performing method

- **Final Model**:
    - `results/<outdir>/rf_model.joblib` - Trained Random Forest model
    - `results/<outdir>/rf_scaler.joblib` - StandardScaler for feature normalization
    - `results/<outdir>/rf_feature_names.csv` - Feature names used by the model
    - `results/<outdir>/lr_binned_10_model.joblib` - Trained Logistic Regression binned prediction model
    - `results/<outdir>/feature_importances.csv` - All feature importances ranked
    - `results/<outdir>/top_40_features.csv` - Top 40 features

**Analyzing Feature Importance**

`scripts/analyze_feature_importance.py` summarizes the Random Forest feature
importances (grouped by source/protein/feature family, plus top-k cumulative
importance):

```bash
python scripts/analyze_feature_importance.py results/final_model/all_features.csv --outdir results/final_model
```

---

# **4. Using the Model for Classification**

Once the model has been trained and finalized, you can use it to classify new query contigs.

## **4.1 Preprocessing Query Contigs**

You must process query contigs through ORF detection, pairwise alignment, and training input preparation. These steps replicate parts of the training pipeline.

### **4.1.1 ORF Detection and Pairwise Alignment**

```bash
# 1. ORF detection and alignment
python scripts/getorfs_pairwise_aln.py ../data/contigs/contigs_all_deduplicated.fasta <out_dir_name> <contig_orientation>
```
**Output File**
- `results/<output_dir>/pairwise_aln.csv` - pairwise alignment metrics table of test data

### **4.1.2 Preprocessing of test data**

Preprocessing of test data (same as training, but with different name parsing and without ground truth)

```bash
# 2. Generate processed input for prediction
python scripts/preprocess.py <path/to/reference_database.xlsx> <path/to/orf.fasta> <path/to/pairwise_aln.csv> <output_dir> --test
```

**Arguments for Test Processing**
| Argument                           | Description                                                     |
|------------------------------------|-----------------------------------------------------------------|
| `<path/to/reference_database.xlsx>`| Path to the Excel file containing reference protein data.       |
| `<path/to/orf.fasta>`              | Path to the ORF FASTA file from previous step.                 |
| `<path/to/pairwise_aln.csv>`       | Path to the CSV file containing pairwise alignment results.     |
| `<output_dir>`                     | Directory name where output files will be saved in results/.    |
| `--test`                           | Flag to specify processing for test/query data.                |

**Output File**
- `results/<output_dir>/testing_input.csv` - Processed data for prediction

## **4.2 Predicting Tobamovirus Contigs**

Once you have the processed input features and the trained models, run the prediction script:

```bash
python scripts/predict_query_contigs.py results/<output_dir>/testing_input.csv results/final_model --outdir predictions --bin-num 10
```

**Arguments**
| Argument                     | Description                                                                 |
|------------------------------|-----------------------------------------------------------------------------|
| `<testing_input_df.csv>`     | Path to the processed input CSV file from preprocessing step.              |
| `<model_dir>`                | Directory containing all model files (RF model, scaler, LR model, etc.).    |
| `--outdir`                   | Name of output directory for prediction results (default: "predictions").   |
| `--bin-num`                  | Number of bins for binned prediction approach (default: 10).               |

The script expects these files in the model directory with standard names:
- `rf_model.joblib` - Trained Random Forest model
- `rf_scaler.joblib` - StandardScaler for feature normalization
- `rf_feature_names.csv` - Feature names used by the model
- `lr_binned_<bin-num>_model.joblib` - Trained Logistic Regression binned prediction model

**Output Files**
After running the script, the following files will be saved in `results/<outdir>/`:

| Filename                   | Description                                                  |
|----------------------------|--------------------------------------------------------------|
| `orf_predictions.csv`      | ORF-level predictions with probability scores.               |
| `contig_predictions.csv`   | Final contig-level predictions with class and probability.   |

---

# **5. Example Usage & External Validation**

This section provides complete examples showing how to run the entire workflow from start to finish.

## **5.1 Complete Training Workflow Example**

Here's a step-by-step example of training a model from scratch:

### **Step 1: Fit Distribution Curve**
```bash
# Run the notebook to fit curve on Snakemake output data
# Input: Snakemake pipeline output data
# Output: results/training/sampling/fitted_curve_lens_freq.json
jupyter notebook notebooks/01_fit_distribution_curve.ipynb
```

### **Step 2: Sample Reference Genomes**
```bash
# Generate training contigs from reference genomes
python scripts/sample_refs.py \
    ../data/tobamo/reference_nukleotidne.fasta \
    training \
    300 \
    30 \
    results/training/sampling/fitted_curve_lens_freq.json
```

### **Step 3: ORF Detection and Pairwise Alignment**
```bash
# Process sampled contigs to find ORFs and perform alignments
time python scripts/getorfs_pairwise_aln.py \
    results/training/sampling/2025-07-11_sampled_contigs_30.fasta \
    training \
    unknown
```

### **Step 4: Preprocess Training Data**
```bash
# Generate training input features
python scripts/preprocess.py \
    ../data/tobamo/reference_database.xlsx \
    results/training/orfs/combined_orfs.fasta \
    results/training/pairwise_aln.csv \
    training \
    --train \
    --contigs results/training/sampling/2025-07-11_sampled_contigs_30.fasta
```

### **Step 5: Model Training (3-stage process)**
```bash
# Stage 1: Model selection via grid search
python scripts/train_model_pipeline.py \
    results/training/training_input.csv \
    ../data/tobamo/reference_database.xlsx \
    results/training/sampling/2025-07-11_sampled_contigs_30.fasta \
    --stage select \
    --outdir model_selection

# Stage 2: Cross-validation evaluation (30 iterations of 5-fold CV)
python scripts/train_model_pipeline.py \
    results/training/training_input.csv \
    ../data/tobamo/reference_database.xlsx \
    results/training/sampling/2025-07-11_sampled_contigs_30.fasta \
    --stage evaluate \
    --iterations 30 \
    --sample_depth 30 \
    --outdir evaluation_results

# Stage 2 alternative: fixed threshold evaluation (recommended for parity with final model threshold=0.5)
python scripts/train_model_pipeline.py \
    results/training/training_input.csv \
    ../data/tobamo/reference_database.xlsx \
    results/training/sampling/2025-07-11_sampled_contigs_30.fasta \
    --stage evaluate \
    --iterations 30 \
    --sample_depth 30 \
    --use_fixed_threshold \
    --threshold 0.5 \
    --outdir evaluation_results_fixed05

# Stage 3: Train final model using RF + LR binned prediction with bins=10
python scripts/train_model_pipeline.py \
    results/training/training_input.csv \
    ../data/tobamo/reference_database.xlsx \
    results/training/sampling/2025-07-11_sampled_contigs_30.fasta \
    --stage final \
    --outdir final_model
```

## **5.2 External Validation**

This mirrors the **External validation** stage of the pipeline schematic
(see [`images/pipeline.png`](../images/pipeline.png)): the finalized model
(Section 3.5.3) is applied to two contig sets it never saw during training
or model selection, neither of which is drawn from the Virgaviridae
reference genomes used to simulate training data.

1. **This study's candidate contigs** (Section 5.2.1) — the actual
   production use case: classifying the Snakemake pipeline's real output
   for this study.
2. **Sanity check: viral sequences outside *Virgaviridae*** (Section 5.2.2)
   — divergent viral sequences that are compositionally similar to the
   training data (real viral ORFs/contigs) but belong to neither the
   tobamoviral nor the non-viral/outgroup class the model was trained to
   separate, used to probe specificity on out-of-distribution input.

### **5.2.1 This Study's Candidate Contigs (Snakemake Output)**

Classifying the 510 curated candidate contigs discovered by this study's
`tobamo-snakemake` pipeline (see the root [`README.md`](../README.md) and
[`manual_curation_annotation/`](../manual_curation_annotation)):

**Preprocessing Query Contigs**
```bash
# Step 1: ORF detection and pairwise alignment for test data
time python scripts/getorfs_pairwise_aln.py \
    ../data/contigs/contigs_all_deduplicated.fasta \
    curated_candidate_contigs \
    unknown

# Step 2: Preprocess test data
python scripts/preprocess.py \
    ../data/tobamo/reference_database.xlsx \
    results/curated_candidate_contigs/orfs/combined_orfs.fasta \
    results/curated_candidate_contigs/pairwise_aln.csv \
    curated_candidate_contigs \
    --test \
    --snakemake
```

**Filter Data (Optional)**
```bash
# Optional: Filter out non-target contigs (keep only 510 non-cellular contigs)
# Use notebook: notebooks/filter_snakemake_pairwise_results.ipynb
```
Note: this notebook does not just filter rows — it re-runs the
aggregation/pivoting step (`aggregate_df` + `pivot_df` in `scripts/utils.py`)
on the filtered subset, so its output
(`pairwise_aln_all_deduplicated_non_cellular_filtered.csv`) is already in
the same wide, per-ORF-feature format as `testing_input.csv`, just for fewer
contigs. Feed it straight to `predict_query_contigs.py` below — **do not**
pass it back through `preprocess.py`, which expects the raw long-format
pairwise alignment table and will fail with `KeyError: 'ref_name'` on an
already-aggregated file.

**Make Predictions**
```bash
# Predict tobamovirus contigs using trained model
python scripts/predict_query_contigs.py \
    results/curated_candidate_contigs/pairwise_aln_all_deduplicated_non_cellular_filtered.csv \
    results/final_model \
    --outdir curated_candidate_contigs
```

### **5.2.2 Viral Sequences Outside Virgaviridae**

Testing the model with divergent viral sequences (outside *Virgaviridae*,
not seen during training and not fitting either the tobamoviral or the
non-viral/outgroup class) to sanity-check specificity on data the model was
never meant to classify confidently one way or the other:

```bash
# Step 1: Preprocess non-Virgaviridae viral sequences
time python scripts/getorfs_pairwise_aln.py \
    ../data/non-virga_representatives/non-virga_tpdb2_diamond_selected.fasta \
    viral_sequences_outside_virgaviridae \
    unknown

python scripts/preprocess.py \
    ../data/tobamo/reference_database.xlsx \
    results/viral_sequences_outside_virgaviridae/orfs/combined_orfs.fasta \
    results/viral_sequences_outside_virgaviridae/pairwise_aln.csv \
    viral_sequences_outside_virgaviridae \
    --test \
    --contig_fasta_path ../data/non-virga_representatives/non-virga_tpdb2_diamond_selected.fasta

# Step 2: Make predictions on non-Virgaviridae sequences
python scripts/predict_query_contigs.py \
    results/viral_sequences_outside_virgaviridae/testing_input.csv \
    results/final_model \
    --outdir viral_sequences_outside_virgaviridae
```

## **5.3 Expected Directory Structure After Training**

After completing the training workflow, your directory structure should look like:

```
results/
├── training/
│   ├── sampling/
│   │   ├── fitted_curve_lens_freq.json
│   │   └── 2025-07-11_sampled_contigs_30.fasta
│   ├── orfs/
│   │   └── combined_orfs.fasta
│   ├── pairwise_aln.csv
│   └── training_input.csv
├── model_selection/
│   ├── all_performance_metrics.csv
│   └── model_selection_summary.txt
├── evaluation_results/
│   ├── extreme_predictions_results.csv
│   ├── binned_predictions_results.csv
│   ├── binned_predictions_results_tuned.csv
│   └── method_comparison_stats.csv
├── final_model/
│   ├── rf_model.joblib
│   ├── rf_scaler.joblib
│   ├── rf_feature_names.csv
│   ├── lr_binned_10_model.joblib
│   ├── feature_importances.csv
│   └── top_40_features.csv
├── curated_candidate_contigs/
│   ├── orfs/
│   │   └── combined_orfs.fasta
│   ├── pairwise_aln.csv
│   ├── testing_input.csv
│   ├── orf_predictions.csv
│   └── contig_predictions.csv
└── viral_sequences_outside_virgaviridae/
    ├── orfs/
    │   └── combined_orfs.fasta
    ├── pairwise_aln.csv
    ├── testing_input.csv
    ├── orf_predictions.csv
    └── contig_predictions.csv
```

## Notes on this migration

To keep this repo git-friendly, several large regeneratable intermediates
were dropped from `results/` (all reproducible via the commands documented
above):

- `results/training/pairwise_aln.csv` (2.2GB) — regenerate with
  `scripts/getorfs_pairwise_aln.py` on the sampled training contigs.
- `results/curated_candidate_contigs/pairwise_aln.csv` (961MB) — same script,
  run on `../data/contigs/contigs_all_deduplicated.fasta`.
- The raw per-fold/per-ORF prediction dumps inside `evaluation_results_tuned/`
  and `eval_once_bins_5_10_15_20_t05/` (each folder's `best_method.txt`,
  `method_comparison_*.csv`, `threshold_summary.txt`, and `*_iteration_metrics.csv`
  are kept as the citable summary, matching Supplementary Figures S3/S4).
- `results/evaluation_results_fixed05/orf_predictions_results.csv` (~170MB,
  exceeds GitHub's 100MB file limit) — needed only by
  `visualizations/supp_fig6_orf-predictionss.ipynb`; regenerate via
  `train_model_pipeline.py --stage evaluate --outdir evaluation_results_fixed05
  --use_fixed_threshold --threshold 0.5`.
- `results/model_selection/` — the reproducibility-bug fix re-run completed
  on 2026-07-24 (5/5 iterations; the bug itself, missing `random_state` on
  non-RandomForest classifiers in `model_selection()`, is fixed in
  `scripts/train_model_pipeline.py`). Overall best model: RandomForest,
  0.7643 ± 0.0212 average accuracy (selected in 20/25 folds) — see
  `results/model_selection/model_selection_summary.txt`.
