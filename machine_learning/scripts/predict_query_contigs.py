#!/usr/bin/env python3
# filepath: machine_learning/scripts/predict_query_contigs.py

import os
import sys
import argparse
import pandas as pd
from joblib import load
from utils import prepare_bin_df


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Predict Tobamovirus contigs using trained models")
    parser.add_argument("input_path", help="Path to the processed test data CSV file")
    parser.add_argument("model_dir", help="Directory containing all model files")
    parser.add_argument("--outdir", help="Output directory name")
    parser.add_argument("--bin-num", type=int, default=10, help="Number of bins for binned prediction (default: 10)")
    return parser.parse_args()


def predict_contigs(input_path, model_dir, outdir, bin_num=10):
    """
    Predict Tobamovirus contigs using trained models

    Parameters:
    -----------
    input_path : str
        Path to the model input CSV file
    model_dir : str
        Directory containing all model files
    outdir : str
        Output directory name
    bin_num : int
        Number of bins for binned prediction approach

    Returns:
    --------
    pd.DataFrame
        DataFrame with contig predictions
    """
    print(f"Starting contig prediction using {bin_num} bins for binned prediction approach...")

    predictions_dir = f"{outdir}/predictions"

    # Construct model file paths
    rf_model_path = os.path.join(model_dir, "rf_model.joblib")
    rf_scaler_path = os.path.join(model_dir, "rf_scaler.joblib")
    lr_model_path = os.path.join(model_dir, f"lr_binned_{bin_num}_model.joblib")
    feature_names_path = os.path.join(model_dir, "rf_feature_names.csv")

    # Check if files exist
    for file_path in [rf_model_path, rf_scaler_path, lr_model_path, feature_names_path]:
        if not os.path.exists(file_path):
            print(f"Error: Required model file not found: {file_path}")
            # Try fallbacks for LR model if specific bin file is not found.
            # Naming evolved: histogram -> stacking -> binned.
            if "lr_binned" in file_path and not os.path.exists(file_path):
                fallback_candidates = [
                    os.path.join(model_dir, f"lr_stacking_{bin_num}_model.joblib"),
                    os.path.join(model_dir, f"lr_histogram_{bin_num}_model.joblib"),
                    os.path.join(model_dir, "lr_binned_model.joblib"),
                    os.path.join(model_dir, "lr_stacking_model.joblib"),
                    os.path.join(model_dir, "lr_histogram_model.joblib"),
                ]
                fallback_path = next((p for p in fallback_candidates if os.path.exists(p)), None)
                if fallback_path is not None:
                    print(f"Using fallback LR model: {fallback_path}")
                    lr_model_path = fallback_path
                else:
                    sys.exit(1)
            else:
                sys.exit(1)

    # Load models and data
    print("Loading models and data...")
    morf = load(rf_model_path)
    sorf = load(rf_scaler_path)
    mc = load(lr_model_path)
    input_df = pd.read_csv(input_path)

    # Create output directory
    os.makedirs(f"results/{predictions_dir}", exist_ok=True)

    # Load feature names
    cols = pd.read_csv(feature_names_path, header=None)[0].tolist()

    # Prepare input data
    print("Processing input data...")
    if input_df.index.name != "orf_name":
        input_df = input_df.set_index("orf_name")

    # Make sure all required features are present
    missing_cols = [col for col in cols if col not in input_df.columns]
    if missing_cols:
        print(f"Warning: Missing columns in input data: {missing_cols}")
        print("Proceeding with available features only.")
        cols = [col for col in cols if col in input_df.columns]

    input_df = input_df[cols]

    # Scale data
    X_test = sorf.transform(input_df)

    # Make ORF predictions
    print("Making ORF-level predictions...")
    y_pred = morf.predict(X_test)
    y_prob = morf.predict_proba(X_test)[:, 1]
    results_df = pd.DataFrame({"orf_name": input_df.index, "prediction": y_pred, "prob_1": y_prob})

    # Extract contig name
    results_df["contig_name"] = results_df["orf_name"].str.extract(r"(.*)(?:_ORF\.\d+|_aa_frame\d+)", expand=False)

    # Save ORF predictions
    results_df.to_csv(f"results/{predictions_dir}/orf_predictions.csv", index=False)
    print(f"ORF predictions saved to results/{predictions_dir}/orf_predictions.csv")

    # Prepare predictions for contig prediction
    print(f"Preparing stacking with {bin_num} bins for contig-level prediction...")
    bin_df = prepare_bin_df(results_df, "stacking", num_bins=bin_num)
    X_test = bin_df

    # Predict contigs
    print("Making contig-level predictions...")
    final_predictions = mc.predict(X_test)
    final_predictions_prob = mc.predict_proba(X_test)[:, 1]
    final_predictions_df = pd.DataFrame(
        {"contig_name": bin_df.index, "predicted_class": final_predictions, "prob_1": final_predictions_prob}
    )

    # Save contig predictions
    final_predictions_df.to_csv(f"results/{predictions_dir}/contig_predictions.csv", index=False)
    print(f"Contig predictions saved to results/{predictions_dir}/contig_predictions.csv")

    # Summary statistics
    tobamo_count = final_predictions_df["predicted_class"].sum()
    total_count = len(final_predictions_df)
    print(
        f"Summary: {tobamo_count} out of {total_count} contigs classified as Tobamovirus ({tobamo_count/total_count:.1%})"
    )

    return final_predictions_df


def main():
    """Main function to handle argument parsing and prediction workflow"""
    args = parse_args()

    # Run prediction
    predict_contigs(args.input_path, args.model_dir, args.outdir, args.bin_num)

    print("Prediction completed successfully!")


if __name__ == "__main__":
    main()
