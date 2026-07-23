import os
import subprocess
import sys
import tempfile
import warnings
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from glob import glob
from itertools import product, combinations
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import mpu
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.interpolate import interp1d
from scipy.stats import sem, t
from Bio import AlignIO, SeqIO
from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from joblib import dump, load
from scipy.optimize import curve_fit
from sklearn.calibration import calibration_curve
from sklearn.cluster import KMeans
from sklearn.discriminant_analysis import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import (
    ConfusionMatrixDisplay,
    accuracy_score,
    adjusted_rand_score,
    auc,
    classification_report,
    confusion_matrix,
    f1_score,
    make_scorer,
    mean_squared_error,
    precision_recall_curve,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import (
    GridSearchCV,
    LeaveOneOut,
    ParameterGrid,
    StratifiedKFold,
    TunedThresholdClassifierCV,
    cross_val_predict,
    train_test_split,
)
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from tqdm import tqdm
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd


def remove_duplicates(fin_path, fout_path, mapper_path):
    """Remove duplicate sequences and create mapping file"""
    with open(fin_path, "r") as fin, open(fout_path, "w") as fout, open(mapper_path, "w") as mapper:
        fin_rec = list(SeqIO.parse(fin_path, "fasta"))

        # Use a dictionary to store unique sequences and their corresponding IDs
        unique_sequences = {}
        for rec in fin_rec:
            seq_str = str(rec.seq)
            if seq_str not in unique_sequences:
                unique_sequences[seq_str] = [rec.id]
            else:
                unique_sequences[seq_str].append(rec.id)

        # Write unique sequences to the output file
        unique_records = [SeqRecord(Seq(seq), id=ids[0], description="") for seq, ids in unique_sequences.items()]
        SeqIO.write(unique_records, fout, "fasta")

        # Write the mapping of identical sequences to the mapper file
        for ids in unique_sequences.values():
            mapper.write(f"{ids[0]}: {', '.join(ids)}\n")

        print(f"    Total records: {len(fin_rec)}")
        print(f"    Unique sequences: {len(unique_sequences)}")
        print(f"    Duplicates removed: {len(fin_rec) - len(unique_sequences)}")

        return unique_sequences


def pairwise_refs(fasta_path):
    """
    Perform pairwise alignment between all sequences in a FASTA file
    Returns a DataFrame with alignment results
    """
    sequences = list(SeqIO.parse(fasta_path, "fasta"))
    results = []

    print(f"    Performing pairwise alignments for {len(sequences)} sequences...")

    if len(sequences) < 2:
        print(f"    Warning: Only {len(sequences)} sequence(s) found. No pairwise comparisons needed.")
        return pd.DataFrame(columns=["orf_name", "ref_name", "identity_score"])

    for i, seq1 in enumerate(tqdm(sequences)):
        for j, seq2 in enumerate(sequences):
            if i != j:  # Don't compare sequence with itself
                # Calculate identity score using alignment
                identity = calculate_identity_simple(str(seq1.seq), str(seq2.seq))
                results.append({"orf_name": seq1.id, "ref_name": seq2.id, "identity_score": identity})

    return pd.DataFrame(results)


def calculate_identity_simple(seq1, seq2):
    """Calculate sequence identity between two sequences using simple pairwise alignment"""
    from Bio.Align import PairwiseAligner

    aligner = PairwiseAligner()
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5

    alignment = aligner.align(seq1, seq2)[0]

    # Calculate identity as matches / alignment_length
    matches = sum(1 for a, b in zip(alignment[0], alignment[1]) if a == b and a != "-")
    alignment_length = len(alignment[0])

    return matches / alignment_length if alignment_length > 0 else 0.0


def prepend_length(ref, length_mapper):
    """Add sequence length prefix to reference ID"""
    length = length_mapper.get(ref, 0)
    return f"{length}_{ref}"


def prep_df(df, length_mapper, selected_col):
    """Prepare dataframe for distance matrix calculation"""
    if df.empty:
        return pd.DataFrame(columns=["ref1", "ref2", "diff"])

    df = df[["orf_name", "ref_name", selected_col]].copy()
    df.rename(columns={"orf_name": "ref1", "ref_name": "ref2"}, inplace=True)

    if selected_col == "identity_score":
        df["diff"] = 1 - df["identity_score"]
        df.drop(columns="identity_score", inplace=True)

    df["ref1"] = df["ref1"].apply(lambda x: prepend_length(x, length_mapper))
    df["ref2"] = df["ref2"].apply(lambda x: prepend_length(x, length_mapper))

    return df


def make_distance_matrix(df, selected_col):
    """Create distance matrix from pairwise alignment data"""
    if df.empty:
        return pd.DataFrame(), pd.DataFrame(), []

    idx = pd.concat([df["ref1"], df["ref2"]]).unique()
    col = "diff" if selected_col == "identity_score" else "clustering_score"

    res = (
        pd.pivot(df, index=["ref1"], columns=["ref2"], values=col)
        .reindex(index=idx, columns=idx)
        .fillna(0)
        .astype(float)
    )
    matrix = res + res.T

    # Ensure diagonal is 0 (distance from sequence to itself)
    np.fill_diagonal(matrix.values, 0)

    return matrix, res, idx


def make_dendrograms_genus(matrix, res, selected_col, genus):
    """Create dendrogram for genus-specific clustering with error handling"""
    try:
        condensed_distance_matrix = ssd.squareform(matrix)

        # Check for negative distances
        if np.any(condensed_distance_matrix < 0):
            print(f"      Warning: Negative distances detected for {genus}. Setting to 0.")
            condensed_distance_matrix = np.maximum(condensed_distance_matrix, 0)

        # Check for insufficient diversity
        if len(np.unique(condensed_distance_matrix)) < 3:
            print(f"      Warning: Insufficient distance diversity for {genus}.")
            return None

        linkage_matrix = sch.linkage(condensed_distance_matrix, method="single", optimal_ordering=True)
        num_labels = len(res.index)
        fig_height = max(10, num_labels * 0.3)

        plt.figure(figsize=(12, fig_height))
        dendrogram = sch.dendrogram(linkage_matrix, labels=res.index, orientation="left")
        plt.title(f"Hierarchical Dendrogram of {genus.upper()} based on {selected_col}")
        plt.xlabel("Distance")

        os.makedirs("results/reference_selection/dendrograms", exist_ok=True)
        plt.savefig(
            f"results/reference_selection/dendrograms/{genus}_RNA1_{selected_col}_dendrogram.png",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close()
        return linkage_matrix

    except Exception as e:
        print(f"      Error creating dendrogram for {genus}: {str(e)}")
        return None


def make_cluster_df_genus(linkage_matrix, idx, selected_col, genus):
    """Create cluster analysis dataframe for genus with error handling"""
    if linkage_matrix is None:
        return pd.DataFrame(columns=["distance", "identity", "clusters"])

    try:
        distances = linkage_matrix[:, 2]
        unique_distances = sorted(np.unique(distances))
        results = []

        for dist in unique_distances:
            clusters = sch.fcluster(linkage_matrix, t=dist, criterion="distance")
            cluster_labels = {i: [] for i in np.unique(clusters)}

            for label, cluster in zip(idx, clusters):
                cluster_labels[cluster].append(label)

            results.append({"distance": dist, "identity": 1 - dist, "clusters": len(cluster_labels)})

        clusters_df = pd.DataFrame(results)
        os.makedirs("results/reference_selection/clusters", exist_ok=True)
        clusters_df.to_csv(
            f"results/reference_selection/clusters/{genus}_RNA1_{selected_col}_clusters.csv", index=False
        )
        return clusters_df

    except Exception as e:
        print(f"      Error creating cluster analysis for {genus}: {str(e)}")
        return pd.DataFrame(columns=["distance", "identity", "clusters"])


# Function to select one species per cluster
def select_one_species(group):
    if group["orig_ref"].sum() > 0:
        return group[group["orig_ref"] == 1].iloc[0]
    else:
        return group.loc[group["length"].idxmax()]


# select species
def select_species(linkage_matrix, distance_threshold, idx, orig_ref, selected_col, genus):
    clusters = sch.fcluster(linkage_matrix, t=distance_threshold, criterion="distance")
    cluster_labels = {i: [] for i in np.unique(clusters)}

    for label, cluster in zip(idx, clusters):
        cluster_labels[cluster].append(label)

    cluster_labels
    data = []
    for cluster, species_list in cluster_labels.items():
        for species in species_list:
            data.append({"cluster": cluster, "species": species})

    # Convert the list of dictionaries to a DataFrame
    dd = pd.DataFrame(data)
    dd["length"] = dd["species"].apply(lambda x: int(x.split("_")[0]))
    dd["orig_ref"] = np.where(dd["species"].isin(orig_ref), 1, 0)

    # Apply the function to each cluster
    selected_species = dd.groupby("cluster").apply(select_one_species, include_groups=False).reset_index(drop=True)
    selected_species_list = selected_species["species"].tolist()

    os.makedirs("results/selected_species", exist_ok=True)
    selected_species.to_csv(f"results/selected_species/{genus}_{selected_col}_RNA1_selected_species.csv", index=False)
    return selected_species_list


def make_distance_matrix(df, selected_col):
    idx = pd.concat([df["ref1"], df["ref2"]]).unique()
    if selected_col == "identity_score":
        col = "diff"
    else:
        col = "clustering_score"
    res = (
        pd.pivot(df, index=["ref1"], columns=["ref2"], values=col)
        .reindex(index=idx, columns=idx)
        .fillna(0)
        .astype(float)
    )
    matrix = res + res.T
    return matrix, res, idx


def make_dendrograms(matrix, res, selected_col, genus):
    condensed_distance_matrix = ssd.squareform(matrix)
    linkage_matrix = sch.linkage(condensed_distance_matrix, method="single", optimal_ordering=True)
    num_labels = len(res.index)
    fig_height = max(10, num_labels * 0.3)  # Adjust the multiplier as needed for better fit

    # Create the dendrogram plot
    plt.figure(figsize=(12, fig_height))
    dendrogram = sch.dendrogram(linkage_matrix, labels=res.index, orientation="left")
    plt.title(f"Hierarchical Dendrogram of {genus} RNA1 base on {selected_col}")
    plt.xlabel("Distance")
    os.makedirs("results/dendrograms", exist_ok=True)
    # plt.savefig(f'results/dendrograms/{genus}_RNA1_{selected_col}_dendrogram.png', dpi=300, bbox_inches='tight')
    plt.close()
    return linkage_matrix


def make_dendrogram_selected_species(matrix, res, selected_col, genus, selected_species_list, distance_threshold):
    # Convert the distance matrix to condensed form
    condensed_distance_matrix = ssd.squareform(matrix)

    # Perform hierarchical clustering
    linkage_matrix = sch.linkage(condensed_distance_matrix, method="single", optimal_ordering=True)

    # Determine the figure size based on the number of labels
    num_labels = len(res.index)
    fig_height = max(10, num_labels * 0.3)  # Adjust the multiplier as needed for better fit

    # Create the dendrogram plot
    plt.figure(figsize=(12, fig_height))
    dendrogram = sch.dendrogram(
        linkage_matrix,
        labels=res.index,
        orientation="left",
        color_threshold=distance_threshold,
        above_threshold_color="grey",
    )
    plt.title(
        f"Hierarchical Dendrogram of {genus} RNA1 based on {selected_col}: identity threshold = {1-distance_threshold:.2f}"
    )
    plt.xlabel("Distance")

    # Get the current axis
    ax = plt.gca()

    # Get the y-axis labels
    y_labels = ax.get_ymajorticklabels()

    # Make labels from selected_species_list in bold text
    for label in y_labels:
        if label.get_text() in selected_species_list:
            label.set_fontweight("bold")
        # if label.get_text() in all_ids_of_original_seqs_:
        #     label.set_color('darkgreen')

    # Plot a vertical line at the specified distance
    plt.axvline(x=distance_threshold, color="black", linestyle="--")

    # Adjust layout to ensure labels are not cut off
    plt.subplots_adjust(left=0.5, right=0.9, top=0.9, bottom=0.1)  # Adjust margins as needed

    # Create the output directory if it doesn't exist
    out_path = "results/reference_selection/dendrograms_selected_species"
    os.makedirs("out_path", exist_ok=True)

    # Save the dendrogram plot
    plt.savefig(f"{out_path}/{genus}_RNA1_{selected_col}_dendrogram.png", dpi=100, bbox_inches="tight")
    # plt.show()
    plt.close()


def make_cluster_df(linkage_matrix, idx, selected_col, genus):
    distances = linkage_matrix[:, 2]
    unique_distances = sorted(np.unique(distances))
    results = []
    for dist in unique_distances:
        clusters = sch.fcluster(linkage_matrix, t=dist, criterion="distance")
        cluster_labels = {i: [] for i in np.unique(clusters)}

        for label, cluster in zip(idx, clusters):
            cluster_labels[cluster].append(label)

        results.append({"distance": dist, "identity": 1 - dist, "clusters": len(cluster_labels)})

    clusters_df = pd.DataFrame(results)
    os.makedirs("results/clusters", exist_ok=True)
    clusters_df.to_csv(f"results/clusters/{genus}_RNA1_{selected_col}_clusters.csv", index=False)
    return clusters_df


def train_lr_final_model(all_predictions, bin_num, fixed_threshold=0.5):
    """Train final logistic regression model with fixed threshold"""
    bin_df = prepare_bin_df_refs_hist_test(all_predictions, "stacking", num_bins=bin_num)
    X = bin_df.drop(columns=["ground_truth"])
    y = bin_df["ground_truth"]

    # Use regular LogisticRegression for final model
    model_lr = LogisticRegression(max_iter=1000)
    model_lr.fit(X, y)

    return {"model": model_lr, "threshold": fixed_threshold, "bin_df": bin_df}  # Return bin_df for saving


def check_file_exists(filepath, filetype):
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"{filetype} file not found: {filepath}")
    if filetype == "Excel" and not filepath.endswith((".xls", ".xlsx")):
        raise ValueError(f"{filetype} file is not in the correct format: {filepath}")
    if filetype == "CSV" and not filepath.endswith(".csv"):
        raise ValueError(f"{filetype} file is not in the correct format: {filepath}")
    if filetype == "FASTA" and not filepath.endswith((".fasta", ".fa", ".fna")):
        raise ValueError(f"{filetype} file is not in the correct format: {filepath}")


def make_mappers(refs):
    """
    Create mappers from the reference DataFrame.

    Parameters:
    refs (pd.DataFrame): DataFrame containing reference information.

    Returns:
    tuple: A tuple containing four dictionaries:
        - aa_id2type: Maps record IDs to their types.
        - aa_id2prot: Maps record IDs to their protein types.
        - inv_type_acc_dict: Maps accession IDs to their types.
        - inv_nt_virusname2id: Maps accession IDs to their virus names.
    """
    aa_id2type = {}
    aa_id2prot = {}
    inv_type_acc_dict = {}
    inv_nt_virusname2id = {}

    for idx, row in refs.iterrows():
        type_value = row["type"]
        for col in ["orf1_record_id", "orf2_record_id", "cp_record_id"]:
            record_id = row[col]
            if pd.notna(record_id):
                aa_id2type[record_id] = type_value

    for idx, row in refs.iterrows():
        for col, prot_type in [("orf1_record_id", "orf1"), ("orf2_record_id", "orf2"), ("cp_record_id", "cp")]:
            record_id = row[col]
            if pd.notna(record_id):
                aa_id2prot[record_id] = prot_type

    for _, row in refs.iterrows():
        inv_type_acc_dict[row["accession"]] = row["type"]
        inv_nt_virusname2id[row["accession"]] = row["virus_name"]

    return aa_id2type, aa_id2prot, inv_type_acc_dict, inv_nt_virusname2id


def filter_extreme_probability(df, idx, refs_: bool):
    """
    Filter the DataFrame to keep only the row with the most extreme probability for each contig_name.

    Parameters:
    df (pd.DataFrame): The input DataFrame containing probabilities.
    idx (int): The fold index.
    refs_ (bool): Whether to include ground truth in the output.

    Returns:
    pd.DataFrame: The filtered DataFrame with the most extreme probabilities.
    """
    # Calculate distance from 0.5 (the midpoint) to find the most extreme probability
    df["extreme_score"] = df["prob_1"].apply(lambda p: abs(p - 0.5))

    # Sort by contig_name and extreme_score (to select the most extreme probability)
    df_sorted = df.sort_values(by=["contig_name", "extreme_score"], ascending=[True, False])

    # Keep only the row with the maximum extreme score for each contig_name
    filtered_df = df_sorted.groupby("contig_name").head(1)

    # Drop the 'extreme_score' column as it's no longer needed
    filtered_df = filtered_df.drop(columns=["extreme_score"])

    # Add idx and mc_name columns
    filtered_df["fold_idx"] = idx
    filtered_df["mc_name"] = "most_extreme"

    # Select columns to include in the output
    if refs_:
        cols = ["fold_idx", "mc_name", "contig_name", "prob_1", "ground_truth"]
    else:
        cols = ["fold_idx", "mc_name", "contig_name", "prob_1"]

    out = filtered_df[cols]

    return out


def predict_orfs(test, morf, sorf, refs: bool = True):
    """
    Predict ORFs using a trained model and scaler.

    Parameters:
    test (pd.DataFrame): The test DataFrame containing ORF information.
    morf (RandomForestClassifier): The trained Random Forest model.
    sorf (StandardScaler): The scaler used to standardize the data.
    refs (bool): Whether to include ground truth in the predictions. Default is True.

    Returns:
    pd.DataFrame: A DataFrame containing predictions and ground truth (if refs is True).
    """
    # Check if the columns to be dropped are present in the DataFrame
    columns_to_drop = ["orf_type", "strand", "virus_name", "accession", "contig_name"]
    columns_present = [col for col in columns_to_drop if col in test.columns]

    # Check if the index is set to 'orf_name', if not set it
    if test.index.name != "orf_name":
        test = test.set_index("orf_name")

    # Drop unnecessary columns
    X_test = test.drop(columns=columns_present)

    # Save the index of the test DataFrame
    X_test_indices = X_test.index

    # Transform the test data
    X_test = sorf.transform(X_test)

    # Make predictions
    predictions = morf.predict_proba(X_test)[:, 1]

    # Create a DataFrame with predictions
    if refs:
        y_test = (test["orf_type"] == "tobamo").astype(int)
        predictions_df = pd.DataFrame(
            {"orf_name": X_test_indices, "prob_1": predictions, "ground_truth": y_test.values}
        )
    else:
        predictions_df = pd.DataFrame({"orf_name": X_test_indices, "prob_1": predictions})

    # Extract contig_name from orf_name
    predictions_df["contig_name"] = predictions_df["orf_name"].str.extract(r"(.*)(?:_ORF\.\d+|_aa_frame\d+)")

    return predictions_df


def train_lr_and_predict(all_predictions, methods=["stacking", "cumsum", "cumsum_reverse"]):
    """
    Train Logistic Regression models for the specified methods.

    Parameters:
    all_predictions (pd.DataFrame): DataFrame containing predictions and ground truth.
    methods (list or str): List of methods or a single method to use for training.
                           Defaults to ["stacking", "cumsum", "cumsum_reverse"].

    Returns:
    dict: A dictionary of trained Logistic Regression models for each method.
    """
    # Ensure methods is a list, even if a single method is passed as a string
    if isinstance(methods, str):
        methods = [methods]

    models = {}

    for method_name in methods:
        bin_df = prepare_bin_df_refs(all_predictions, method_name, num_bins=10)
        X = bin_df.drop(columns=["ground_truth"])
        y = bin_df["ground_truth"]
        model_lr = LogisticRegression(max_iter=1000)
        model_lr.fit(X, y)

        models[method_name] = model_lr

    return models


def train_lr_and_predict_binned_test(all_predictions, nums: list):
    models = {}

    for num in nums:
        bin_df = prepare_bin_df_refs_hist_test(all_predictions, "stacking", num_bins=num)
        X = bin_df.drop(columns=["ground_truth"])
        y = bin_df["ground_truth"]

        scorer = make_scorer(accuracy_score)
        base_model = LogisticRegression(max_iter=1000)
        model_lr = TunedThresholdClassifierCV(base_model, scoring=scorer, cv=StratifiedKFold(n_splits=10))

        model_lr.fit(X, y)

        # Store both the model and the best threshold
        models[f"stacking_{num}"] = {"model": model_lr, "best_threshold": model_lr.best_threshold_}

    return models


# def train_lr_final_model(all_predictions, num_bins, fixed_threshold=0.5):
#     """Train final logistic regression model with fixed threshold"""
#     bin_df = prepare_bin_df_refs_hist_test(all_predictions, "stacking", num_bins=num_bins)
#     X = bin_df.drop(columns=["ground_truth"])
#     y = bin_df["ground_truth"]

#     # Use regular LogisticRegression for final model
#     model_lr = LogisticRegression(max_iter=1000)
#     model_lr.fit(X, y)

#     return {
#         'model': model_lr,
#         'threshold': fixed_threshold
#     }


def prepare_bin_df_refs_hist_test(all_predictions, method_name, num_bins):
    """
    Prepare the bin_df DataFrame from all_predictions.

    Parameters:
    all_predictions (pd.DataFrame): DataFrame containing predictions and ground truth.
    num_bins (int): Number of bins for the stacking.
    method (str): Method for stacking calculation ('stacking', 'cumsum', 'cumsum_reverse').

    Returns:
    pd.DataFrame: The bin_df DataFrame with binned probabilities and ground truth.
    """
    # Raise error if method is not valid
    if method_name not in ["stacking", "cumsum", "cumsum_reverse"]:
        raise ValueError(
            f"Unknown method: {method_name}. Please choose from 'stacking', 'cumsum', or 'cumsum_reverse'."
        )

    # Select the appropriate stacking calculation function
    if method_name == "stacking":
        histogram_func = calculate_histogram
    elif method_name == "cumsum":
        histogram_func = calculate_histogram_cumsum
    elif method_name == "cumsum_reverse":
        histogram_func = calculate_histogram_cumsum_reverse

    # Create a ground truth mapper
    ground_truth_mapper = dict(all_predictions.groupby("contig_name")["ground_truth"].first())

    # Calculate histogram for each contig_name
    bin_df = all_predictions.groupby("contig_name").apply(histogram_func, num_bins=num_bins, include_groups=False)

    # # Reset index and map ground truth
    bin_df = bin_df.reset_index()
    bin_df["ground_truth"] = bin_df["contig_name"].map(ground_truth_mapper)
    bin_df.set_index("contig_name", inplace=True)

    return bin_df


def prepare_bin_df_refs(all_predictions, method_name, num_bins=10):
    """
    Prepare the bin_df DataFrame from all_predictions.

    Parameters:
    all_predictions (pd.DataFrame): DataFrame containing predictions and ground truth.
    num_bins (int): Number of bins for the stacking.
    method (str): Method for stacking calculation ('stacking', 'cumsum', 'cumsum_reverse').

    Returns:
    pd.DataFrame: The bin_df DataFrame with binned probabilities and ground truth.
    """
    # Raise error if method is not valid
    if method_name not in ["stacking", "cumsum", "cumsum_reverse"]:
        raise ValueError(
            f"Unknown method: {method_name}. Please choose from 'stacking', 'cumsum', or 'cumsum_reverse'."
        )

    # Select the appropriate stacking calculation function
    if method_name == "stacking":
        histogram_func = calculate_histogram
    elif method_name == "cumsum":
        histogram_func = calculate_histogram_cumsum
    elif method_name == "cumsum_reverse":
        histogram_func = calculate_histogram_cumsum_reverse

    # Create a ground truth mapper
    ground_truth_mapper = dict(all_predictions.groupby("contig_name")["ground_truth"].first())

    # Calculate stacking for each contig_name
    bin_df = all_predictions.groupby("contig_name").apply(histogram_func, num_bins=num_bins, include_groups=False)

    # # Reset index and map ground truth
    bin_df = bin_df.reset_index()
    bin_df["ground_truth"] = bin_df["contig_name"].map(ground_truth_mapper)
    bin_df.set_index("contig_name", inplace=True)

    return bin_df


# Function to prepare bin_df
def prepare_bin_df(all_predictions, method_name, num_bins=10):
    """
    Prepare the bin_df DataFrame from all_predictions.

    Parameters:
    all_predictions (pd.DataFrame): DataFrame containing predictions and ground truth.
    num_bins (int): Number of bins for the stacking.
    method (str): Method for stacking calculation ('stacking', 'cumsum', 'cumsum_reverse').

    Returns:
    pd.DataFrame: The bin_df DataFrame with binned probabilities and ground truth.
    """
    # Raise error if method is not valid
    if method_name not in ["stacking", "cumsum", "cumsum_reverse"]:
        raise ValueError(
            f"Unknown method: {method_name}. Please choose from 'stacking', 'cumsum', or 'cumsum_reverse'."
        )

    # Select the appropriate stacking calculation function
    if method_name == "stacking":
        histogram_func = calculate_histogram
    elif method_name == "cumsum":
        histogram_func = calculate_histogram_cumsum
    elif method_name == "cumsum_reverse":
        histogram_func = calculate_histogram_cumsum_reverse

    # Calculate stacking for each contig_name
    bin_df = all_predictions.groupby("contig_name").apply(histogram_func, num_bins=num_bins, include_groups=False)

    return bin_df


def generate_threshold_summary(binned_results, outdir):
    """Generate a concise threshold analysis summary"""

    # Calculate basic statistics
    overall_mean = binned_results["threshold"].mean()
    overall_std = binned_results["threshold"].std()
    overall_min = binned_results["threshold"].min()
    overall_max = binned_results["threshold"].max()

    # Per-iteration means
    iter_means = binned_results.groupby("iteration")["threshold"].mean()
    iter_mean_avg = iter_means.mean()
    iter_mean_std = iter_means.std()

    # Simple recommendation
    if abs(overall_mean - 0.5) < 0.05:
        recommendation = "Use default threshold 0.5"
    else:
        recommendation = f"Consider CV threshold: {overall_mean:.3f}"

    # Create short report
    report = f"""
        THRESHOLD ANALYSIS SUMMARY
        ==========================
        Overall statistics:
        Mean: {overall_mean:.4f} ± {overall_std:.4f}
        Range: {overall_min:.4f} - {overall_max:.4f}

        Per-iteration means:
        Mean: {iter_mean_avg:.4f} ± {iter_mean_std:.4f}
        
        Recommendation: {recommendation}
        """

    # Save to file
    with open(f"results/{outdir}/threshold_summary.txt", "w") as f:
        f.write(report)

    # Print to console
    print(report)

    return overall_mean, overall_std


def predict_contigs(predictions, mc, idx, mc_name, refs_: bool):

    if refs_:
        bin_df = prepare_bin_df_refs(predictions, mc_name, num_bins=10)

        X_test = bin_df.drop(columns=["ground_truth"])
        final_predictions = mc.predict_proba(X_test)[:, 1]

        y_test = bin_df["ground_truth"]
        final_predictions_df = pd.DataFrame(
            {
                "fold_idx": idx,
                "mc_name": mc_name,
                "contig_name": bin_df.index,
                "prob_1": final_predictions,
                "ground_truth": y_test,
            }
        )

    else:
        bin_df = prepare_bin_df(predictions, mc_name, num_bins=10)
        X_test = bin_df
        final_predictions = mc.predict_proba(X_test)[:, 1]
        final_predictions_df = pd.DataFrame(
            {"fold_idx": idx, "mc_name": mc_name, "contig_name": bin_df.index, "prob_1": final_predictions}
        )

    return final_predictions_df


def predict_contigs_binned_test(predictions, mc, idx, mc_name, num: int, refs_: bool):

    if refs_:
        bin_df = prepare_bin_df_refs_hist_test(predictions, mc_name, num_bins=num)

        X_test = bin_df.drop(columns=["ground_truth"])
        final_predictions = mc.predict_proba(X_test)[:, 1]

        y_test = bin_df["ground_truth"]
        final_predictions_df = pd.DataFrame(
            {
                "fold_idx": idx,
                "mc_name": mc_name,
                "contig_name": bin_df.index,
                "prob_1": final_predictions,
                "ground_truth": y_test,
            }
        )

    else:
        bin_df = prepare_bin_df(predictions, mc_name, num_bins=num)
        X_test = bin_df
        final_predictions = mc.predict_proba(X_test)[:, 1]
        final_predictions_df = pd.DataFrame(
            {"fold_idx": idx, "mc_name": mc_name, "contig_name": bin_df.index, "prob_1": final_predictions}
        )

    return final_predictions_df


def train_rf_and_predict(train, selected_model):
    """
    Train a Random Forest model and predict probabilities for each virus name.

    Parameters:
    train (pd.DataFrame): The training DataFrame.
    selected_model (dict): The selected model configuration.

    Returns:
    pd.DataFrame: A DataFrame containing predictions and ground truth.
    """

    unique_virusnames = train["virus_name"].unique()
    all_predictions = pd.DataFrame()

    # Perform LOOCV for each unique virus name to make ORF predictions
    for virusname in unique_virusnames:
        # Split data
        sub_train = train[train["virus_name"] != virusname]
        sub_test = train[train["virus_name"] == virusname]

        X_train = sub_train.drop(columns=["orf_type", "strand", "virus_name", "accession", "contig_name"])
        y_train = (sub_train["orf_type"] == "tobamo").astype(int)
        X_test = sub_test.drop(columns=["orf_type", "strand", "virus_name", "accession", "contig_name"])
        y_test = (sub_test["orf_type"] == "tobamo").astype(int)

        # Save original indices of X_test
        X_test_indices = X_test.index

        # Standardize training and test data
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)

        # Initialize the model
        model = selected_model["model"]

        # Train the model
        model.fit(X_train, y_train)

        # Predict and save probabilities of the test set
        predictions = model.predict_proba(X_test)[:, 1]
        predictions_df = pd.DataFrame(
            {"orf_name": X_test_indices, "prob_1": predictions, "ground_truth": y_test.values}
        )

        predictions_df["contig_name"] = predictions_df["orf_name"].str.extract(r"(.*)(?:_ORF\.\d+|_aa_frame\d+)")
        # Append the predictions to the all_predictions dataframe
        all_predictions = pd.concat([all_predictions, predictions_df], ignore_index=True)

    return all_predictions


def train_rf_on_all_data(train, selected_model):
    """
    Train a Random Forest model on all data and save the model and scaler.

    Parameters:
    train (pd.DataFrame): The training DataFrame.
    model_path (str): Path to save the trained model.
    scaler_path (str): Path to save the scaler.

    Returns:
    tuple: A tuple containing the trained model and scaler.
    """
    X_train = train.drop(columns=["orf_type", "strand", "virus_name", "accession", "contig_name"])
    y_train = (train["orf_type"] == "tobamo").astype(int)
    feature_names = X_train.columns.tolist()

    # Standardize training data
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)

    # Initialize and train the model
    model = selected_model["model"]
    model.fit(X_train, y_train)

    return model, scaler, feature_names


def add_info_to_training_input_df(input_df, orf_fasta_path: str, inv_nt_virusname2id, inv_type_acc_dict):
    """
    Add additional information to the training input DataFrame from various sources.

    Parameters:
    input_df (pd.DataFrame): The input DataFrame.
    orf_fasta_path (str): Path to the ORF FASTA file.
    inv_nt_virusname2id_path (str): Path to the file containing nucleotide to virus name mapping.
    inv_type_acc_dict_path (str): Path to the file containing type to accession mapping.

    Returns:
    pd.DataFrame: The input DataFrame with additional information.
    """
    # Extract accession from orf_name
    input_df["accession"] = input_df["orf_name"].str.extract(r"^(.*\.\d)\_")

    # add species info
    input_df["virus_name"] = input_df["accession"].map({k: v for k, v in inv_nt_virusname2id.items()})

    # Add orf_type info
    input_df["orf_type"] = input_df["accession"].map(inv_type_acc_dict)

    # Add strand orientation info
    with open(orf_fasta_path, "r") as file:
        orf_dict = {
            record.description.split()[0]: "REVERSE" if "(-)" in record.description else "FORWARD"
            for record in SeqIO.parse(file, "fasta")
        }
    input_df["strand"] = input_df["orf_name"].map(orf_dict)

    return input_df


def add_info_basic(
    input_df, orf_fasta_path: str, snakemake: bool, contig_fasta_path=None
):  # , contig_fasta_path: str):
    """
    Add additional information to the input DataFrame from ORF and contig FASTA files.

    Parameters:
    input_df (pd.DataFrame): The input DataFrame.
    orf_fasta_path (str): Path to the ORF FASTA file.
    contig_fasta_path (str): Path to the contig FASTA file.

    Returns:
    pd.DataFrame: The input DataFrame with additional information.
    """

    if snakemake:
        # replace "=" with "_" and extract contig_name
        input_df["orf_name"] = input_df["orf_name"].str.replace("=", "_")
        input_df["contig_name"] = input_df["orf_name"].str.extract(r"(.*)(?:_ORF\.\d+|_aa_frame\d+)", expand=False)
        # add orf_len
        orf_lens = {
            k.replace("=", "_"): len(v.seq) for k, v in SeqIO.to_dict(SeqIO.parse(orf_fasta_path, "fasta")).items()
        }
        input_df["orf_len"] = input_df["orf_name"].map(orf_lens)
        # add contig_length
        input_df["contig_length"] = input_df["contig_name"].str.extract(r"(?:len_|length_)(\d+)").astype(int)
        # add stop_codons
        all_orfs_stops = {
            k.replace("=", "_"): int(v.seq.count("*"))
            for k, v in SeqIO.to_dict(SeqIO.parse(orf_fasta_path, "fasta")).items()
        }
        input_df["stop_codons"] = input_df["orf_name"].map(all_orfs_stops)

    else:
        # add contig_name
        input_df["contig_name"] = input_df["orf_name"].str.extract(r"(.*)(?:_ORF\.\d+|_aa_frame\d+)")
        # add orf_len
        orf_lens = {k: len(v.seq) for k, v in SeqIO.to_dict(SeqIO.parse(orf_fasta_path, "fasta")).items()}
        input_df["orf_len"] = input_df["orf_name"].map(orf_lens)
        # add contig_length
        contig_lens = {rec.id: len(rec.seq) for rec in SeqIO.parse(contig_fasta_path, "fasta")}
        input_df["contig_length"] = input_df["contig_name"].map(contig_lens)
        # add stop_codons
        all_orfs_stops = {
            k: int(v.seq.count("*")) for k, v in SeqIO.to_dict(SeqIO.parse(orf_fasta_path, "fasta")).items()
        }
        input_df["stop_codons"] = input_df["orf_name"].map(all_orfs_stops)

    # add orf_acc_method info
    input_df["orf_acc_method"] = np.where(input_df["orf_name"].str.contains("frame"), "translate", "orfipy")
    orf_acc_mapping = {"orfipy": 0, "translate": 1}
    input_df["orf_acc_method"] = input_df["orf_acc_method"].map(orf_acc_mapping)

    # Check if all dictionary keys are in the 'orf_name' column of input_df
    missing_keys = set(all_orfs_stops.keys()) - set(input_df["orf_name"])

    if missing_keys:
        print(f"Error: The following ORF names are missing in the 'orf_name' column of input_df: {missing_keys[0:10]}")
        sys.exit(1)

    return input_df


def aggregate_df(df, aa_id2type, aa_id2prot):
    df = df.copy()
    # Apply mappings to create new columns
    df["ref_type"] = df["ref_name"].map(aa_id2type)
    df["ref_prot_type"] = df["ref_name"].map(aa_id2prot)

    # Group by 'orf_name', 'ref_type', and 'ref_prot_type' and aggregate min and max values
    grouped = df.groupby(["orf_name", "ref_type", "ref_prot_type"]).agg(
        {
            "identity_score": ["min", "max", "mean"],
            "gap_openings": ["min", "max", "mean"],
            "gap_ratio": ["min", "max", "mean"],
            "N/aln_len": ["min", "max", "mean"],
            "aln_orf_len": ["min", "max", "mean"],
            "M": ["min", "max", "mean"],
            "N": ["min", "max", "mean"],
            "aln_len": ["min", "max", "mean"],
        }
    )

    # Flatten the multi-index column names
    grouped.columns = ["{}_{}".format(col[0], col[1]) for col in grouped.columns]

    # Reset index to make 'orf_name' a column again
    new_grouped = grouped.reset_index()

    assert not new_grouped.isnull().values.any(), "There are NaN values in the dataframe"

    # return new_grouped
    return new_grouped


def pivot_df(
    df,
    orf_name_col="orf_name",
    ref_type_col="ref_type",
    ref_prot_type_col="ref_prot_type",
):
    ref_types = df[ref_type_col].unique()
    ref_prot_types = df[ref_prot_type_col].unique()
    groups = list(product(ref_types, ref_prot_types))

    df_groups = {}
    for g in groups:
        temp_df = df.copy()
        temp_df = temp_df[(temp_df[ref_type_col] == g[0]) & (temp_df[ref_prot_type_col] == g[1])]
        temp_df = temp_df.set_index(orf_name_col)
        temp_df = temp_df.drop(columns=[ref_type_col, ref_prot_type_col])
        temp_df.columns = [f"{g[0]}_{g[1]}_{c}" for c in temp_df.columns]
        df_groups[g] = temp_df.copy()

    df_pivoted = df_groups[groups[0]]
    for g in groups[1:]:
        df_pivoted = df_pivoted.join(df_groups[g], how="outer")

    df_pivoted = df_pivoted.reset_index()

    return df_pivoted


def prepare_train_test(input_df, train_refs, test_refs):
    """
    Prepare training and testing datasets from the input DataFrame and reference DataFrames.

    Parameters:
    input_df (pd.DataFrame): The input DataFrame containing ORF information.
    train_refs (pd.DataFrame): The DataFrame containing training reference information.
    test_refs (pd.DataFrame): The DataFrame containing testing reference information.

    Returns:
    tuple: A tuple containing the training and testing DataFrames.
    """
    # Set 'orf_name' as the index if it's not already set
    if input_df.index.name != "orf_name":
        input_df = input_df.set_index("orf_name")

    # Get lists of nucleotide IDs for training and testing
    train_nt_list = train_refs["accession"].tolist()
    test_nt_list = test_refs["accession"].tolist()

    # Keep only selected contigs for training and testing
    train = input_df[input_df["accession"].isin(train_nt_list)]
    test = input_df[input_df["accession"].isin(test_nt_list)]

    return train, test


def stratified_kfold_split(refs, n_splits: int, random_state: int):
    """
    Perform a stratified k-fold split on the refs DataFrame based on unique virus names and their types.

    Parameters:
    refs (pd.DataFrame): The input DataFrame containing 'virus_name' and 'type' columns.
    n_splits (int): The number of splits for k-fold cross-validation.
    random_state (int): The random seed for reproducibility.

    Returns:
    list: A list of tuples containing training and testing DataFrames for each fold.
    """
    # Extract unique virus names and their corresponding types
    unique_virus_names = refs[["virus_name", "type"]].drop_duplicates()

    # Initialize StratifiedKFold with the specified number of splits
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)

    # Perform the splits on unique virus names
    folds = []
    for train_index, test_index in skf.split(unique_virus_names["virus_name"], unique_virus_names["type"]):
        train_virus_names = unique_virus_names.iloc[train_index]["virus_name"]
        test_virus_names = unique_virus_names.iloc[test_index]["virus_name"]

        # Map the split back to the original DataFrame
        train_refs = refs[refs["virus_name"].isin(train_virus_names)].copy()
        test_refs = refs[refs["virus_name"].isin(test_virus_names)].copy()

        # Calculate sampling_prob for each test and train reference group
        train_refs = train_refs.loc[train_refs["training"] == 1]  # keep only training references
        train_refs.loc[:, "sampling_prob"] = train_refs["length"] / train_refs.groupby("virus_name")[
            "length"
        ].transform("sum")
        test_refs.loc[:, "sampling_prob"] = test_refs["length"] / test_refs.groupby("virus_name")["length"].transform(
            "sum"
        )

        folds.append((train_refs, test_refs))

    return folds


def make_mapper(refs):
    """
    Create a mapper from nucleotide accession to associated amino acid record IDs.

    Args:
        refs (pd.DataFrame): DataFrame containing reference information.

    Returns:
        dict: A dictionary mapping nucleotide accession to a list of associated amino acid record IDs.
    """
    nt_aa_mapper = {}
    for idx, row in refs.iterrows():
        record_id = row["accession"]
        orf1_record_id = row["orf1_record_id"]
        orf2_record_id = row["orf2_record_id"]
        cp_record_id = row["cp_record_id"]
        values = [
            val for val in (orf1_record_id, orf2_record_id, cp_record_id) if pd.notna(val)
        ]  # Filter out NaN values
        nt_aa_mapper[record_id] = values
    return nt_aa_mapper


def remove_associated_aa_refs(dataframe, nt_aa_mapper):
    """
    Filter out rows in the DataFrame where the accession matches and ref_name is in the list of associated amino acid record IDs.

    Args:
        dataframe (pd.DataFrame): DataFrame to be filtered.
        nt_aa_mapper (dict): Dictionary mapping nucleotide accession to associated amino acid record IDs.

    Returns:
        tuple: A tuple containing the filtered DataFrame and the removed references DataFrame.
    """
    # Check if the required columns exist in the DataFrame
    if "accession" not in dataframe.columns or "ref_name" not in dataframe.columns:
        raise KeyError("The DataFrame must contain 'accession' and 'ref_name' columns.")

    # Create a boolean mask for filtering
    mask = pd.Series([True] * len(dataframe))

    # Iterate over the nt_aa_mapper dictionary
    for accession, ref_names in nt_aa_mapper.items():
        # Update the mask to filter out rows where accession matches and ref_name is in the list of ref_names
        mask &= ~((dataframe["accession"] == accession) & (dataframe["ref_name"].isin(ref_names)))

    # Apply the mask to the DataFrame
    filtered_df = dataframe[mask]
    removed_refs = dataframe[~mask]

    return filtered_df, removed_refs


def filter_pairwise(pairwise, refs):
    """
    Adds contig_name and accession to dataframe and filters out parent amino acid sequences.

    Args:
        pairwise (pd.DataFrame): DataFrame containing pairwise information.
        refs (pd.DataFrame): DataFrame containing reference information.

    Returns:
        tuple: A tuple containing the filtered pairwise DataFrame and the removed references DataFrame.
    """
    pairwise = pairwise.copy()
    # Parse the contig_name and accession from orf_name
    pairwise["contig_name"] = pairwise["orf_name"].str.extract(r"(.*)(?:_ORF\.\d+|_aa_frame\d+)")
    pairwise["accession"] = pairwise["orf_name"].str.extract(r"^(.*\.\d)\_")

    # Filter out the non-aa sequences
    nt_aa_mapper = make_mapper(refs)
    pairwise_minus, removed_refs = remove_associated_aa_refs(pairwise, nt_aa_mapper)

    return pairwise_minus, removed_refs


# make pairs and run pairwise
def pairwise_refs(ref_path):
    ref_seqs = list(SeqIO.parse(ref_path, "fasta"))
    args_list = []
    for seq1, seq2 in combinations(ref_seqs, 2):
        args_list.append([seq1, seq2])
    results = process_in_chunks(args_list, 1000, 32)
    df = pd.DataFrame(results)
    return df


def subsample_contigs(refs: pd.DataFrame, contigs, num: int, output_dir: str, random_seed):
    """
    Subsample contigs for each virus name and save the results as JSON files.

    Args:
        pairwise_minus (pd.DataFrame): DataFrame containing filtered pairwise information.
        refs (pd.DataFrame): DataFrame containing reference information.
        contigs (dict): Dictionary containing contig information.
        num (int): Number of contigs to sample for each virus name. Default is 30.
        output_dir (str): Directory to save the output JSON files.
        random_seed (int): Random seed for reproducibility.

    Returns:
        tuple: A tuple containing the filtered contigs dictionary and the sampling counts dictionary.
    """

    # Set the random seed for reproducibility
    np.random.seed(random_seed)

    all_filtered_contigs = {}
    all_sampling_counts = {}

    for virus_name, accessions in (
        refs.groupby("virus_name")[["accession", "sampling_prob"]]
        .apply(lambda x: list(x.itertuples(index=False, name=None)))
        .items()
    ):
        species_contigs = [contig.id for contig in contigs.values() if any(acc in contig.id for acc, _ in accessions)]
        filtered_contigs_ids = []
        sampling_counts = defaultdict(int)

        # Ensure one untrimmed contig (_NT) for each accession
        untrimmed_contig_ids = []
        for accession, _ in accessions:
            untrimmed_contigs_for_accession = [
                contig_id for contig_id in species_contigs if "_NT" in contig_id and accession in contig_id
            ]
            untrimmed_contig_ids.extend(untrimmed_contigs_for_accession)
            sampling_counts[accession] += len(untrimmed_contigs_for_accession)

        # Add untrimmed contigs to filtered_contigs
        filtered_contigs_ids.extend(untrimmed_contig_ids)

        # Calculate the number of remaining contigs needed for this species
        remaining_contigs_needed = num - len(untrimmed_contig_ids)

        # Perform weighted sampling for the remaining contigs
        remaining_contigs = [contig for contig in species_contigs if "_NT" not in contig]
        if remaining_contigs_needed > 0 and remaining_contigs:
            weights = [
                prob
                for acc, prob in accessions
                for _ in range(len([contig for contig in remaining_contigs if acc in contig]))
            ]
            fixed_weights = [w / sum(weights) for w in weights]
            sampled_indices = np.random.choice(
                len(remaining_contigs),
                size=min(remaining_contigs_needed, len(remaining_contigs)),
                replace=False,
                p=fixed_weights,
            )
            sampled_contigs = [remaining_contigs[i] for i in sampled_indices]
            filtered_contigs_ids.extend(sampled_contigs)
            for contig in sampled_contigs:
                for accession, _ in accessions:
                    if accession in contig:
                        sampling_counts[accession] += 1

        # Store the filtered contigs and sampling counts for the current virus name
        all_filtered_contigs[virus_name] = filtered_contigs_ids
        all_sampling_counts[virus_name] = sampling_counts

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Save all_filtered_contigs as a JSON file using mpu.io
    filtered_contigs_path = os.path.join(output_dir, "filtered_contigs.json")
    mpu.io.write(filtered_contigs_path, all_filtered_contigs)

    # Save all_sampling_counts as a JSON file using mpu.io
    sampling_counts_path = os.path.join(output_dir, "sampling_counts.json")
    mpu.io.write(sampling_counts_path, all_sampling_counts)

    selected_contigs = [contig for sublist in all_filtered_contigs.values() for contig in sublist]

    return selected_contigs


def calculate_identity(seq1, seq2):
    """Calculate percentage identity between two sequences."""
    matches = sum(a == b for a, b in zip(seq1, seq2) if a != "-" and b != "-")
    length = min(len(seq1.replace("-", "")), len(seq2.replace("-", "")))
    return (matches / length) * 100 if length > 0 else 0


def reverse_complement_record(record):
    """
    Returns the reverse complement of a SeqRecord, retaining metadata.

    Parameters:
        record (SeqRecord): The input SeqRecord object.

    Returns:
        SeqRecord: A new SeqRecord object representing the reverse complement.
    """
    rc = record.reverse_complement()
    rc.id = record.id
    rc.name = record.name
    rc.description = record.description
    return rc


def translate_seq(seq_rec: SeqRecord):
    """
    Translate a nucleotide sequence in all three reading frames.

    Parameters:
    seq_rec (SeqRecord): The input SeqRecord containing the nucleotide sequence.

    Returns:
    dict: A dictionary containing the translated protein sequences in all three reading frames.
    """
    seq = seq_rec.seq

    # First reading frame
    prot1 = seq[: len(seq) - len(seq) % 3].translate()

    # Second reading frame
    prot2 = seq[1 : len(seq) - (len(seq) - 1) % 3].translate()

    # Third reading frame
    prot3 = seq[2 : len(seq) - (len(seq) - 2) % 3].translate()

    protein1 = SeqRecord(prot1, id=f"{seq_rec.id}_aa_frame1", description=seq_rec.description)
    protein2 = SeqRecord(prot2, id=f"{seq_rec.id}_aa_frame2", description=seq_rec.description)
    protein3 = SeqRecord(prot3, id=f"{seq_rec.id}_aa_frame3", description=seq_rec.description)

    return {protein1.id: protein1, protein2.id: protein2, protein3.id: protein3}


def translate_seq_6frames(seq_rec: SeqRecord):
    """
    Translate a nucleotide sequence in all six reading frames (three forward and three reverse).

    Parameters:
    seq_rec (SeqRecord): The input SeqRecord containing the nucleotide sequence.

    Returns:
    dict: A dictionary containing the translated protein sequences in all six reading frames.
    """
    seq = seq_rec.seq

    # Forward strand translations
    prot1 = seq[: len(seq) - len(seq) % 3].translate()  # First reading frame
    prot2 = seq[1 : len(seq) - (len(seq) - 1) % 3].translate()  # Second reading frame
    prot3 = seq[2 : len(seq) - (len(seq) - 2) % 3].translate()  # Third reading frame

    # Reverse strand translations
    rev_seq = reverse_complement_record(seq_rec).seq
    rev_prot1 = rev_seq[: len(rev_seq) - len(rev_seq) % 3].translate()  # First reading frame
    rev_prot2 = rev_seq[1 : len(rev_seq) - (len(rev_seq) - 1) % 3].translate()  # Second reading frame
    rev_prot3 = rev_seq[2 : len(rev_seq) - (len(rev_seq) - 2) % 3].translate()  # Third reading frame

    protein1 = SeqRecord(prot1, id=f"{seq_rec.id}_aa_frame1", description=seq_rec.description)
    protein2 = SeqRecord(prot2, id=f"{seq_rec.id}_aa_frame2", description=seq_rec.description)
    protein3 = SeqRecord(prot3, id=f"{seq_rec.id}_aa_frame3", description=seq_rec.description)
    rev_protein1 = SeqRecord(rev_prot1, id=f"{seq_rec.id}_aa_frame4", description=seq_rec.description)
    rev_protein2 = SeqRecord(rev_prot2, id=f"{seq_rec.id}_aa_frame5", description=seq_rec.description)
    rev_protein3 = SeqRecord(rev_prot3, id=f"{seq_rec.id}_aa_frame6", description=seq_rec.description)

    return {
        protein1.id: protein1,
        protein2.id: protein2,
        protein3.id: protein3,
        rev_protein1.id: rev_protein1,
        rev_protein2.id: rev_protein2,
        rev_protein3.id: rev_protein3,
    }


def run_orfipy(contig, log_file_dir):
    """
    Run ORFipy on a given contig and return the predicted ORFs.

    Parameters:
    contig (SeqRecord): The input contig sequence.
    log_file_dir (str): Directory to save the log file.

    Returns:
    dict: A dictionary of predicted ORFs.
    """
    with tempfile.TemporaryDirectory() as tempdir:
        # Write the contig to a temporary FASTA file
        fasta_path = os.path.join(tempdir, "input.fasta")
        SeqIO.write(contig, fasta_path, "fasta")

        # Ensure the log file directory exists
        os.makedirs(log_file_dir, exist_ok=True)
        log_file_path = os.path.join(log_file_dir, "log_file")

        # Run ORFipy and capture the output in the log file
        orfipy_command = f"orfipy {fasta_path} --min 300 --pep PEP --bed BED --between-stops --outdir {tempdir} > {log_file_path} 2>&1"
        os.system(orfipy_command)

        # Parse the predicted ORFs from the output file
        orfs = SeqIO.to_dict(SeqIO.parse(os.path.join(tempdir, "PEP"), "fasta"))

    return orfs


def align_mafft(s1, s2):
    """
    Align two sequences using MAFFT and return the alignment.

    Parameters:
    s1 (str): The first sequence to align.
    s2 (str): The second sequence to align.

    Returns:
    list: A list containing the alignment result.
    """
    try:
        with tempfile.TemporaryDirectory() as tempdir:
            fin_name = os.path.join(tempdir, "input.fasta")
            fout_name = os.path.join(tempdir, "output.fasta")
            records = [SeqRecord(Seq(s1), id="seqA"), SeqRecord(Seq(s2), id="seqB")]
            SeqIO.write(records, fin_name, "fasta")

            mafft_command = f"mafft --quiet {fin_name} > {fout_name}"
            result = subprocess.run(mafft_command, shell=True, capture_output=True, text=True)

            if result.returncode != 0:
                # Log the full input sequences to a file for further inspection
                with open("mafft_error_input.fasta", "w") as error_file:
                    SeqIO.write(records, error_file, "fasta")
                raise RuntimeError(f"Error in MAFFT: {result.stderr}")

            aln = AlignIO.read(fout_name, "fasta")
            aln.seqA = str(aln._records[0].seq.upper())
            aln.seqB = str(aln._records[1].seq.upper())
        return [aln]

    except Exception as e:
        print(f"Error aligning sequences: {s1[:30]}..., {s2[:30]}...")
        print(f"Exception: {e}")
        raise


def compute_identity_score_orfxref(args):
    """
    Compute the identity score between two sequences using MAFFT alignment.

    Parameters:
    args (tuple): A tuple containing two SeqRecord objects.

    Returns:
    dict: A dictionary containing the identity score and other alignment metrics.
    """
    s1, s2 = args

    for aln in align_mafft(s1.seq, s2.seq):
        score, M, N, aln_len, seqA_len, seqB_len, n_aln_len, gaps, gap_ratio, aln_orf_len = score_alignment(aln)
        break

    return {
        "orf_name": s1.id.replace("=", "_"),
        "ref_name": s2.id,
        "identity_score": score,
        "M": M,
        "N": N,
        "aln_len": aln_len,
        "orf_len": seqA_len,
        "ref_len": seqB_len,
        "N/aln_len": n_aln_len,
        "gap_openings": gaps,
        "gap_ratio": gap_ratio,
        "aln_orf_len": aln_orf_len,
    }


def process_in_chunks(args_list, chunk_size: int, max_workers: int):
    """
    Process a list of arguments in chunks using a process pool executor.

    Parameters:
    args_list (list): List of arguments to process.
    chunk_size (int): Size of each chunk.
    max_workers (int): Maximum number of workers for the process pool.

    Returns:
    list: A list of results from processing the chunks.
    """
    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        for i in range(0, len(args_list), chunk_size):
            chunk = args_list[i : i + chunk_size]
            results.extend(tqdm(executor.map(compute_identity_score_orfxref, chunk), total=len(chunk)))
    return results


def gap_openings(aln):
    """
    Calculate the number of gap openings in an alignment.

    Parameters:
    aln (MultipleSeqAlignment): The alignment object containing the sequences.

    Returns:
    int: The total number of gap openings in the alignment.
    """
    gapsA = sum(1 for i in range(1, len(aln.seqA)) if aln.seqA[i] == "-" and aln.seqA[i - 1] != "-")
    gapsB = sum(1 for i in range(1, len(aln.seqB)) if aln.seqB[i] == "-" and aln.seqB[i - 1] != "-")

    # Adjust for gaps at the end of the sequences
    if aln.seqA.endswith("-") and gapsA > 0:
        gapsA -= 1
    if aln.seqB.endswith("-") and gapsB > 0:
        gapsB -= 1

    return gapsA + gapsB


def score_alignment(aln):
    """
    Calculate the alignment score and various metrics for a given alignment.

    Parameters:
    aln (MultipleSeqAlignment): The alignment object containing the sequences.

    Returns:
    tuple: A tuple containing the alignment score and various metrics.
    """
    M = sum(n1 != n2 for n1, n2 in zip(aln.seqA, aln.seqB) if n1 != "-" and n2 != "-")
    N = sum(1 for n1, n2 in zip(aln.seqA, aln.seqB) if n1 != "-" and n2 != "-")
    seqA_len = sum(1 for n1 in aln.seqA if n1 != "-")
    seqB_len = sum(1 for n2 in aln.seqB if n2 != "-")
    aln_len, aln_orf_len = calculate_alignment_length(aln)
    gaps = gap_openings(aln)
    score = 1 - M / N
    return score, M, N, aln_len, seqA_len, seqB_len, N / aln_len, gaps, gaps / aln_len, aln_orf_len


def calculate_alignment_length(aln):
    """
    Calculate the alignment length and ORF alignment length for a given alignment.

    Parameters:
    aln (MultipleSeqAlignment): The alignment object containing the sequences.

    Returns:
    tuple: A tuple containing the alignment length and ORF alignment length.
    """
    match_positions = [i for i, (n1, n2) in enumerate(zip(aln.seqA, aln.seqB)) if n1 != "-" and n2 != "-"]
    if match_positions:
        first_match = match_positions[0]
        last_match = match_positions[-1] + 1
        aln_len = last_match - first_match
    else:
        aln_len = 0
    aln_orf_len = aln_len - sum(1 for n1 in aln.seqA[first_match:last_match] if n1 == "-")
    return aln_len, aln_orf_len


def calculate_histogram(group, num_bins: int):
    """
    Calculate histogram counts for a given group.

    Parameters:
    group (pd.DataFrame): The input DataFrame containing probabilities.
    num_bins (int): The number of bins for the histogram.

    Returns:
    pd.Series: A Series containing the histogram counts.
    """
    bins = np.linspace(0, 1, num_bins + 1)
    counts, _ = np.histogram(group["prob_1"], bins=bins)
    return pd.Series(counts, index=bins[:-1])


def calculate_histogram_cumsum(group, num_bins: int):
    """
    Calculate cumulative sum of histogram counts for a given group.

    Parameters:
    group (pd.DataFrame): The input DataFrame containing probabilities.
    num_bins (int): The number of bins for the histogram.

    Returns:
    pd.Series: A Series containing the cumulative sum of histogram counts.
    """
    bins = np.linspace(0, 1, num_bins + 1)
    counts, _ = np.histogram(group["prob_1"], bins=bins)
    return pd.Series(counts, index=bins[:-1]).cumsum()


def calculate_histogram_cumsum_reverse(group, num_bins: int):
    """
    Calculate reverse cumulative sum of histogram counts for a given group.

    Parameters:
    group (pd.DataFrame): The input DataFrame containing probabilities.
    num_bins (int): The number of bins for the histogram.

    Returns:
    pd.Series: A Series containing the reverse cumulative sum of histogram counts.
    """
    bins = np.linspace(0, 1, num_bins + 1)
    counts, _ = np.histogram(group["prob_1"], bins=bins)

    # Compute right-to-left cumulative sum
    reversed_cumsum = np.cumsum(counts[::-1])[::-1]

    return pd.Series(reversed_cumsum, index=bins[:-1])


def sampling(seq_record, n, lens_freq_path, rng: np.random.Generator):
    """
    Generate a sample of trimmed sequences based on a frequency distribution of lengths.

    Parameters:
    - seq_record: SeqRecord object, the sequence to be sampled.
    - n: Number of sequences to generate.
    - rng: np.random.Generator, random number generator for reproducibility.

    Returns:
    - contigs_dict: Dictionary of sampled sequences.
    - num_contigs: Number of sampled sequences.
    """
    # Load the frequency distribution of lengths
    lens_freq = mpu.io.read(lens_freq_path)
    trimmed_seqs = []

    # Get the sequence length
    seq_length = len(seq_record.seq)

    # Check if the sequence is larger than 8000 nucleotides
    if seq_length <= 8000:
        # Add one copy of the full untrimmed sequence if length <= 8000
        untrimmed_seq = SeqRecord(
            seq_record.seq,
            id=f"{seq_record.id}_start-0_len-{seq_length}_NT",
            description=seq_record.description,
        )
        trimmed_seqs.append(untrimmed_seq)  # Add the untrimmed copy first

    # Generate random lengths and trim sequences accordingly
    random_lengths = rng.choice(a=list(lens_freq.keys()), p=list(lens_freq.values()), size=n).astype(int)
    for rand_length in random_lengths:
        if rand_length < seq_length:
            start = rng.integers(0, seq_length - rand_length)
            trimmed_seq = seq_record.seq[start : start + rand_length]
            trimmed_seqs.append(
                SeqRecord(
                    trimmed_seq,
                    id=f"{seq_record.id}_start-{start}_len-{rand_length}",
                    description=seq_record.description,
                )
            )

    # Convert list of SeqRecords to a dictionary where keys are the sequence IDs
    contigs_dict = {trimmed_seq.id: trimmed_seq for trimmed_seq in trimmed_seqs}

    return contigs_dict, len(contigs_dict)


def subsample_contigs_dicts(outdir, list_of_dicts: list, selected_n: int, rng: np.random.Generator):
    """
    Subsample contigs to the same number as the genome with the least contigs.

    Parameters:
    outdir (str): Directory to save the report.
    list_of_dicts (list): List of dictionaries containing contigs.
    selected_n (int): Number of contigs to sample.
    rng (np.random.Generator): Random number generator for reproducibility.

    Returns:
    tuple: A tuple containing the combined dictionary of sampled contigs and the minimum dictionary size.
    """
    # Create a combined dictionary to hold the final result
    date = datetime.now().strftime("%Y-%m-%d")
    combined_dict = {}
    report_filename = f"results/{outdir}/sampling/{date}_sampled_contigs_report.txt"

    # Open a file to write the report
    with open(report_filename, "w") as report_file:
        # Write header for the report
        report_file.write("Contig Subsampling Report\n")
        report_file.write("=========================\n")

        # List to store the sizes of dictionaries
        sizes = []

        # Iterate over each dictionary and subsample
        for i, d in enumerate(list_of_dicts):
            dict_size = len(d)  # Get the size of the current dictionary
            sizes.append(dict_size)  # Add the size to the list

            # Calculate the sample size (ensure it doesn't exceed the dictionary size)
            sample_size = min(selected_n, dict_size)

            # Ensure at least one untrimmed contig is included
            untrimmed_contig_key = None
            for key in d.keys():
                if "_NT" in key:  # Identify the untrimmed contig by its identifier
                    untrimmed_contig_key = key
                    break

            sampled_keys = set()
            if untrimmed_contig_key:
                combined_dict[untrimmed_contig_key] = d[untrimmed_contig_key]
                sampled_keys.add(untrimmed_contig_key)
                sample_size -= 1  # Subtract one from the total number of contigs to sample

            # Randomly sample the remaining keys (excluding the untrimmed contig)
            available_keys = [key for key in d.keys() if key not in sampled_keys]
            if available_keys:  # Only sample if there are keys left to sample
                additional_sampled_keys = rng.choice(available_keys, size=sample_size, replace=False)
                for key in additional_sampled_keys:
                    combined_dict[key] = d[key]

        # Report the overall sampling information
        report_file.write(f"Selected sampling_n: {selected_n}\n")
        report_file.write(f"All contig lengths: {sizes}\n")
        report_file.write(f"Median dictionary size: {np.median(sizes)}\n")
        report_file.write(f"Mean dictionary size: {np.mean(sizes)}\n")
        report_file.write(f"Minimum dictionary size: {min(sizes)}\n")

    return combined_dict, min(sizes)


def plot_sampled_contigs_lens(outdir, subsampled, n):
    """
    Plot the length distribution of sampled contigs.

    Parameters:
    outdir (str): Directory to save the plot.
    subsampled (dict): Dictionary containing the sampled contigs.
    n (int): Number of sampled contigs.

    Returns:
    None
    """
    contig_lengths = [len(seq.seq) for seq in subsampled.values()]

    # Plot the lengths of the sampled contigs
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(contig_lengths, bins=20, color="skyblue", edgecolor="black")
    ax.set_title(f"Length distribution of sampled contigs (n={n})", fontsize=16)
    ax.set_xlabel("Contig lengths", fontsize=14)
    ax.set_ylabel("Frequency", fontsize=14)
    ax.grid(axis="y", linestyle="--", alpha=0.7)
    plt.tight_layout()

    # Save the plot
    date = datetime.now().strftime("%Y-%m-%d")
    save_path = f"results/{outdir}/sampling/{date}_sampled_contigs_{n}_length_distribution.png"
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path)
    plt.close(fig)


# def train_lr_and_predict_tuned(morf_predictions):
#     bin_df = prepare_bin_df_refs(morf_predictions, "stacking", num_bins=10)
#     X = bin_df.drop(columns=["ground_truth"])
#     y = bin_df["ground_truth"]

#     scorer = make_scorer(accuracy_score)
#     base_model = LogisticRegression(max_iter=1000)
#     model = TunedThresholdClassifierCV(base_model, scoring=scorer, cv=StratifiedKFold(n_splits=10))

#     model.fit(X, y)
#     best_score = model.best_threshold_

#     return "stacking", model, best_score


# def create_confusion_matrix(df, label_col="label", prediction_col="prediction"):
#     """
#     Create and display a confusion matrix for binary classification.

#     Parameters:
#     - df (pd.DataFrame): The input DataFrame with 'label' and 'prediction' columns.
#     - label_col (str): The name of the label column. ('tobamo' = 1, 'outgroup' & 'random' = 0)
#     - prediction_col (str): The name of the prediction column, containing binary predictions (0 or 1).

#     Returns:
#     - None: Displays the confusion matrix.
#     """
#     # Map the label column to binary values: 'tobamo' = 1, 'outgroup' and 'random' = 0
#     label_mapping = {"tobamo": 1, "outgroup": 0, "random": 0}
#     df["binary_label"] = df[label_col].map(label_mapping)

#     # Extract the true labels and predictions
#     y_true = df["binary_label"]
#     y_pred = df[prediction_col]

#     # Generate the confusion matrix
#     cm = confusion_matrix(y_true, y_pred)

#     # Display the confusion matrix
#     disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=["Outgroup/Random", "Tobamo"])
#     disp.plot(cmap="Blues", values_format="d")

# def filter_by_identity(input_file, output_file, threshold, min_length):
#     """
#     Filter sequences with less than the specified percent identity
#     and remove sequences shorter than the minimum length.
#     """
#     records = list(SeqIO.parse(input_file, "fasta"))
#     filtered = []

#     for rec1 in records:
#         # Skip sequences shorter than the minimum length
#         if len(rec1.seq) < min_length:
#             continue

#         unique = True
#         for rec2 in filtered:
#             identity = calculate_identity(str(rec1.seq), str(rec2.seq))
#             if identity > threshold:
#                 unique = False
#                 break
#         if unique:
#             filtered.append(rec1)

#     SeqIO.write(filtered, output_file, "fasta")
#     print(
#         f"After filtering by {threshold}% identity and removing sequences shorter than {min_length} nucleotides, "
#         f"{len(filtered)} out of {len(records)} sequences remain."
#     )

# def train_lr_and_predict(all_predictions):
#     models = {}
#     methods = ["stacking", "cumsum", "cumsum_reverse"]

#     for method_name in methods:
#         bin_df = prepare_bin_df_refs(all_predictions, method_name, num_bins=10)
#         X = bin_df.drop(columns=["ground_truth"])
#         y = bin_df["ground_truth"]
#         model_lr = LogisticRegression(max_iter=1000)
#         model_lr.fit(X, y)

#         models[method_name] = model_lr

#     return models
