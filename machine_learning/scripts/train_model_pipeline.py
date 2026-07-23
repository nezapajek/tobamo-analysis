#!/usr/bin/env python3

from utils import *
import argparse
import os
import sys
import pandas as pd
import numpy as np
from joblib import dump, load
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score, recall_score, precision_score

try:
    from joblib.externals.loky.process_executor import TerminatedWorkerError
except Exception:
    TerminatedWorkerError = None


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Train and evaluate tobamovirus classification models")
    parser.add_argument("input_df", help="Path to input DataFrame CSV")
    parser.add_argument("refs", help="Path to references Excel file")
    parser.add_argument("contigs", help="Path to contigs FASTA file")
    parser.add_argument("--outdir", default="default", help="Output directory name")
    parser.add_argument(
        "--stage", choices=["select", "evaluate", "final"], default="final", help="Pipeline stage to run"
    )
    parser.add_argument("--iterations", type=int, default=30, help="Number of iterations for cross-validation")
    parser.add_argument("--sample_depth", type=int, default=30, help="Number of contigs to sample per species")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=2,
        help="Number of parallel workers for model fitting (-1 uses all cores, default: 2)",
    )
    parser.add_argument("--bin-num", type=int, default=10, help="Number of bins for stacking model (default: 10)")
    parser.add_argument(
        "--threshold", type=float, default=None, help="Custom classification threshold (overrides CV results)"
    )
    parser.add_argument(
        "--use_fixed_threshold",
        action="store_true",
        help="In evaluate stage, use fixed threshold (--threshold or 0.5) instead of CV-optimized thresholds",
    )
    parser.add_argument(
        "--bins",
        type=int,
        nargs="+",
        default=[10],
        help="Bin counts to evaluate in --stage evaluate (default: 10)",
    )
    return parser.parse_args()


def check_inputs(input_df_path, refs_path, contig_fasta_path):
    """Check if input files exist and are in correct format"""
    check_file_exists(input_df_path, "CSV")
    check_file_exists(refs_path, "Excel")
    check_file_exists(contig_fasta_path, "FASTA")

    return (
        pd.read_csv(input_df_path, index_col=0),
        pd.read_excel(refs_path),
        SeqIO.to_dict(SeqIO.parse(contig_fasta_path, "fasta")),
    )


def run_grid_search_with_retry(model, params, X_train, y_train, n_jobs, pre_dispatch, cv=5, scoring="accuracy"):
    """Run GridSearchCV and fall back to serial execution if loky workers terminate unexpectedly."""
    grid_search = GridSearchCV(
        model,
        param_grid=params,
        n_jobs=n_jobs,
        pre_dispatch=pre_dispatch,
        cv=cv,
        scoring=scoring,
    )

    try:
        grid_search.fit(X_train, y_train)
        return grid_search
    except Exception as exc:
        is_terminated_worker_error = (
            TerminatedWorkerError is not None and isinstance(exc, TerminatedWorkerError)
        ) or "A worker process managed by the executor was unexpectedly terminated" in str(exc)

        if not is_terminated_worker_error or n_jobs == 1:
            raise

        print("    Warning: loky worker terminated during GridSearchCV. Retrying this model with n_jobs=1...")

        fallback_grid_search = GridSearchCV(
            model,
            param_grid=params,
            n_jobs=1,
            pre_dispatch=1,
            cv=cv,
            scoring=scoring,
        )
        fallback_grid_search.fit(X_train, y_train)
        return fallback_grid_search


def model_selection_old(input_df, refs, contigs, outdir="model_selection", random_seed=42, n_jobs=2):
    """Perform model selection and hyperparameter tuning"""
    print("Starting model selection and hyperparameter tuning...")
    print(f"Parallel workers: {n_jobs}")

    grid_pre_dispatch = n_jobs if n_jobs > 0 else "2*n_jobs"

    # Define the models and their parameter grids
    models = {
        "LogisticRegression": {
            "model": LogisticRegression(max_iter=10000),
            "params": {
                "C": [0.01, 0.1, 1, 10],
                "solver": ["liblinear", "saga"],
                "penalty": ["l1", "l2"],
                "class_weight": ["balanced", None],
            },
        },
        "RandomForest": {
            "model": RandomForestClassifier(n_jobs=1, random_state=random_seed),
            "params": {"n_estimators": [50, 100, 150, 200, 300], "max_depth": [5, 10, 20, 40, 50, None]},
        },
        "SVM": {"model": SVC(), "params": {"C": [10, 50, 100], "kernel": ["linear", "rbf", "poly"]}},
        "DecisionTree": {"model": DecisionTreeClassifier(), "params": {"max_depth": [None, 5, 10]}},
        "NaiveBayes": {"model": GaussianNB(), "params": {}},
        "KNN": {"model": KNeighborsClassifier(), "params": {"n_neighbors": [3, 5, 7, 15]}},
    }

    os.makedirs(f"results/{outdir}/report", exist_ok=True)

    # Create 5-fold CV split
    folds = stratified_kfold_split(refs, n_splits=5, random_state=random_seed)

    performance_metrics_list = []

    for idx, (train_refs, test_refs) in enumerate(folds):
        print(f"Processing fold {idx+1}/5")

        # Prepare TRAIN and TEST data
        train_all, test_all = prepare_train_test(input_df, train_refs, test_refs)

        # Subsample contigs
        selected_training_contigs = subsample_contigs(
            train_refs,
            contigs,
            num=30,
            output_dir=f"results/{outdir}/report/{idx}_train_report",
            random_seed=random_seed,
        )
        selected_test_contigs = subsample_contigs(
            test_refs, contigs, num=30, output_dir=f"results/{outdir}/report/{idx}_test_report", random_seed=random_seed
        )

        # Filter datasets to selected contigs
        train = train_all[train_all["contig_name"].isin(selected_training_contigs)]
        test = test_all[test_all["contig_name"].isin(selected_test_contigs)]

        # Use only forward ORFs for training
        train = train[train["strand"] == "FORWARD"]

        # Prepare features and target
        X_train = train.drop(columns=["orf_type", "strand", "virus_name", "accession", "contig_name"])
        y_train = (train["orf_type"] == "tobamo").astype(int)
        X_test = test.drop(columns=["orf_type", "strand", "virus_name", "accession", "contig_name"])
        y_test = (test["orf_type"] == "tobamo").astype(int)

        # Standardize data
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)

        # Find best model
        best_model, best_score, best_params, best_model_name = None, -float("inf"), None, None

        for model_name, model_info in models.items():
            print(f"  Testing {model_name}...")
            model, params = model_info["model"], model_info["params"]

            grid_search = run_grid_search_with_retry(
                model,
                params,
                X_train,
                y_train,
                n_jobs=n_jobs,
                pre_dispatch=grid_pre_dispatch,
                cv=5,
                scoring="accuracy",
            )

            if grid_search.best_score_ > best_score:
                best_model = grid_search.best_estimator_
                best_score = grid_search.best_score_
                best_params = grid_search.best_params_
                best_model_name = model_name

        # Evaluate best model
        y_pred = best_model.predict(X_test)

        metrics = {
            "fold": idx,
            "model": best_model_name,
            "accuracy": accuracy_score(y_test, y_pred),
            "auc_roc": roc_auc_score(y_test, y_pred),
            "f1_score": f1_score(y_test, y_pred),
            "recall": recall_score(y_test, y_pred),
            "specificity": recall_score(y_test, y_pred, pos_label=0),
            "precision": precision_score(y_test, y_pred),
            "negative_predictive_value": precision_score(y_test, y_pred, pos_label=0),
            "best_params": best_params,
        }

        performance_metrics_list.append(metrics)

    # Save results
    pd.DataFrame(performance_metrics_list).to_csv(f"results/{outdir}/performance_metrics.csv")

    # Determine overall best model
    model_performances = {}
    for metric in performance_metrics_list:
        model_name = metric["model"]
        if model_name not in model_performances:
            model_performances[model_name] = []
        model_performances[model_name].append(metric["accuracy"])

    # Calculate average accuracy for each model
    avg_performances = {model: np.mean(scores) for model, scores in model_performances.items()}
    best_model = max(avg_performances, key=avg_performances.get)

    # Save best model information
    with open(f"results/{outdir}/best_model.txt", "w") as f:
        f.write(f"Best model: {best_model}\n")
        f.write(f"Average accuracy: {avg_performances[best_model]:.4f}\n")

    print(f"Model selection completed! Best model: {best_model}")


def model_selection(input_df, refs, contigs, outdir="model_selection", random_seed=42, iterations=5, n_jobs=2):
    """Perform model selection and hyperparameter tuning with multiple iterations"""
    print(f"Starting model selection and hyperparameter tuning with {iterations} iterations...")
    print(f"Parallel workers: {n_jobs}")

    grid_pre_dispatch = n_jobs if n_jobs > 0 else "2*n_jobs"

    # Define the models and their parameter grids
    models = {
        "LogisticRegression": {
            "model": LogisticRegression(max_iter=10000, random_state=random_seed),
            "params": {
                "C": [0.01, 0.1, 1, 10],
                "solver": ["liblinear", "saga"],
                "penalty": ["l1", "l2"],
                "class_weight": ["balanced", None],
            },
        },
        "RandomForest": {
            "model": RandomForestClassifier(n_jobs=1, random_state=random_seed),
            "params": {"n_estimators": [50, 100, 150, 200, 300], "max_depth": [5, 10, 20, 40, 50, None]},
        },
        "SVM": {"model": SVC(random_state=random_seed), "params": {"C": [10, 50, 100], "kernel": ["linear", "rbf", "poly"]}},
        "DecisionTree": {"model": DecisionTreeClassifier(random_state=random_seed), "params": {"max_depth": [None, 5, 10]}},
        "NaiveBayes": {"model": GaussianNB(), "params": {}},
        "KNN": {"model": KNeighborsClassifier(), "params": {"n_neighbors": [3, 5, 7, 15]}},
    }

    os.makedirs(f"results/{outdir}/report", exist_ok=True)

    # List to store results across all iterations
    all_performance_metrics = []

    # Run multiple iterations
    for iteration in range(iterations):
        print(f"\nStarting iteration {iteration+1}/{iterations}")
        iter_seed = random_seed + iteration

        # Create 5-fold CV split for this iteration
        folds = stratified_kfold_split(refs, n_splits=5, random_state=iter_seed)

        performance_metrics_list = []

        for idx, (train_refs, test_refs) in enumerate(folds):
            print(f"Processing fold {idx+1}/5")

            # Prepare TRAIN and TEST data
            train_all, test_all = prepare_train_test(input_df, train_refs, test_refs)

            # Subsample contigs
            selected_training_contigs = subsample_contigs(
                train_refs,
                contigs,
                num=30,
                output_dir=f"results/{outdir}/report/iter{iteration}_fold{idx}_train_report",
                random_seed=iter_seed,
            )
            selected_test_contigs = subsample_contigs(
                test_refs,
                contigs,
                num=30,
                output_dir=f"results/{outdir}/report/iter{iteration}_fold{idx}_test_report",
                random_seed=iter_seed,
            )

            # Filter datasets to selected contigs
            train = train_all[train_all["contig_name"].isin(selected_training_contigs)]
            test = test_all[test_all["contig_name"].isin(selected_test_contigs)]

            # Use only forward ORFs for training
            train = train[train["strand"] == "FORWARD"]

            # Prepare features and target
            X_train = train.drop(columns=["orf_type", "strand", "virus_name", "accession", "contig_name"])
            y_train = (train["orf_type"] == "tobamo").astype(int)
            X_test = test.drop(columns=["orf_type", "strand", "virus_name", "accession", "contig_name"])
            y_test = (test["orf_type"] == "tobamo").astype(int)

            # Standardize data
            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)
            X_test = scaler.transform(X_test)

            # Find best model
            best_model, best_score, best_params, best_model_name = None, -float("inf"), None, None

            for model_name, model_info in models.items():
                print(f"  Testing {model_name}...")
                model, params = model_info["model"], model_info["params"]

                grid_search = run_grid_search_with_retry(
                    model,
                    params,
                    X_train,
                    y_train,
                    n_jobs=n_jobs,
                    pre_dispatch=grid_pre_dispatch,
                    cv=5,
                    scoring="accuracy",
                )

                if grid_search.best_score_ > best_score:
                    best_model = grid_search.best_estimator_
                    best_score = grid_search.best_score_
                    best_params = grid_search.best_params_
                    best_model_name = model_name

            # Evaluate best model
            y_pred = best_model.predict(X_test)

            metrics = {
                "iteration": iteration,
                "fold": idx,
                "model": best_model_name,
                "accuracy": accuracy_score(y_test, y_pred),
                "auc_roc": roc_auc_score(y_test, y_pred),
                "f1_score": f1_score(y_test, y_pred),
                "recall": recall_score(y_test, y_pred),
                "specificity": recall_score(y_test, y_pred, pos_label=0),
                "precision": precision_score(y_test, y_pred),
                "negative_predictive_value": precision_score(y_test, y_pred, pos_label=0),
                "best_params": best_params,
            }

            performance_metrics_list.append(metrics)

        # Save this iteration's results
        iter_df = pd.DataFrame(performance_metrics_list)
        iter_df.to_csv(f"results/{outdir}/iter_{iteration}_performance_metrics.csv")

        # Append to overall results
        all_performance_metrics.extend(performance_metrics_list)

    # Save combined results from all iterations
    pd.DataFrame(all_performance_metrics).to_csv(f"results/{outdir}/all_performance_metrics.csv")

    # Calculate overall model performance across all iterations
    model_performances = {}
    for metric in all_performance_metrics:
        model_name = metric["model"]
        if model_name not in model_performances:
            model_performances[model_name] = []
        model_performances[model_name].append(metric["accuracy"])

    # Calculate average and std of accuracy for each model
    avg_performances = {
        model: {"mean": np.mean(scores), "std": np.std(scores), "count": len(scores)}
        for model, scores in model_performances.items()
    }

    # Find best model based on average performance
    best_model = max(avg_performances.keys(), key=lambda m: avg_performances[m]["mean"])

    # Count how many times each model was selected as best in individual folds
    model_selection_counts = {}
    for metric in all_performance_metrics:
        model_name = metric["model"]
        if model_name not in model_selection_counts:
            model_selection_counts[model_name] = 0
        model_selection_counts[model_name] += 1

    # Save summary information
    with open(f"results/{outdir}/model_selection_summary.txt", "w") as f:
        f.write(f"Model Selection Summary (across {iterations} iterations)\n")
        f.write("=" * 50 + "\n\n")

        f.write("Overall Best Model: {}\n".format(best_model))
        f.write(
            "Average accuracy: {:.4f} ± {:.4f}\n\n".format(
                avg_performances[best_model]["mean"], avg_performances[best_model]["std"]
            )
        )

        f.write("Performance by Model:\n")
        f.write("-" * 40 + "\n")
        for model, stats in avg_performances.items():
            f.write(
                "{}: {:.4f} ± {:.4f} (selected in {}/{} folds)\n".format(
                    model,
                    stats["mean"],
                    stats["std"],
                    model_selection_counts[model],
                    iterations * 5,  # total number of folds across all iterations
                )
            )

    print(f"Model selection completed across {iterations} iterations! Best model: {best_model}")
    return best_model, avg_performances


def train_and_evaluate(
    input_df,
    refs,
    contigs,
    outdir="cv_evaluation",
    iterations=30,
    sample_depth=30,
    random_seed=42,
    selected_model=None,
    use_fixed_threshold=False,
    fixed_threshold=0.5,
    bins_to_evaluate=None,
):
    """Perform extensive cross-validation with both prediction methods"""
    print(f"Starting comprehensive evaluation with {iterations} iterations...")
    threshold_mode = "fixed" if use_fixed_threshold else "tuned"
    print(f"Threshold mode: {threshold_mode} ({fixed_threshold if use_fixed_threshold else 'CV-optimized'})")
    bins_to_evaluate = [10] if bins_to_evaluate is None else list(bins_to_evaluate)
    print(f"Evaluating binned prediction bins: {bins_to_evaluate}")

    os.makedirs(f"results/{outdir}", exist_ok=True)

    # Lists to store all results
    all_extreme_predictions = []
    all_binned_predictions = []
    all_orf_predictions = []

    for iteration in range(iterations):
        print(f"Starting iteration {iteration+1}/{iterations}")
        iter_seed = random_seed + iteration

        # Create 5-fold CV split
        folds = stratified_kfold_split(refs, n_splits=5, random_state=iter_seed)

        for idx, (train_refs, test_refs) in enumerate(folds):
            fold_dir = f"results/{outdir}/report/iter{iteration}_fold{idx}"
            os.makedirs(fold_dir, exist_ok=True)

            # Prepare TRAIN and TEST data
            train_all, test_all = prepare_train_test(input_df, train_refs, test_refs)

            # Subsample contigs
            selected_training_contigs = subsample_contigs(
                train_refs, contigs, num=sample_depth, output_dir=f"{fold_dir}/train_report", random_seed=iter_seed
            )
            selected_test_contigs = subsample_contigs(
                test_refs, contigs, num=sample_depth, output_dir=f"{fold_dir}/test_report", random_seed=iter_seed
            )

            # Filter datasets
            train = train_all[train_all["contig_name"].isin(selected_training_contigs)]
            test = test_all[test_all["contig_name"].isin(selected_test_contigs)]

            # Use only forward ORFs for training
            train = train[train["strand"] == "FORWARD"]

            ################################################ TRAIN PHASE

            # Make ORF predictions: Direct Random Forest
            morf_predictions = train_rf_and_predict(train, selected_model)  # make predictions using LOOCV

            # Train Logistic Regression model on LOOCV prediction
            mc_dict = train_lr_and_predict_binned_test(morf_predictions, bins_to_evaluate)

            # Train RandomForest model on all training data
            morf, sorf, _ = train_rf_on_all_data(train, selected_model)

            ################################################ TEST phase

            # Predict test ORFs using the trained RandomForest model
            test_orf_predictions = predict_orfs(test, morf, sorf, refs=True)

            test_orf_predictions["iteration"] = iteration
            test_orf_predictions["fold"] = idx
            all_orf_predictions.append(test_orf_predictions)

            # Method 1: Select contig class based on most extreme ORF probability
            extreme_predictions = filter_extreme_probability(test_orf_predictions, idx, refs_=True)
            extreme_predictions["iteration"] = iteration
            all_extreme_predictions.append(extreme_predictions)

            # Method 2: Predict contigs using logistic regression (based on binned ORF predictions)
            for mc_name_n, mc_info in mc_dict.items():
                mc_name = mc_name_n.split("_")[0]
                num = int(mc_name_n.split("_")[1])
                mc = mc_info["model"]  # Extract the model
                best_threshold = fixed_threshold if use_fixed_threshold else mc_info["best_threshold"]

                final_predictions = predict_contigs_binned_test(test_orf_predictions, mc, idx, mc_name, num, refs_=True)
                final_predictions["n"] = num
                final_predictions["threshold"] = best_threshold
                final_predictions["threshold_mode"] = threshold_mode
                final_predictions["threshold_source"] = "fixed_user" if use_fixed_threshold else "cv_optimized"
                final_predictions["random_seed"] = random_seed
                final_predictions["iteration"] = iteration
                all_binned_predictions.append(final_predictions)

    # Combine and save all results
    extreme_results = pd.concat(all_extreme_predictions)
    extreme_results.to_csv(f"results/{outdir}/extreme_predictions_results.csv", index=False)

    binned_results = pd.concat(all_binned_predictions)
    binned_results.to_csv(f"results/{outdir}/binned_predictions_results.csv", index=False)
    threshold_suffix = f"fixed_{str(fixed_threshold).replace('.', 'p')}" if use_fixed_threshold else "tuned"
    binned_results.to_csv(f"results/{outdir}/binned_predictions_results_{threshold_suffix}.csv", index=False)

    orf_results = pd.concat(all_orf_predictions)
    orf_results.to_csv(f"results/{outdir}/orf_predictions_results.csv", index=False)

    # Calculate summary performance metrics.
    # With multiple bins, evaluate each stacking bin as its own method.
    methods = {"extreme": extreme_results}
    if not binned_results.empty:
        for n in sorted(binned_results["n"].unique()):
            methods[f"binned_{int(n)}"] = binned_results[binned_results["n"] == n].copy()

    summary_metrics = []

    for method_name, results_df in methods.items():
        # Add predicted class based on method threshold.
        # Stacking predictions include per-row threshold; extreme method defaults to 0.5.
        if "threshold" in results_df.columns:
            results_df["predicted_class"] = np.where(results_df["prob_1"] >= results_df["threshold"], 1, 0)
        else:
            results_df["predicted_class"] = np.where(results_df["prob_1"] >= 0.5, 1, 0)

        # Group by iteration to calculate per-iteration metrics
        iteration_metrics = []

        for iteration, group in results_df.groupby("iteration"):
            iter_metrics = {
                "iteration": iteration,
                "method": method_name,
                "accuracy": accuracy_score(group["ground_truth"], group["predicted_class"]),
                "f1": f1_score(group["ground_truth"], group["predicted_class"]),
                "precision": precision_score(group["ground_truth"], group["predicted_class"]),
                "recall": recall_score(group["ground_truth"], group["predicted_class"]),
                "auc": roc_auc_score(group["ground_truth"], group["prob_1"]),
            }
            iteration_metrics.append(iter_metrics)

        # Create DataFrame of per-iteration metrics
        iter_metrics_df = pd.DataFrame(iteration_metrics)

        # Calculate statistics across iterations
        metrics_stats = {
            "method": method_name,
            "accuracy_mean": iter_metrics_df["accuracy"].mean(),
            "accuracy_std": iter_metrics_df["accuracy"].std(),
            "accuracy_ci95_low": iter_metrics_df["accuracy"].mean()
            - 1.96 * iter_metrics_df["accuracy"].std() / np.sqrt(len(iter_metrics_df)),
            "accuracy_ci95_high": iter_metrics_df["accuracy"].mean()
            + 1.96 * iter_metrics_df["accuracy"].std() / np.sqrt(len(iter_metrics_df)),
            "f1_mean": iter_metrics_df["f1"].mean(),
            "f1_std": iter_metrics_df["f1"].std(),
            "f1_ci95_low": iter_metrics_df["f1"].mean()
            - 1.96 * iter_metrics_df["f1"].std() / np.sqrt(len(iter_metrics_df)),
            "f1_ci95_high": iter_metrics_df["f1"].mean()
            + 1.96 * iter_metrics_df["f1"].std() / np.sqrt(len(iter_metrics_df)),
            "precision_mean": iter_metrics_df["precision"].mean(),
            "precision_std": iter_metrics_df["precision"].std(),
            "precision_ci95_low": iter_metrics_df["precision"].mean()
            - 1.96 * iter_metrics_df["precision"].std() / np.sqrt(len(iter_metrics_df)),
            "precision_ci95_high": iter_metrics_df["precision"].mean()
            + 1.96 * iter_metrics_df["precision"].std() / np.sqrt(len(iter_metrics_df)),
            "recall_mean": iter_metrics_df["recall"].mean(),
            "recall_std": iter_metrics_df["recall"].std(),
            "recall_ci95_low": iter_metrics_df["recall"].mean()
            - 1.96 * iter_metrics_df["recall"].std() / np.sqrt(len(iter_metrics_df)),
            "recall_ci95_high": iter_metrics_df["recall"].mean()
            + 1.96 * iter_metrics_df["recall"].std() / np.sqrt(len(iter_metrics_df)),
            "auc_mean": iter_metrics_df["auc"].mean(),
            "auc_std": iter_metrics_df["auc"].std(),
            "auc_ci95_low": iter_metrics_df["auc"].mean()
            - 1.96 * iter_metrics_df["auc"].std() / np.sqrt(len(iter_metrics_df)),
            "auc_ci95_high": iter_metrics_df["auc"].mean()
            + 1.96 * iter_metrics_df["auc"].std() / np.sqrt(len(iter_metrics_df)),
        }

        # Save the per-iteration metrics
        iter_metrics_df.to_csv(f"results/{outdir}/{method_name}_iteration_metrics.csv", index=False)

        # Add summary metrics
        summary_metrics.append(metrics_stats)

    # Save summary metrics
    summary_df = pd.DataFrame(summary_metrics)
    summary_df.to_csv(f"results/{outdir}/method_comparison_stats.csv", index=False)

    # Create a simplified summary with mean ± std
    simplified_summary = []
    for _, row in summary_df.iterrows():
        method = row["method"]
        simplified = {
            "method": method,
            "accuracy": f"{row['accuracy_mean']:.3f} ± {row['accuracy_std']:.3f}",
            "f1": f"{row['f1_mean']:.3f} ± {row['f1_std']:.3f}",
            "precision": f"{row['precision_mean']:.3f} ± {row['precision_std']:.3f}",
            "recall": f"{row['recall_mean']:.3f} ± {row['recall_std']:.3f}",
            "auc": f"{row['auc_mean']:.3f} ± {row['auc_std']:.3f}",
        }
        simplified_summary.append(simplified)

    # Save simplified summary
    pd.DataFrame(simplified_summary).to_csv(f"results/{outdir}/method_comparison_simplified.csv", index=False)

    # Determine best method based on mean accuracy
    best_method = summary_df.loc[summary_df["accuracy_mean"].idxmax()]["method"]

    # Save best method information
    with open(f"results/{outdir}/best_method.txt", "w") as f:
        f.write(f"Best method: {best_method}\n")
        best_metrics = summary_df[summary_df["method"] == best_method].iloc[0]
        for metric in ["accuracy", "f1", "precision", "recall", "auc"]:
            mean = best_metrics[f"{metric}_mean"]
            std = best_metrics[f"{metric}_std"]
            ci_low = best_metrics[f"{metric}_ci95_low"]
            ci_high = best_metrics[f"{metric}_ci95_high"]
            f.write(f"{metric}: {mean:.4f} ± {std:.4f} (95% CI: {ci_low:.4f}-{ci_high:.4f})\n")

    # Add threshold analysis
    print("\nAnalyzing threshold optimization results...")
    generate_threshold_summary(binned_results, outdir)

    print(f"\nEvaluation completed! Best method: {best_method}")
    print(f"Results saved to: results/{outdir}/")


def train_final_model(
    input_df, refs, contigs, outdir="final_model", bin_num=10, fixed_threshold=0.5, random_seed=42, selected_model=None
):
    """Train the final RF+LR stacking model on all training data

    Parameters:
    -----------
    input_df : pandas.DataFrame
        Input DataFrame with features
    refs : pandas.DataFrame
        References dataframe
    contigs : dict
        Dictionary of contigs from FASTA
    outdir : str, default="final_model"
        Output directory name
    bin_num : int, default=10
        Number of bins to use for the stacking model (determined in evaluation)
    fixed_threshold : float, default=0.5
        Classification threshold. If 0.5, uses default. Otherwise uses value from CV.
    random_seed : int, default=42
        Random seed for reproducibility
    selected_model : dict, optional
        Selected model configuration
    """

    print(f"Training final model using RF + LR stacking approach with {bin_num} bins...")
    print(f"Using classification threshold: {fixed_threshold}")

    os.makedirs(f"results/{outdir}", exist_ok=True)

    # Keep only training data
    train_refs = refs.loc[refs["training"] == 1]
    train_refs.loc[:, "sampling_prob"] = train_refs["length"] / train_refs.groupby("virus_name")["length"].transform(
        "sum"
    )

    # Subsample contigs
    selected_contigs = subsample_contigs(
        train_refs, contigs, num=30, output_dir=f"results/{outdir}/sampling_report", random_seed=random_seed
    )

    # Filter to selected contigs
    train = input_df[input_df["contig_name"].isin(selected_contigs)]

    # Use only forward ORFs for training
    train_fwd = train[train["strand"] == "FORWARD"]

    # Train RandomForest model
    morf, sorf, features = train_rf_on_all_data(train_fwd, selected_model)

    # Save RF model and features
    dump(morf, f"results/{outdir}/rf_model.joblib")
    dump(sorf, f"results/{outdir}/rf_scaler.joblib")
    pd.Series(features).to_csv(f"results/{outdir}/rf_feature_names.csv", index=False, header=False)

    # Train RF on subset and get predictions with LOOCV for LR model training
    morf_predictions = train_rf_and_predict(train_fwd, selected_model)

    # Save morf_predictions for reproducibility and analysis
    morf_predictions.to_csv(f"results/{outdir}/morf_predictions.csv", index=False)
    print(f"Saved morf_predictions with {len(morf_predictions)} ORF predictions")

    # Train Logistic Regression stacking model with selected bins num
    mc_dict = train_lr_final_model(morf_predictions, bin_num, fixed_threshold)
    binned_model = mc_dict["model"]
    threshold = mc_dict["threshold"]
    bin_df = mc_dict["bin_df"]

    # Save bin_df for reproducibility and analysis
    bin_df.to_csv(f"results/{outdir}/bin_df_{bin_num}.csv", index=False)
    print(f"Saved bin_df with {len(bin_df)} binned contig predictions")

    # Save LR stacking model
    dump(binned_model, f"results/{outdir}/lr_binned_{bin_num}_model.joblib")
    print(f"Saved binned prediction model with {bin_num} bins")

    # Save threshold info
    with open(f"results/{outdir}/threshold_info.txt", "w") as f:
        f.write(f"Classification threshold: {threshold}\n")
        f.write(f"Source: {'default' if threshold == 0.5 else 'cross_validation'}\n")

    print(f"Saved binned prediction model with {bin_num} bins and threshold {threshold}")

    # Save feature importances
    feature_importances = morf.feature_importances_
    feature_names = train_fwd.drop(columns=["orf_type", "strand", "virus_name", "accession", "contig_name"]).columns
    feature_importances_df = pd.DataFrame({"Feature": feature_names, "Importance": feature_importances})
    feature_importances_df = feature_importances_df.sort_values("Importance", ascending=False)
    feature_importances_df.to_csv(f"results/{outdir}/feature_importances.csv", index=False)

    # Save all and top 40 features for easy reference
    feature_importances_df.to_csv(f"results/{outdir}/all_features.csv", index=False)
    top_features = feature_importances_df.head(40)
    top_features.to_csv(f"results/{outdir}/top_40_features.csv", index=False)

    print(f"Final models (RF + LR stacking with {bin_num} bins) trained and saved!")


def main():
    args = parse_args()
    print(f"Starting pipeline at stage: {args.stage}")

    if args.n_jobs == 0:
        raise ValueError("--n-jobs cannot be 0. Use -1 for all cores or a positive integer.")

    # Check inputs and load data
    input_df, refs, contigs = check_inputs(args.input_df, args.refs, args.contigs)
    print("Input data loaded successfully")

    # Default selected model
    default_selected_model = {
        "model": RandomForestClassifier(n_estimators=200, max_depth=50, n_jobs=args.n_jobs, random_state=args.seed)
    }

    # Run the requested pipeline stage
    if args.stage == "select":
        model_selection(
            input_df,
            refs,
            contigs,
            outdir=args.outdir,
            random_seed=args.seed,
            iterations=args.iterations,
            n_jobs=args.n_jobs,
        )

    elif args.stage == "evaluate":
        fixed_threshold = args.threshold if args.threshold is not None else 0.5
        if not 0 <= fixed_threshold <= 1:
            raise ValueError("--threshold must be between 0 and 1")
        bins_to_evaluate = sorted(set(args.bins))
        if any(n <= 0 for n in bins_to_evaluate):
            raise ValueError("--bins values must be positive integers")

        train_and_evaluate(
            input_df,
            refs,
            contigs,
            outdir=args.outdir,
            iterations=args.iterations,
            sample_depth=args.sample_depth,
            random_seed=args.seed,
            selected_model=default_selected_model,
            use_fixed_threshold=args.use_fixed_threshold,
            fixed_threshold=fixed_threshold,
            bins_to_evaluate=bins_to_evaluate,
        )

    elif args.stage == "final":
        train_final_model(
            input_df,
            refs,
            contigs,
            outdir=args.outdir,
            bin_num=args.bin_num,
            fixed_threshold=args.threshold if args.threshold is not None else 0.5,
            random_seed=args.seed,
            selected_model=default_selected_model,
        )

    print("Pipeline completed successfully!")


if __name__ == "__main__":
    main()
