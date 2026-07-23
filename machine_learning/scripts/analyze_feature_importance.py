#!/usr/bin/env python3

import argparse
import json
import os
from pathlib import Path

import pandas as pd


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Summarize Random Forest feature importance results from all_features.csv"
    )
    parser.add_argument("importance_csv", help="Path to all_features.csv (columns: Feature,Importance)")
    parser.add_argument(
        "--outdir",
        default=None,
        help="Optional output directory for summary files (default: no files written)",
    )
    parser.add_argument(
        "--top-k",
        nargs="+",
        type=int,
        default=[5, 10, 20],
        help="Top-k cutoffs for cumulative importance summary (default: 5 10 20)",
    )
    parser.add_argument(
        "--print-top",
        type=int,
        default=10,
        help="Number of top features to print and export (default: 10)",
    )
    return parser.parse_args()


def validate_input(df):
    """Validate required columns and value types."""
    required_columns = {"Feature", "Importance"}
    missing = required_columns - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    if df.empty:
        raise ValueError("Input file has no feature rows")

    df["Importance"] = pd.to_numeric(df["Importance"], errors="raise")
    if (df["Importance"] < 0).any():
        raise ValueError("Importance values must be non-negative")


def make_share(raw_value, total):
    """Return raw value and normalized share (0-1)."""
    if total <= 0:
        return raw_value, 0.0
    return raw_value, raw_value / total


def summarize_feature_importance(df, top_k):
    """Compute global and grouped feature importance summaries."""
    df_sorted = df.sort_values("Importance", ascending=False).reset_index(drop=True)
    total_importance = float(df_sorted["Importance"].sum())

    summary = {
        "n_features": int(len(df_sorted)),
        "total_importance": total_importance,
        "top_k_cumulative": {},
        "source_groups": {},
        "protein_groups": {},
        "feature_families": {},
    }

    for k in sorted(set(top_k)):
        if k <= 0:
            continue
        cumulative = float(df_sorted.head(k)["Importance"].sum())
        raw, share = make_share(cumulative, total_importance)
        summary["top_k_cumulative"][str(k)] = {
            "importance_sum": raw,
            "share": share,
        }

    source_masks = {
        "tobamo": df_sorted["Feature"].str.startswith("tobamo_"),
        "outgroup": df_sorted["Feature"].str.startswith("outgroup_"),
    }
    generic_mask = ~(source_masks["tobamo"] | source_masks["outgroup"])
    source_masks["generic"] = generic_mask

    for name, mask in source_masks.items():
        raw = float(df_sorted.loc[mask, "Importance"].sum())
        raw, share = make_share(raw, total_importance)
        summary["source_groups"][name] = {"importance_sum": raw, "share": share}

    protein_masks = {
        "orf1": df_sorted["Feature"].str.contains(r"^(?:tobamo|outgroup)_orf1_", regex=True),
        "orf2": df_sorted["Feature"].str.contains(r"^(?:tobamo|outgroup)_orf2_", regex=True),
        "cp": df_sorted["Feature"].str.contains(r"^(?:tobamo|outgroup)_cp_", regex=True),
    }

    for name, mask in protein_masks.items():
        raw = float(df_sorted.loc[mask, "Importance"].sum())
        raw, share = make_share(raw, total_importance)
        summary["protein_groups"][name] = {"importance_sum": raw, "share": share}

    family_masks = {
        "identity_score": df_sorted["Feature"].str.contains("identity_score", regex=False),
        "gap_metrics": df_sorted["Feature"].str.contains(r"gap_ratio|gap_openings", regex=True),
        "N_over_aln_len": df_sorted["Feature"].str.contains("N/aln_len", regex=False),
        "alignment_length_any": df_sorted["Feature"].str.contains(r"aln_len|aln_orf_len", regex=True),
    }

    for name, mask in family_masks.items():
        raw = float(df_sorted.loc[mask, "Importance"].sum())
        raw, share = make_share(raw, total_importance)
        summary["feature_families"][name] = {"importance_sum": raw, "share": share}

    return df_sorted, summary


def print_summary(summary, top_df):
    """Print concise human-readable summary."""
    print("Feature Importance Summary")
    print("=" * 28)
    print(f"n_features: {summary['n_features']}")
    print(f"total_importance: {summary['total_importance']:.6f}")

    print("\nTop-k cumulative importance:")
    for k, values in summary["top_k_cumulative"].items():
        print(f"  top{k}: {values['importance_sum']:.6f} ({values['share']:.2%})")

    print("\nSource groups:")
    for name, values in summary["source_groups"].items():
        print(f"  {name}: {values['importance_sum']:.6f} ({values['share']:.2%})")

    print("\nProtein groups:")
    for name, values in summary["protein_groups"].items():
        print(f"  {name}: {values['importance_sum']:.6f} ({values['share']:.2%})")

    print("\nFeature families:")
    for name, values in summary["feature_families"].items():
        print(f"  {name}: {values['importance_sum']:.6f} ({values['share']:.2%})")

    print("\nTop features:")
    for idx, row in top_df.iterrows():
        print(f"  {idx + 1}. {row['Feature']} ({row['Importance']:.6f})")


def write_outputs(outdir, summary, top_df):
    """Write summary files for downstream reporting."""
    out_path = Path(outdir)
    out_path.mkdir(parents=True, exist_ok=True)

    summary_json = out_path / "feature_importance_summary.json"
    with open(summary_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    top_csv = out_path / "top_features_summary.csv"
    top_df.to_csv(top_csv, index=False)

    group_rows = []
    for group_type in ["source_groups", "protein_groups", "feature_families"]:
        for group_name, values in summary[group_type].items():
            group_rows.append(
                {
                    "group_type": group_type,
                    "group_name": group_name,
                    "importance_sum": values["importance_sum"],
                    "share": values["share"],
                }
            )
    pd.DataFrame(group_rows).to_csv(out_path / "feature_group_summary.csv", index=False)

    print(f"\nSaved summary files to: {out_path}")


def main():
    args = parse_args()

    if not os.path.exists(args.importance_csv):
        raise FileNotFoundError(f"Input file not found: {args.importance_csv}")

    df = pd.read_csv(args.importance_csv)
    validate_input(df)

    df_sorted, summary = summarize_feature_importance(df, args.top_k)
    top_df = df_sorted.head(max(args.print_top, 1)).copy()

    print_summary(summary, top_df)

    if args.outdir:
        write_outputs(args.outdir, summary, top_df)


if __name__ == "__main__":
    main()
