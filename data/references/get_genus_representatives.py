from __future__ import annotations

import argparse
import math
import os
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


GENERA = [
    "furovirus",
    "goravirus",
    "hordeivirus",
    "pecluvirus",
    "pomovirus",
    "tobamovirus",
    "tobravirus",
]


@dataclass
class PipelineConfig:
    input_dir: Path
    output_root: Path
    identity_threshold: float
    max_workers: int
    overwrite: bool
    genera: list[str]


def parse_args() -> PipelineConfig:
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent.parent
    default_input = repo_root / "data/domain_sci_input/ncbi_virgaviridae_20250120"
    default_output = script_dir / "results/virga"

    parser = argparse.ArgumentParser(
        description=(
            "Process Virgaviridae representative sequences by genus: deduplicate, "
            "compute pairwise identity, cluster, and generate dendrograms."
        )
    )
    parser.add_argument("--input-dir", type=Path, default=default_input)
    parser.add_argument("--output-root", type=Path, default=default_output)
    parser.add_argument("--identity-threshold", type=float, default=0.90)
    parser.add_argument("--max-workers", type=int, default=max(1, min(32, (os.cpu_count() or 1))))
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--genera", nargs="+", default=GENERA)

    args = parser.parse_args()
    return PipelineConfig(
        input_dir=args.input_dir,
        output_root=args.output_root,
        identity_threshold=args.identity_threshold,
        max_workers=max(1, args.max_workers),
        overwrite=args.overwrite,
        genera=args.genera,
    )


def ensure_dirs(root: Path) -> dict[str, Path]:
    paths = {
        "data": root / "data",
        "pairwise_aln": root / "pairwise_aln",
        "clusters": root / "clusters",
        "selected_species": root / "selected_species",
        "dendrograms": root / "dendrograms",
        "dendrograms_selected_species": root / "dendrograms_selected_species",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def is_non_empty(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def genus_completed(genus: str, paths: dict[str, Path]) -> bool:
    unique_fasta = paths["data"] / f"{genus}_RNA1_unique.fasta"
    mapper_file = paths["data"] / f"{genus}_RNA1_mapper.txt"
    pairwise_csv = paths["pairwise_aln"] / f"pairwise_aln_{genus}_RNA1_sequences.csv"
    cluster_csv = paths["clusters"] / f"{genus}_RNA1_identity_score_clusters.csv"
    selected_csv = paths["selected_species"] / f"{genus}_identity_score_RNA1_selected_species.csv"
    dendrogram_png = paths["dendrograms"] / f"{genus}_RNA1_identity_score_dendrogram.png"
    dendrogram_selected_png = paths["dendrograms_selected_species"] / f"{genus}_RNA1_identity_score_dendrogram.png"

    full_outputs = [
        unique_fasta,
        mapper_file,
        pairwise_csv,
        cluster_csv,
        selected_csv,
        dendrogram_png,
        dendrogram_selected_png,
    ]
    if all(is_non_empty(path) for path in full_outputs):
        return True

    short_path_outputs = [
        unique_fasta,
        mapper_file,
        selected_csv,
    ]
    return all(is_non_empty(path) for path in short_path_outputs)


def check_mafft_available() -> None:
    try:
        subprocess.run(["mafft", "--version"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception as error:
        raise RuntimeError("MAFFT is required but not available in PATH.") from error


def remove_duplicates(fin_path: Path, fout_path: Path, mapper_path: Path) -> list[SeqRecord]:
    records = list(SeqIO.parse(str(fin_path), "fasta"))
    unique_sequences: dict[str, list[str]] = {}
    for record in records:
        seq = str(record.seq)
        unique_sequences.setdefault(seq, []).append(record.id)

    unique_records = [SeqRecord(Seq(sequence), id=ids[0], description="") for sequence, ids in unique_sequences.items()]
    SeqIO.write(unique_records, str(fout_path), "fasta")

    with mapper_path.open("w") as mapper_handle:
        for ids in unique_sequences.values():
            mapper_handle.write(f"{ids[0]}: {', '.join(ids)}\n")

    return unique_records


def align_mafft(seq_a: str, seq_b: str) -> tuple[str, str]:
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        fin_name = temp_path / "input.fasta"
        fout_name = temp_path / "output.fasta"
        records = [SeqRecord(Seq(seq_a), id="seqA"), SeqRecord(Seq(seq_b), id="seqB")]
        SeqIO.write(records, str(fin_name), "fasta")

        with fout_name.open("w") as output_handle:
            subprocess.run(
                ["mafft", "--quiet", str(fin_name)], check=True, stdout=output_handle, stderr=subprocess.DEVNULL
            )

        aln = AlignIO.read(str(fout_name), "fasta")
        return str(aln[0].seq).upper(), str(aln[1].seq).upper()


def calculate_alignment_length(seq_a: str, seq_b: str) -> tuple[int, int]:
    match_positions = [index for index, (a, b) in enumerate(zip(seq_a, seq_b)) if a != "-" and b != "-"]
    if not match_positions:
        return 0, 0

    first_match = match_positions[0]
    last_match = match_positions[-1] + 1
    aln_len = last_match - first_match
    aln_orf_len = aln_len - sum(1 for value in seq_a[first_match:last_match] if value == "-")
    return aln_len, aln_orf_len


def gap_openings(seq_a: str, seq_b: str) -> int:
    gaps_a = sum(1 for i in range(1, len(seq_a)) if seq_a[i] == "-" and seq_a[i - 1] != "-")
    gaps_b = sum(1 for i in range(1, len(seq_b)) if seq_b[i] == "-" and seq_b[i - 1] != "-")

    if seq_a.endswith("-") and gaps_a > 0:
        gaps_a -= 1
    if seq_b.endswith("-") and gaps_b > 0:
        gaps_b -= 1
    return gaps_a + gaps_b


def score_alignment(seq_a: str, seq_b: str) -> dict[str, float]:
    mismatches = sum(1 for a, b in zip(seq_a, seq_b) if a != b and a != "-" and b != "-")
    comparable = sum(1 for a, b in zip(seq_a, seq_b) if a != "-" and b != "-")
    seq_a_len = sum(1 for value in seq_a if value != "-")
    seq_b_len = sum(1 for value in seq_b if value != "-")
    aln_len, aln_orf_len = calculate_alignment_length(seq_a, seq_b)
    gap_count = gap_openings(seq_a, seq_b)

    if comparable == 0:
        identity_score = 0.0
    else:
        identity_score = 1 - (mismatches / comparable)

    if aln_len == 0:
        n_aln_len = 0.0
        gap_ratio = 0.0
        clustering_score = 0.0
    else:
        n_aln_len = comparable / aln_len
        gap_ratio = gap_count / aln_len
        clustering_score = identity_score * (1 - gap_ratio)

    return {
        "identity_score": identity_score,
        "M": mismatches,
        "N": comparable,
        "aln_len": aln_len,
        "orf_len": seq_a_len,
        "ref_len": seq_b_len,
        "N/aln_len": n_aln_len,
        "gap_openings": gap_count,
        "gap_ratio": gap_ratio,
        "aln_orf_len": aln_orf_len,
        "clustering_score": clustering_score,
    }


def compute_pair(args: tuple[str, str, str, str]) -> dict[str, object]:
    seq1_id, seq1, seq2_id, seq2 = args
    aligned_a, aligned_b = align_mafft(seq1, seq2)
    scores = score_alignment(aligned_a, aligned_b)
    return {
        "orf_name": seq1_id.replace("=", "_"),
        "ref_name": seq2_id,
        **scores,
    }


def pairwise_refs(records: list[SeqRecord], max_workers: int) -> pd.DataFrame:
    record_pairs = [(a.id, str(a.seq), b.id, str(b.seq)) for a, b in combinations(records, 2)]
    if not record_pairs:
        return pd.DataFrame(
            columns=[
                "orf_name",
                "ref_name",
                "identity_score",
                "M",
                "N",
                "aln_len",
                "orf_len",
                "ref_len",
                "N/aln_len",
                "gap_openings",
                "gap_ratio",
                "aln_orf_len",
                "clustering_score",
            ]
        )

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        rows = list(executor.map(compute_pair, record_pairs))
    return pd.DataFrame(rows)


def prepend_length(identifier: str, length_mapper: dict[str, int]) -> str:
    return f"{length_mapper.get(identifier, 0)}_{identifier}"


def prep_df(pairwise_df: pd.DataFrame) -> pd.DataFrame:
    df = pairwise_df[["orf_name", "ref_name", "identity_score"]].copy()
    df.rename(columns={"orf_name": "ref1", "ref_name": "ref2"}, inplace=True)
    df["diff"] = 1 - df["identity_score"]
    return df


def make_distance_matrix(
    df: pd.DataFrame, length_mapper: dict[str, int]
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    df = df.copy()
    df["ref1"] = df["ref1"].apply(lambda value: prepend_length(value, length_mapper))
    df["ref2"] = df["ref2"].apply(lambda value: prepend_length(value, length_mapper))
    idx = pd.concat([df["ref1"], df["ref2"]]).unique()

    res = (
        pd.pivot(df, index=["ref1"], columns=["ref2"], values="diff")
        .reindex(index=idx, columns=idx)
        .fillna(0)
        .astype(float)
    )
    matrix = (res + res.T).to_numpy()
    return res, matrix, idx


def make_linkage(matrix: np.ndarray) -> np.ndarray:
    condensed_distance_matrix = ssd.squareform(matrix)
    return sch.linkage(condensed_distance_matrix, method="single", optimal_ordering=True)


def build_cluster_df(linkage_matrix: np.ndarray, idx: np.ndarray) -> pd.DataFrame:
    distances = linkage_matrix[:, 2]
    unique_distances = sorted(np.unique(distances))
    rows = []

    for distance in unique_distances:
        clusters = sch.fcluster(linkage_matrix, t=distance, criterion="distance")
        cluster_labels = {cluster_id: [] for cluster_id in np.unique(clusters)}
        for label, cluster_id in zip(idx, clusters):
            cluster_labels[cluster_id].append(label)
        rows.append({"distance": distance, "identity": 1 - distance, "clusters": len(cluster_labels)})

    return pd.DataFrame(rows)


def select_longest_species_per_cluster(
    linkage_matrix: np.ndarray,
    idx: np.ndarray,
    distance_threshold: float,
) -> pd.DataFrame:
    clusters = sch.fcluster(linkage_matrix, t=distance_threshold, criterion="distance")
    cluster_labels = {cluster_id: [] for cluster_id in np.unique(clusters)}
    for label, cluster_id in zip(idx, clusters):
        cluster_labels[cluster_id].append(label)

    data = []
    for cluster_id, species_list in cluster_labels.items():
        for species in species_list:
            length = int(species.split("_", 1)[0]) if "_" in species else 0
            data.append({"cluster": cluster_id, "species": species, "length": length})

    table = pd.DataFrame(data)
    selected = table.loc[table.groupby("cluster")["length"].idxmax()].sort_values("cluster").reset_index(drop=True)
    return selected


def save_dendrogram(
    linkage_matrix: np.ndarray,
    labels: list[str],
    output_path: Path,
    title: str,
    distance_threshold: float | None = None,
    highlighted_labels: set[str] | None = None,
) -> None:
    num_labels = len(labels)
    fig_height = max(10, num_labels * 0.3)
    plt.figure(figsize=(12, fig_height))

    if distance_threshold is None:
        sch.dendrogram(linkage_matrix, labels=labels, orientation="left")
    else:
        sch.dendrogram(
            linkage_matrix,
            labels=labels,
            orientation="left",
            color_threshold=distance_threshold,
            above_threshold_color="grey",
        )
        plt.axvline(x=distance_threshold, color="black", linestyle="--")

    plt.title(title)
    plt.xlabel("Distance")

    if highlighted_labels:
        axis = plt.gca()
        for tick in axis.get_ymajorticklabels():
            if tick.get_text() in highlighted_labels:
                tick.set_fontweight("bold")

    plt.subplots_adjust(left=0.5, right=0.9, top=0.9, bottom=0.1)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def process_genus(genus: str, config: PipelineConfig, paths: dict[str, Path]) -> None:
    input_fasta = config.input_dir / f"{genus}RNA1_sequences.fasta"
    if not input_fasta.exists():
        print(f"[{genus}] Missing input file: {input_fasta}")
        return

    if not config.overwrite and genus_completed(genus, paths):
        print(f"[{genus}] Already completed successfully. Skipping.")
        return

    print(f"[{genus}] Processing...")
    unique_fasta = paths["data"] / f"{genus}_RNA1_unique.fasta"
    mapper_file = paths["data"] / f"{genus}_RNA1_mapper.txt"
    pairwise_csv = paths["pairwise_aln"] / f"pairwise_aln_{genus}_RNA1_sequences.csv"
    cluster_csv = paths["clusters"] / f"{genus}_RNA1_identity_score_clusters.csv"
    selected_csv = paths["selected_species"] / f"{genus}_identity_score_RNA1_selected_species.csv"
    dendrogram_png = paths["dendrograms"] / f"{genus}_RNA1_identity_score_dendrogram.png"
    dendrogram_selected_png = paths["dendrograms_selected_species"] / f"{genus}_RNA1_identity_score_dendrogram.png"

    unique_records = remove_duplicates(input_fasta, unique_fasta, mapper_file)
    length_mapper = {record.id: len(record.seq) for record in unique_records}

    if len(unique_records) < 2:
        print(f"[{genus}] Skipping pairwise/clustering (need at least 2 sequences, got {len(unique_records)}).")
        pd.DataFrame(columns=["cluster", "species", "length"]).to_csv(selected_csv, index=False)
        return

    if pairwise_csv.exists() and not config.overwrite:
        pairwise_df = pd.read_csv(pairwise_csv)
        print(f"[{genus}] Using existing pairwise file: {pairwise_csv}")
    else:
        pairwise_df = pairwise_refs(unique_records, max_workers=config.max_workers)
        pairwise_df.to_csv(pairwise_csv, index=False)

    if pairwise_df.empty:
        print(f"[{genus}] Pairwise result is empty, skipping clustering.")
        return

    prep = prep_df(pairwise_df)
    res, matrix, idx = make_distance_matrix(prep, length_mapper)
    linkage_matrix = make_linkage(matrix)

    cluster_df = build_cluster_df(linkage_matrix, idx)
    cluster_df.to_csv(cluster_csv, index=False)

    eligible = cluster_df.loc[cluster_df["identity"] >= config.identity_threshold, "distance"]
    if eligible.empty:
        distance_threshold = float(cluster_df["distance"].median())
    else:
        distance_threshold = float(eligible.max())

    selected = select_longest_species_per_cluster(linkage_matrix, idx, distance_threshold)
    selected.to_csv(selected_csv, index=False)
    selected_labels = set(selected["species"].tolist())

    labels = list(res.index)
    save_dendrogram(
        linkage_matrix=linkage_matrix,
        labels=labels,
        output_path=dendrogram_png,
        title=f"Hierarchical Dendrogram of {genus} RNA1 (identity score)",
        distance_threshold=distance_threshold,
    )

    save_dendrogram(
        linkage_matrix=linkage_matrix,
        labels=labels,
        output_path=dendrogram_selected_png,
        title=(
            f"Hierarchical Dendrogram of {genus} RNA1 (selected representatives, "
            f"identity threshold={1 - distance_threshold:.2f})"
        ),
        distance_threshold=distance_threshold,
        highlighted_labels=selected_labels,
    )

    print(
        f"[{genus}] Done. Unique={len(unique_records)}, pairs={len(pairwise_df)}, "
        f"selected={len(selected)}, threshold={distance_threshold:.6f}"
    )


def main() -> None:
    config = parse_args()
    if not config.input_dir.exists():
        raise FileNotFoundError(f"Input directory does not exist: {config.input_dir}")

    if not (0.0 < config.identity_threshold <= 1.0):
        raise ValueError("--identity-threshold must be in (0, 1].")

    check_mafft_available()
    paths = ensure_dirs(config.output_root)

    print(f"Input dir: {config.input_dir}")
    print(f"Output root: {config.output_root}")
    print(f"Genera: {', '.join(config.genera)}")
    print(f"Identity threshold: {config.identity_threshold}")
    print(f"Max workers: {config.max_workers}")

    for genus in config.genera:
        process_genus(genus=genus.lower(), config=config, paths=paths)

    print("Finished processing all genera.")


if __name__ == "__main__":
    main()
