import shlex
import subprocess
from pathlib import Path


def run_command(command):
    """Run a shell command and check for errors."""
    result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return result.stdout, result.stderr


def main():
    references_dir = Path(__file__).resolve().parent
    repo_root = references_dir.parents[1]
    # NOTE: this script spans two repos. diamond_db and select_contigs_script live in a sibling
    # tobamo-snakemake checkout -- diamond_db (resources/tpdb2.dmnd) is a build artifact not
    # committed to either repo (build via `diamond makedb` from the tobamoviral protein set, see
    # the paper's Methods); adjust snakemake_root below if your checkout lives elsewhere.
    snakemake_root = repo_root.parent / "tobamo-snakemake"

    input_fasta = (
        repo_root / "data/non-virga_representatives/ncbi_virus_refseq_30k_dedupl_seqs_20250116.fasta"
    )
    results_dir = references_dir / "results/non-virga"
    output_fasta = results_dir / "non-virga_tpdb2_diamond_selected.fasta"
    daa_file = results_dir / "non-virga_tpdb2.daa"
    diamond_info_tsv = results_dir / "non-virga_tpdb2_diamond_info.tsv"
    diamond_db = snakemake_root / "resources/tpdb2.dmnd"
    select_contigs_script = snakemake_root / "workflow/scripts/select_contigs.py"

    for required_path in (input_fasta, diamond_db, select_contigs_script):
        if not required_path.exists():
            raise FileNotFoundError(f"Required file not found: {required_path}")

    results_dir.mkdir(parents=True, exist_ok=True)

    # Run diamond blastx
    print("Running diamond blastx...")
    blastx_command = (
        "diamond blastx --threads 32 "
        f"-d {shlex.quote(str(diamond_db))} "
        f"-k 20 -q {shlex.quote(str(input_fasta))} "
        f"--daa {shlex.quote(str(daa_file))}"
    )
    run_command(blastx_command)

    # Convert .daa file to tabular format
    print("Converting .daa file to tabular format...")
    view_command = (
        f"diamond view --daa {shlex.quote(str(daa_file))} " f"--outfmt 6 > {shlex.quote(str(diamond_info_tsv))}"
    )
    run_command(view_command)

    # Run select_contigs.py script
    print("Running select_contigs.py script...")
    select_contigs_command = (
        f"python {shlex.quote(str(select_contigs_script))} "
        f"{shlex.quote(str(input_fasta))} "
        f"{shlex.quote(str(diamond_info_tsv))} "
        f"{shlex.quote(str(output_fasta))}"
    )
    run_command(select_contigs_command)

    print(f"Pipeline completed. Output FASTA file: {output_fasta}")


if __name__ == "__main__":
    main()
