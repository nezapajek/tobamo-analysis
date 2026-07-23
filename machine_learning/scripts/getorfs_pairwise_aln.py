from utils import *

contigs = sys.argv[1]
outdir = sys.argv[2]
orientation = sys.argv[3]

aa_refs = "../data/tobamo/all_proteins.fasta"

# Check if the contigs file is a FASTA file
if not contigs.endswith((".fasta", ".fa", ".fna")):
    print("Error: The contigs file must be a FASTA file with extension .fasta, .fa, or .fna")
    sys.exit(1)

if orientation not in ["forward", "unknown"]:
    print("Error: The orientation must be either 'forward' or 'unknown'")
    sys.exit(1)

# Check if the output already exists
combined_orfs_path = f"results/{outdir}/orfs/combined_orfs.fasta"
if os.path.exists(combined_orfs_path):
    print(f"Skipping processing as {combined_orfs_path} already exists.")
else:
    print("Processing...")

    # Create the output directory
    os.makedirs(f"results/{outdir}/orfs", exist_ok=True)

    ### GET ORFS FROM CONTIGS
    # RUN ORFIPY
    os.system(
        f"orfipy {contigs} --min 300 --pep PEP --bed BED --between-stops --outdir results/{outdir}/orfs > results/{outdir}/orfs/logfile 2>&1"
    )

    # TRANSLATE contigs
    contigs = SeqIO.to_dict(SeqIO.parse(contigs, "fasta"))
    orfs = {}

    # Translate ORFs to 3 reading frames if forward or 6 reading frames if unknown
    if orientation == "forward":
        for contig in contigs.values():
            orfs.update(translate_seq(contig))
    elif orientation == "unknown":
        for contig in contigs.values():
            orfs.update(translate_seq_6frames(contig))

    out = f"results/{outdir}/orfs/translated_orfs.fasta"
    with open(out, "w") as f:
        SeqIO.write(orfs.values(), f, "fasta")

    # COMBINE ORFS
    os.system(f"cat results/{outdir}/orfs/PEP {out} > {combined_orfs_path}")
    print("ORFs extracted and translated successfully")


pairwise_aln = f"results/{outdir}/pairwise_aln.csv"
if os.path.exists(pairwise_aln):
    print(f"Skipping processing as {pairwise_aln} already exists.")
else:
    ### PAIRWISE ALIGNMENTS
    all_orfs = SeqIO.to_dict(SeqIO.parse(combined_orfs_path, "fasta"))
    aa = SeqIO.to_dict(SeqIO.parse(aa_refs, "fasta"))

    args_list = []
    for orf in all_orfs.values():
        for aa_ref in aa.values():
            args_list.append([orf, aa_ref])

    print(f"Computing pairwise alignments for {len(args_list)} pairs")

    results = process_in_chunks(args_list, 1000, 32)

    df = pd.DataFrame(results)
    df.to_csv(pairwise_aln, index=False)
    print("Pairwise alignments computed successfully")
