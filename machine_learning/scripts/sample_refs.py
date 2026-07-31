from utils import *

refs_path = sys.argv[1]
outdir = sys.argv[2]
sampling_size = int(sys.argv[3])
subsampling_n = int(sys.argv[4])
lens_freq_path = sys.argv[5]
rng = np.random.default_rng(42)

refs = SeqIO.to_dict(SeqIO.parse(refs_path, "fasta"))
os.makedirs(f"results/{outdir}/sampling", exist_ok=True)
date = datetime.now().strftime("%Y-%m-%d")

# Check if sequences are too short and filter them out
short_sequences = []
filtered_refs = {}
message_printed = False

for seq_id, seq in refs.items():
    if len(seq) < 800:
        short_sequences.append((seq_id, len(seq)))
        if not message_printed:
            print(
                f"Some sequences were shorter than 800bp and were removed. Removed sequences listed in short_sequences_report.txt"
            )
            message_printed = True
    else:
        filtered_refs[seq_id] = seq

# Write short sequences to a report
with open(f"results/{outdir}/sampling/{date}_short_sequences_report.txt", "w") as report_file:
    report_file.write("Sequence ID\tLength\n")
    for seq_id, length in short_sequences:
        report_file.write(f"{seq_id}\t{length}\n")

# Proceed with the filtered sequences
all_sampled_contigs = []
selected_sampling_contigs = {}
contigs_lens = []

for seq in filtered_refs.values():
    contigs, contigs_no = sampling(seq, sampling_size, lens_freq_path, rng)
    all_sampled_contigs.append(contigs)
    contigs_lens.append(contigs_no)

subsampled, min_contigs = subsample_contigs_dicts(outdir, all_sampled_contigs, subsampling_n, rng)


if min_contigs < subsampling_n:
    print(
        f"Warning: sampling size too small. {min_contigs} contigs were sampled. If you want to sample {subsampling_n} contigs, consider increasing the sampling size."
    )
    n = min_contigs
else:
    print(f"Sampling {subsampling_n} contigs.")
    n = subsampling_n

# make a plot of the sampled contigs lengths distribution
plot_sampled_contigs_lens(outdir, subsampled, n)

# write the sampled contigs to a file
output_file = f"results/{outdir}/sampling/{date}_sampled_contigs_{n}.fasta"
with open(output_file, "w") as f:
    SeqIO.write(subsampled.values(), f, "fasta")
