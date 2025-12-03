#!/usr/bin/env python3
import random
from pathlib import Path

def read_fasta(path):
    seq = []
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            seq.append(line.strip().upper())
    return ''.join(seq)

def write_fastq(reads, qualities, out_path):
    with open(out_path, 'w') as f:
        for i, (hdr, seq) in enumerate(reads, start=1):
            qual = qualities.get(hdr, 'I' * len(seq))  # Q40
            f.write(f"@{hdr}\n{seq}\n+\n{qual}\n")

def compute_num_reads(genome_len, read_len, target_coverage):
    total_bases = target_coverage * genome_len
    return max(1, total_bases // read_len)

def circular_subseq(seq, start, length):
    n = len(seq)
    if start + length <= n:
        return seq[start:start+length]
    # wrap-around
    k = (start + length) % n
    return seq[start:] + seq[:k]

def simulate_circular(ref_seq, read_len=5000, coverage=20, seed=42, include_origin_spanners=True):
    random.seed(seed)
    n = len(ref_seq)
    num_reads = compute_num_reads(n, read_len, coverage)

    reads = []
    qualities = {}
    # Uniformly tile starts across the circle
    step = max(1, n // num_reads)
    starts = [(i * step) % n for i in range(num_reads)]

    for idx, s in enumerate(starts):
        rseq = circular_subseq(ref_seq, s, read_len)
        hdr = f"circ_{idx}_{s}"
        reads.append((hdr, rseq))
        qualities[hdr] = 'I' * len(rseq)

    if include_origin_spanners:
        # Force some origin-spanning reads around index 0 and near end
        for off in (0, n - (read_len // 2), n - 100, n - 1):
            rseq = circular_subseq(ref_seq, off % n, read_len)
            hdr = f"circ_origin_{off % n}"
            reads.append((hdr, rseq))
            qualities[hdr] = 'I' * len(rseq)

    return reads, qualities

def simulate_linearized(ref_seq, cut_index, read_len=5000, coverage=20, seed=123, include_cut_flank=True):
    random.seed(seed)
    # Linearize: ref[cut:] + ref[:cut]
    linear = ref_seq[cut_index:] + ref_seq[:cut_index]
    n = len(linear)
    num_reads = compute_num_reads(n, read_len, coverage)

    reads = []
    qualities = {}
    step = max(1, n // num_reads)
    starts = [(i * step) for i in range(num_reads - 2)]  # reserve slots for flanks

    for idx, s in enumerate(starts):
        if s + read_len > n:
            s = max(0, n - read_len)  # avoid truncation
        rseq = linear[s:s+read_len]
        hdr = f"lin_{idx}_{s}"
        reads.append((hdr, rseq))
        qualities[hdr] = 'I' * len(rseq)

    if include_cut_flank:
        # Reads that abut the cut ends (start of linearized and its end)
        flank_positions = [0, max(0, n - read_len)]
        for j, s in enumerate(flank_positions):
            rseq = linear[s:s+read_len]
            hdr = f"lin_cutflank_{j}_{s}"
            reads.append((hdr, rseq))
            qualities[hdr] = 'I' * len(rseq)

    return reads, qualities, linear

def main():
    ref_fasta = "plasmid.fasta"   # set your path
    out_dir = Path("synthetic_reads")
    out_dir.mkdir(exist_ok=True)

    ref_seq = read_fasta(ref_fasta)
    print(f"Loaded reference length: {len(ref_seq)} bp")

    # Parameters you can tune
    read_len = 5000       # typical long-read length
    coverage = 20         # total coverage across the plasmid
    cut_index = 1234      # index of restriction cut in reference coordinates

    # Circular dataset
    circ_reads, circ_quals = simulate_circular(ref_seq, read_len=read_len, coverage=coverage, seed=42)
    write_fastq(circ_reads, circ_quals, out_dir / "plasmid_circular.fastq")
    print(f"Wrote {len(circ_reads)} circular reads")

    # Linearized dataset
    lin_reads, lin_quals, linear_seq = simulate_linearized(ref_seq, cut_index=cut_index,
                                                           read_len=read_len, coverage=coverage, seed=123)
    write_fastq(lin_reads, lin_quals, out_dir / "plasmid_linearized.fastq")
    print(f"Wrote {len(lin_reads)} linearized reads")

    # Also write the linearized reference for alignment
    with open(out_dir / "plasmid_linearized.fasta", "w") as f:
        f.write(">plasmid_linearized\n")
        # Write single-line sequence without wrapping
        f.write(linear_seq + "\n")
    print("Wrote linearized reference FASTA")

if __name__ == "__main__":
    main()
