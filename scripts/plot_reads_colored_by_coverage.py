import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pysam

# --------
# Load reference sequence
# --------
def read_fasta_one_seq(path):
    seq = []
    with open(path, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq.append(line.strip().upper())
    return "".join(seq)

# --------
# Compute per-base coverage from BAM
# --------
def compute_coverage(bam_path, ref_name, ref_length):
    bam = pysam.AlignmentFile(bam_path, "rb")
    coverage = np.zeros(ref_length, dtype=np.int32)
    for pileupcolumn in bam.pileup(ref_name, 0, ref_length, truncate=True):
        coverage[pileupcolumn.reference_pos] = pileupcolumn.nsegments
    bam.close()
    return coverage

def normalize_coverage(cov):
    return (cov / cov.max()) * 100.0 if cov.max() > 0 else np.zeros_like(cov)

# --------
# Map coverage to colors
# --------
def build_colormap():
    return mcolors.LinearSegmentedColormap.from_list("red_blue", ["red", "blue"])

def coverage_to_color_array(coverage, cmap):
    norm = mcolors.Normalize(vmin=0, vmax=100)
    return cmap(norm(coverage))

# --------
# Plot each read as an arc colored by coverage
# --------
def plot_reads_by_coverage(bam_path, ref_seq, coverage_colors, output_path):
    bam = pysam.AlignmentFile(bam_path, "rb")
    ref_length = len(ref_seq)
    theta = np.linspace(0, 2 * np.pi, ref_length, endpoint=False)
    x = np.cos(theta)
    y = np.sin(theta)

    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw={"aspect": "equal"})
    ax.set_axis_off()

    for read in bam.fetch():
        if read.is_unmapped:
            continue
        ref_positions = read.get_reference_positions()
        for pos in ref_positions:
            if pos < 0 or pos >= ref_length:
                continue
            color = coverage_colors[pos]
            ax.plot([0, x[pos]], [0, y[pos]], color=color, linewidth=0.5)

    ax.set_title("Reads colored by coverage (0 → red, 100 → blue)", pad=12)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved plot: {output_path}")
    bam.close()

# --------
# Main
# --------
def main():
    fasta_path = "synthetic_reads/pChlamy_L4417_PCR_Product.fasta"
    bam_path   = "out/linear_run.sorted.bam"
    ref_name   = "pChlamy_SmaI-_L4417_PCR_Product"  # must match the FASTA header
    output_path = "out/reads_by_coverage.png"


    ref_seq = read_fasta_one_seq(fasta_path)
    ref_length = len(ref_seq)

    coverage = compute_coverage(bam_path, ref_name, ref_length)
    coverage_norm = normalize_coverage(coverage)
    cmap = build_colormap()
    coverage_colors = coverage_to_color_array(coverage_norm, cmap)

    plot_reads_by_coverage(bam_path, ref_seq, coverage_colors, output_path)

if __name__ == "__main__":
    main()
