# scripts/plot_plasmid_colormap.py
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# --------
# 1) Load FASTA sequence
# --------
def read_fasta_one_seq(path):
    seq = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue
            seq.append(line.upper())
    return "".join(seq)

# --------
# 2) Build a per-base variable in [0, 100]
#    Replace 'variable' with real coverage later if desired
# --------
def synth_variable(n):
    # Example synthetic variable: sum of smooth waves normalized to 0–100
    x = np.linspace(0, 2 * np.pi, n, endpoint=False)
    v = (
        40 * (0.5 + 0.5 * np.sin(1.0 * x)) +
        35 * (0.5 + 0.5 * np.sin(3.0 * x + 0.8)) +
        25 * (0.5 + 0.5 * np.sin(7.0 * x + 2.2))
    )
    # Clip to 0–100 range for safety
    v = np.clip(v, 0, 100)
    return v

# --------
# 3) Map variable -> colors (red to blue)
# --------
def build_colormap():
    # Red -> Blue linear gradient
    return mcolors.LinearSegmentedColormap.from_list("red_blue", ["red", "blue"])

def variable_to_colors(variable, cmap):
    norm = mcolors.Normalize(vmin=0, vmax=100)
    return cmap(norm(variable))

# --------
# 4) Plot circular ring with per-base color
# --------
def plot_ring(variable_colors, radius=1.0, line_width=2.0, figure_size=6.0):
    n = len(variable_colors)
    theta = np.linspace(0, 2 * np.pi, n + 1)  # n segments
    # Compute x,y points on the circle
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)

    fig, ax = plt.subplots(figsize=(figure_size, figure_size), subplot_kw={"aspect": "equal"})
    ax.set_axis_off()

    # Draw colored segments around the circle
    for i in range(n):
        ax.plot(
            [x[i], x[i + 1]],
            [y[i], y[i + 1]],
            color=variable_colors[i],
            linewidth=line_width,
            solid_capstyle="butt"
        )

    # Optional title to indicate what the color encodes
    ax.set_title("Circular plasmid colored by variable (0 → red, 100 → blue)", pad=12)
    return fig, ax

# --------
# 5) Main
# --------
def main():
    fasta_path = os.path.join("synthetic_reads", "pChlamy_L4417_PCR_Product.fasta")
    out_dir = "out"
    out_png = os.path.join(out_dir, "pChlamy_L4417_PCR_Product_colormap.png")

    os.makedirs(out_dir, exist_ok=True)

    seq = read_fasta_one_seq(fasta_path)
    n = len(seq)
    print(f"Loaded sequence length: {n}")

    variable = synth_variable(n)              # Replace with real coverage later
    cmap = build_colormap()
    colors = variable_to_colors(variable, cmap)

    fig, ax = plot_ring(colors, radius=1.0, line_width=2.0, figure_size=6.0)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"Saved figure: {out_png}")

if __name__ == "__main__":
    main()
