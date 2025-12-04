import math
from collections.abc import Iterable
from enum import Enum
from pathlib import Path

from matplotlib import cm
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import numpy.typing as npt

from concatmap.struct import PolarCoordinate
from concatmap.struct import PolarLineSegment
from concatmap.utils import PositionToAngleConverter

Array1D = npt.NDArray[np.float64]


class OutputFormat(Enum):  # pylint: disable=invalid-name
    eps = '.eps'
    jpeg = '.jpeg'
    jpg = '.jpg'
    pdf = '.pdf'
    pgf = '.pgf'
    png = '.png'
    ps = '.ps'
    raw = '.raw'
    rgba = '.rgba'
    svg = '.svg'
    svgz = '.svgz'
    tif = '.tif'
    tiff = '.tiff'


# --- Coverage colormap utilities ---
_cmap = mcolors.LinearSegmentedColormap.from_list("red_blue", ["red", "blue"])
_norm = mcolors.Normalize(vmin=0, vmax=100)

def coverage_to_color(value: float):
    """Map a coverage value (0–100) to a color."""
    return _cmap(_norm(value))


def _plot_line_segment(ax: plt.Axes,
                       line_segment: PolarLineSegment,
                       color: str = "grey",
                       linewidth: float = 0.5) -> None:
    """Draw a single line segment in polar coordinates."""
    thetas = [line_segment.start_coord.radians, line_segment.end_coord.radians]
    radii = [line_segment.start_coord.radius, line_segment.end_coord.radius]
    ax.plot(thetas, radii, color=color, linewidth=linewidth)


def plot(*,
         line_segments: Iterable[PolarLineSegment] | Iterable[list[tuple[PolarLineSegment, float]]],
         fig_size: float,
         line_spacing: float,
         line_width: float,
         circle_size: float,
         include_clipped_reads: bool,
         figure_file: Path,
         coverage: list[float] = None) -> None:
    """
    Plot reads in polar coordinates.
    If coverage is provided, each read is drawn as per-base colored segments.
    """
    with plt.style.context('ggplot'):
        fig = plt.figure(figsize=(fig_size,) * 2)
        ax = fig.add_subplot(111, polar=True)
        ax.grid(False)
        ax.set_rticks([])
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_theta_zero_location('N')
        ax.set_facecolor('white')
        ax.axis('off')

        # Draw baseline circle
        basis_radius = circle_size - 2 * line_spacing
        basis_curve = PolarLineSegment(
            PolarCoordinate(0, basis_radius),
            PolarCoordinate(2 * np.pi, basis_radius))

        _plot_line_segment(ax, basis_curve, color="black", linewidth=2)

        # Draw reads
        if coverage is None:
            # Grey mode
            for line_segment in line_segments:
                _plot_line_segment(ax, line_segment,
                                   color="grey",
                                   linewidth=line_width)
        else:
            # Coverage mode: each read is a list of (segment, cov_val)
            for read_segments in line_segments:
                for segment, cov_val in read_segments:
                    color = coverage_to_color(cov_val)
                    _plot_line_segment(ax, segment,
                                       color=color,
                                       linewidth=line_width)

    plt.savefig(figure_file, bbox_inches='tight')
