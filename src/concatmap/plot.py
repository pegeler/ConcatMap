from collections.abc import Iterable
from enum import Enum
from pathlib import Path

from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

from concatmap.struct import PolarCoordinate
from concatmap.struct import PolarLineSegment


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


def _plot_line_segment(
        ax: plt.Axes,
        line_segment: PolarLineSegment,
        n_points: int = 500,
        coverage: list[float] = None,
        *args,
        **kwargs
) -> None:
    """
    Plot a line segment in polar coordinates.
    If coverage is provided, color each tiny sub-segment individually.
    """
    thetas = np.linspace(
        line_segment.start_coord.radians,
        line_segment.end_coord.radians,
        n_points)
    radii = np.linspace(
        line_segment.start_coord.radius,
        line_segment.end_coord.radius,
        n_points)

    if coverage is None:
        # Default: single grey line
        ax.plot(thetas, radii, *args, **kwargs)
    else:
        # Per-base coloring: draw each tiny segment separately
        for i in range(len(thetas) - 1):
            cov_val = coverage[i % len(coverage)]
            color = coverage_to_color(cov_val)
            ax.plot(
                [thetas[i], thetas[i+1]],
                [radii[i], radii[i+1]],
                color=color,
                linewidth=kwargs.get("linewidth", 0.5)
            )


def plot(
        *,
        line_segments: Iterable[PolarLineSegment] | Iterable[list[tuple[PolarLineSegment, float]]],
        fig_size: float,
        line_spacing: float,
        line_width: float,
        circle_size: float,
        clip: bool,
        figure_file: Path,
        coverage: list[float] = None,
) -> None:
    """
    Plot reads in polar coordinates.
    If coverage is provided, each read is drawn as per-base colored segments.
    """
    with plt.style.context('ggplot'):
        fig = plt.figure(figsize=(fig_size, ) * 2)
        ax = fig.add_subplot(111, polar=True)
        ax.grid(False)
        ax.set_rticks([])
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_theta_zero_location('N')
        ax.set_facecolor('white')
        ax.axis('off')

        basis_radius = circle_size - 2 * line_spacing
        basis_curve = PolarLineSegment(
            PolarCoordinate(0, basis_radius),
            PolarCoordinate(360, basis_radius))
        _plot_line_segment(ax, basis_curve, linewidth=5)

        # TODO: Draw clipped reads if needed

        if coverage is None:
            # Default: grey lines
            for line_segment in line_segments:
                _plot_line_segment(ax, line_segment,
                                   color='grey',
                                   linewidth=line_width)
        else:
            # Coverage-aware: each read is a list of (segment, cov_value)
            for read_segments in line_segments:
                for segment, cov_val in read_segments:
                    color = coverage_to_color(cov_val)
                    _plot_line_segment(ax, segment,
                                       color=color,
                                       linewidth=line_width)

    # TODO: Flip image so it is cw instead of ccw?
    plt.savefig(figure_file, bbox_inches='tight')
