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

type Array1D = npt.NDArray[tuple[int], np.float64]


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
        coverage: list[int] | None = None,
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

    if coverage is not None:
        # FIXME: This doesn't work.
        coverage_interpolator = _CoverageInterpolator(coverage)
        color = cm.plasma(coverage_interpolator(thetas))
        for i, (t, r) in enumerate(zip(thetas, radii)):
            kwargs['color'] = color[i, :]
            ax.plot(t, r, *args, **kwargs)
            return

    ax.plot(thetas, radii, *args, **kwargs)


class _CoverageInterpolator:

    def __init__(self, coverage: list[float]):
        self.coverage = coverage

        conv = PositionToAngleConverter(len(coverage))
        self.angles = [conv(i) for i in range(len(coverage))]

    def __call__(self, angles: Array1D) -> Array1D:
        return np.interp(angles % (2 * math.pi), self.angles, self.coverage)


def plot(
        *,
        line_segments: Iterable[PolarLineSegment] | Iterable[list[tuple[PolarLineSegment, float]]],
        fig_size: float,
        line_spacing: float,
        line_width: float,
        circle_size: float,
        include_clipped_reads: bool,
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

        for line_segment in line_segments:
            _plot_line_segment(
                ax,
                line_segment,
                color='grey',
                coverage=coverage,
                linewidth=line_width
            )

    # TODO: Flip image so it is cw instead of ccw?
    plt.savefig(figure_file, bbox_inches='tight')
