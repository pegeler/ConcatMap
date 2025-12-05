import math
from collections.abc import Iterable
from enum import Enum
from itertools import pairwise
from pathlib import Path
from typing import Callable

from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
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


def _plot_line_segment(
        ax: plt.Axes,
        line_segment: PolarLineSegment,
        n_points: int = 200,
        coverage_interpolator: Callable[[Array1D], Array1D] | None = None,
        *args,
        **kwargs
) -> None:
    thetas = np.linspace(*line_segment.thetas, n_points)
    radii = np.linspace(*line_segment.radii, n_points)

    if coverage_interpolator is not None:
        lines = list(pairwise(zip(thetas, radii)))
        midpoints = np.array([(a + b) / 2 for (a, _), (b, _) in lines])
        segments = LineCollection(
            lines,
            linewidths=kwargs['linewidth'],
            colors=cm.plasma(coverage_interpolator(midpoints)),
        )
        ax.add_collection(segments)
        ax.set_rmax(radii[0])
    else:
        ax.plot(thetas, radii, *args, **kwargs)


class _CoverageInterpolator:

    def __init__(self, coverage: list[float]):
        self.coverage = coverage

        conv = PositionToAngleConverter(len(coverage))
        self.angles = [conv(i) for i in range(len(coverage))]

    def __call__(self, angles: Array1D) -> Array1D:
        return np.interp(angles % math.tau, self.angles, self.coverage)


def plot(
        *,
        line_segments: Iterable[PolarLineSegment],
        fig_size: float,
        line_spacing: float,
        line_width: float,
        circle_size: float,
        include_clipped_reads: bool,
        figure_file: Path,
        coverage: list[float] | None = None,
) -> None:
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
            PolarCoordinate(math.tau, basis_radius))
        _plot_line_segment(ax, basis_curve, color='black', linewidth=5)

        # TODO: Draw clipped reads

        cov_interp = None
        if coverage is not None:
            cov_interp = _CoverageInterpolator(coverage)

        for line_segment in line_segments:
            _plot_line_segment(
                ax,
                line_segment,
                color='grey',
                coverage_interpolator=cov_interp,
                linewidth=line_width
            )

    # TODO: Flip image so it is cw instead of ccw?
    plt.savefig(figure_file, bbox_inches='tight')
