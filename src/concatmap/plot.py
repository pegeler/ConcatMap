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

from concatmap.struct import PolarCoordinate
from concatmap.struct import PolarLineSegment
from concatmap.typing import Array1D


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


class Plotter:

    def __init__(
            self,
            *,
            line_segments: Iterable[PolarLineSegment],
            fig_size: float,
            line_spacing: float,
            line_width: float,
            circle_size: float,
            include_clipped_reads: bool,
            figure_file: Path,
            coverage_interpolator: Callable[[Array1D], Array1D] | None = None,
    ) -> None:
        self.line_segments = line_segments
        self.fig_size = fig_size
        self.line_spacing = line_spacing
        self.line_width = line_width
        self.circle_size = circle_size
        self.include_clipped_reads = include_clipped_reads
        self.figure_file = figure_file
        self.coverage_interpolator = coverage_interpolator

    def plot(self) -> None:
        with plt.style.context('ggplot'):
            ax = self._setup()
            self._drawBasisCircle(ax)
            if self.include_clipped_reads:
                self._drawClippedReads(ax)
            self._drawReads(ax)
            self._saveFigure()

    def _setup(self) -> plt.Axes:
        fig = plt.figure(figsize=(self.fig_size,) * 2)
        ax = fig.add_subplot(111, polar=True)
        ax.grid(False)
        ax.set_rticks([])
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_theta_zero_location('N')
        ax.set_facecolor('white')
        ax.axis('off')
        return ax

    def _drawBasisCircle(self, ax: plt.Axes) -> None:
        basis_radius = self.circle_size - 2 * self.line_spacing
        basis_curve = PolarLineSegment(
            PolarCoordinate(0, basis_radius),
            PolarCoordinate(math.tau, basis_radius))
        thetas, radii = self._linearize(basis_curve)
        ax.plot(thetas, radii, color='black', linewidth=5)

    def _drawClippedReads(self, ax: plt.Axes) -> None:
        ...  # TODO

    def _drawReads(self, ax: plt.Axes) -> None:
        for line_segment in self.line_segments:
            thetas, radii = self._linearize(line_segment)
            if self.coverage_interpolator is not None:
                self._plotCoverageLineSegment(
                    ax,
                    thetas,
                    radii,
                    linewidth=self.line_width
                )
            else:
                ax.plot(thetas, radii, color='grey', linewidth=self.line_width)

    def _linearize(
            self,
            line_segment: PolarLineSegment,
            n_points: int = 200,
    ) -> tuple[Array1D, Array1D]:
        thetas = np.linspace(*line_segment.thetas, n_points)
        radii = np.linspace(*line_segment.radii, n_points)
        return thetas, radii

    def _plotCoverageLineSegment(
            self,
            ax: plt.Axes,
            thetas: Array1D,
            radii: Array1D,
            **kwargs
    ) -> None:
        lines = list(pairwise(zip(thetas, radii)))
        midpoints = np.array([(a + b) / 2 for (a, _), (b, _) in lines])
        segments = LineCollection(
            lines,
            linewidths=kwargs['linewidth'],
            colors=cm.plasma(self.coverage_interpolator(midpoints)),
        )
        ax.add_collection(segments)
        ax.set_rmax(radii[0])  # TODO: just set this once at the end

    def _saveFigure(self) -> None:
        # TODO: Flip image so it is cw instead of ccw?
        plt.savefig(self.figure_file, bbox_inches='tight')
