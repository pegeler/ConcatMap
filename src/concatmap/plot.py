import abc
import math
from collections.abc import Iterator
from enum import Enum
from itertools import pairwise
from pathlib import Path
from typing import Callable

import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

from concatmap.struct import PlacedRead
from concatmap.struct import PolarCoordinate
from concatmap.struct import PolarLineSegment
from concatmap.struct import ReadSegmentType
from concatmap.struct import SamFileRead
from concatmap.typing import Array1D
from concatmap.utils import AngularCoordinatesInterpolator
from concatmap.utils import Debug
from concatmap.utils import PositionToAngleConverter
from concatmap.utils import normalize


class OutputFormat(Enum):
    # Members mirror file-extension names, so they are intentionally lowercase
    # rather than UPPER_CASE class constants.
    # pylint: disable=invalid-name
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


class AbstractPlotter(abc.ABC):

    BASIS_LINEWIDTH = 5
    BASIS_COLOR = 'red'
    CLIPPED_COLOR = 'red'
    # Reads are drawn as a dense radial grating (one arc per read). At the
    # default ~100 dpi the line density outruns the pixel density and the arcs
    # alias into a moire; a higher save resolution suppresses it.
    DPI = 300

    def __init__(
            self,
            *,
            reads: list[SamFileRead],
            reference_length: int,
            fig_size: float,
            line_spacing: float,
            line_width: float,
            circle_size: float,
            include_clipped_reads: bool,
            figure_file: Path,
    ) -> None:
        self.reads = reads
        self.reference_length = reference_length
        self.conv = PositionToAngleConverter(reference_length)
        self.fig_size = fig_size
        self.line_spacing = line_spacing
        self.line_width = line_width
        self.circle_size = circle_size
        self.include_clipped_reads = include_clipped_reads
        self.figure_file = figure_file

    @abc.abstractmethod
    def _drawLineSegment(
            self,
            ax: plt.Axes,
            thetas: Array1D,
            radii: Array1D,
            read: SamFileRead,
    ) -> None:
        ...

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
        ax.plot(thetas, radii, color=self.BASIS_COLOR, linewidth=self.BASIS_LINEWIDTH)

    def _drawClippedReads(self, ax: plt.Axes) -> None:
        placed_reads = self._convertReadsToLineSegments(
            self.reads,
            self.line_spacing,
            self.circle_size,
            ReadSegmentType.CLIPPED,
        )
        for placed in placed_reads:
            thetas, radii = self._linearize(placed.curve)
            # A clip extension longer than the reference sweeps past a full
            # turn and self-overlaps into a misleading ring; skip it.
            if abs(thetas[-1] - thetas[0]) > math.tau:
                continue
            ax.plot(thetas, radii, color=self.CLIPPED_COLOR, linewidth=self.line_width)

    def _drawReads(self, ax: plt.Axes) -> None:
        placed_reads = self._convertReadsToLineSegments(
            self.reads,
            self.line_spacing,
            self.circle_size,
        )
        for placed in placed_reads:
            thetas, radii = self._linearize(placed.curve)
            self._drawLineSegment(ax, thetas, radii, placed.read)

    def _saveFigure(self) -> None:
        # TODO: Flip image so it is cw instead of ccw?
        plt.savefig(self.figure_file, bbox_inches='tight', dpi=self.DPI)

    def _convertReadsToLineSegments(
            self,
            reads: list[SamFileRead],
            line_spacing: float,
            basis_radius: float,
            segment_type: ReadSegmentType = ReadSegmentType.MAPPED,
    ) -> Iterator[PlacedRead]:
        for i, read in enumerate(reads, 1):
            radius = basis_radius + line_spacing * i
            for start, end in read.getSegments(segment_type):
                yield PlacedRead(read, PolarLineSegment(
                    PolarCoordinate(self.conv(start), radius),
                    PolarCoordinate(self.conv(end), radius),
                ))

    @staticmethod
    def _linearize(
            line_segment: PolarLineSegment,
            n_points: int = 200,
    ) -> tuple[Array1D, Array1D]:
        if Debug().is_debug:
            n_points = 100
        thetas = np.linspace(*line_segment.thetas, n_points)
        radii = np.linspace(*line_segment.radii, n_points)
        return thetas, radii


class DefaultPlotter(AbstractPlotter):

    LINE_COLOR = 'grey'

    def _drawLineSegment(
            self,
            ax: plt.Axes,
            thetas: Array1D,
            radii: Array1D,
            read: SamFileRead,
    ) -> None:
        ax.plot(thetas, radii, color=self.LINE_COLOR, linewidth=self.line_width)


class MulticolorLinePlotter(AbstractPlotter):

    BASIS_COLOR = 'black'
    CLIPPED_COLOR = 'grey'

    def __init__(
            self,
            values: Callable[[Array1D], Array1D] | list[float],
            **kwargs,
    ) -> None:
        """
        :param values: A callable that will provide a normalized color map
                value (on the interval [0, 1]), given an angle in radians.
                Possibly and instance of ``AngularCoordinatesInterpolator``.
                Or a list of values for each reference sequence position to be
                normalized and projected on the color map.
        """
        self.interpolator = (
            values
            if callable(values)
            else AngularCoordinatesInterpolator(normalize(values))
        )
        super().__init__(**kwargs)

    def _drawLineSegment(
            self,
            ax: plt.Axes,
            thetas: Array1D,
            radii: Array1D,
            read: SamFileRead,
    ) -> None:
        lines = list(pairwise(zip(thetas, radii)))
        midpoints = np.array([(a + b) / 2 for (a, _), (b, _) in lines])
        segments = LineCollection(
            lines,
            linewidths=self.line_width,
            colors=cm.plasma(self.interpolator(midpoints)),
        )
        ax.add_collection(segments)
        ax.set_rmax(radii[0])  # TODO: just set this once at the end


IGV_BASE_COLORS: dict[str, str] = {
    'A': '#00C800',  # green
    'C': '#0000C8',  # blue
    'G': '#D17105',  # brown/orange
    'T': '#FF0000',  # red
}


class MismatchPlotter(AbstractPlotter):
    """
    Color each read gray where it matches the reference and with IGV nucleotide
    colors at substituted bases (IGV/MSA alignment-track style).

    A read is drawn as a single gray arc, then each substitution is overpainted
    as a short arc spanning exactly that base's angular slot. Work is
    proportional to the number of substitutions, so clean reads cost the same as
    ``DefaultPlotter``.
    """

    BASIS_COLOR = 'black'        # not red: red is the T substitution color
    CLIPPED_COLOR = 'lightgrey'  # neutral; none of the four base colors
    LINE_COLOR = 'grey'
    BASE_COLORS = IGV_BASE_COLORS

    # One base subtends a negligible angle, so a straight two-point chord is
    # visually indistinguishable from the arc; no dense linearization needed.
    _MISMATCH_ARC_POINTS = 2

    def _drawLineSegment(
            self,
            ax: plt.Axes,
            thetas: Array1D,
            radii: Array1D,
            read: SamFileRead,
    ) -> None:
        ax.plot(thetas, radii, color=self.LINE_COLOR, linewidth=self.line_width)
        radius = radii[0]  # constant along a mapped read segment
        for mismatch in read.mismatches:
            color = self.BASE_COLORS.get(mismatch.read_base)
            if color is None:
                continue  # ambiguous base (e.g. N): leave the gray body showing
            base_arc = PolarLineSegment(
                PolarCoordinate(self.conv(mismatch.position), radius),
                PolarCoordinate(self.conv(mismatch.position + 1), radius),
            )
            arc_thetas, arc_radii = self._linearize(base_arc, self._MISMATCH_ARC_POINTS)
            ax.plot(arc_thetas, arc_radii, color=color, linewidth=self.line_width)
