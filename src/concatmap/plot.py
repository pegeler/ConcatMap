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
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

from concatmap.struct import PlacedRead
from concatmap.struct import PolarCoordinate
from concatmap.struct import PolarLineSegment
from concatmap.struct import ReadSegmentType
from concatmap.struct import SamFileRead
from concatmap.typing import Array1D
from concatmap.utils import AngularCoordinatesInterpolator
from concatmap.utils import Debug
from concatmap.utils import PositionToAngleConverter
from concatmap.utils import minmax
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

    # Legend and tick-label text is set in points, which don't scale with
    # fig_size (inches) the way the plot geometry does; at large fig_size the
    # default matplotlib font shrinks to illegible relative to the figure. Scale
    # it off the fig_size at which the default matplotlib font size (10pt)
    # looks right.
    _REFERENCE_FIG_SIZE = 10.0
    _BASE_FONTSIZE = 10.0

    # Position ticks ring the outside of the read stack. Lengths are fractions
    # of the stack's outer radius, because the stack grows with read count and
    # a fixed radial length would vanish on deep datasets.
    _TICK_TARGET_COUNT = 12
    _TICK_LENGTH_FRACTION = 0.015
    _TICK_LABEL_PAD_FRACTION = 0.01
    _TICK_LINEWIDTH = 1.0

    @property
    def _text_scale(self) -> float:
        return self.fig_size / self._REFERENCE_FIG_SIZE

    @property
    def _text_font_size(self) -> float:
        return self._BASE_FONTSIZE * self._text_scale

    @property
    def _outer_radius(self) -> float:
        # One spacing past the last read, matching _convertReadsToLineSegments.
        return self.circle_size + self.line_spacing * (len(self.reads) + 1)

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
            legend: bool = False,
            ticks: bool = False,
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
        self.legend = legend
        self.ticks = ticks

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
            if self.ticks:
                self._drawTicks(ax)
            if self.legend:
                self._drawLegend(ax)
            self._saveFigure()

    def _setup(self) -> plt.Axes:
        fig = plt.figure(figsize=(self.fig_size,) * 2)
        ax = fig.add_subplot(111, polar=True)
        ax.grid(False)
        ax.set_rticks([])
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
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

    def _drawTicks(self, ax: plt.Axes) -> None:
        # Geometry (conv, radii) is shared by every plotter, so ticks live in
        # the base class rather than a per-plotter hook. Position 0 sits at 12
        # o'clock, matching the origin of the angle converter.
        tick_length = self._TICK_LENGTH_FRACTION * self._outer_radius
        label_pad = self._TICK_LABEL_PAD_FRACTION * self._outer_radius
        inner = self._outer_radius
        outer = inner + tick_length
        for position in self._tickPositions():
            theta = self.conv(position)
            # Ticks fall outside the polar axes' data limits: don't clip them,
            # and don't let them rescale the radial axis (which would shrink
            # the reads relative to a plot without ticks).
            ax.plot(
                [theta, theta],
                [inner, outer],
                color=self.BASIS_COLOR,
                linewidth=self._TICK_LINEWIDTH,
                clip_on=False,
                scalex=False,
                scaley=False,
            )
            ha, va = self._tickLabelAlignment(theta)
            ax.text(
                theta,
                outer + label_pad,
                f'{position:,}',
                ha=ha,
                va=va,
                fontsize=self._text_font_size,
            )

    def _tickPositions(self) -> list[int]:
        # MaxNLocator picks round 1/2/2.5/5 x 10^k steps, so labels read as
        # 0, 500, 1000, ... whatever the reference length. The end of the
        # range coincides with position 0 on a circle, so it is dropped.
        locator = MaxNLocator(
            nbins=self._TICK_TARGET_COUNT,
            steps=[1, 2, 2.5, 5, 10],
            integer=True,
        )
        return [
            int(position)
            for position in locator.tick_values(0, self.reference_length)
            if 0 <= position < self.reference_length
        ]

    @staticmethod
    def _tickLabelAlignment(theta: float) -> tuple[str, str]:
        # Anchor each label on the side facing away from the circle so text
        # grows outward instead of back over the tick. theta runs clockwise
        # from 12 o'clock, so x = sin(theta) and y = cos(theta).
        x, y = math.sin(theta), math.cos(theta)
        ha = 'left' if x > 0.3 else 'right' if x < -0.3 else 'center'
        va = 'bottom' if y > 0.3 else 'top' if y < -0.3 else 'center'
        return ha, va

    def _saveFigure(self) -> None:
        plt.savefig(self.figure_file, bbox_inches='tight', dpi=self.DPI)

    def _drawLegend(self, ax: plt.Axes) -> None:
        """Hook: draw a legend/key. No-op unless a subclass provides one."""

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
        if callable(values):
            self.interpolator = values
            self._value_range = None            # no absolute domain to label
        else:
            self._value_range = minmax(values)  # (vmin, vmax) in real depth units
            self.interpolator = AngularCoordinatesInterpolator(normalize(values))
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

    def _drawLegend(self, ax: plt.Axes) -> None:
        # Rebuild the depth -> color mapping the line segments use
        # (plasma((d - vmin) / (vmax - vmin))) as a standalone colorbar. The bar
        # is scaled to this plot's own min/max depth, so its ticks are absolute
        # coverage values but the scale is relative between plots.
        if self._value_range is None:
            return  # callable values: no absolute depth scale to label
        vmin, vmax = self._value_range
        mappable = cm.ScalarMappable(norm=Normalize(vmin, vmax), cmap=cm.plasma)
        cbar = ax.figure.colorbar(
            mappable,
            ax=ax,
            fraction=0.046,
            pad=0.04,
            shrink=0.6,
        )
        cbar.set_label(
            'Read depth',
            rotation=270,
            labelpad=15 * self._text_scale,
            fontsize=self._text_font_size,
        )
        cbar.ax.tick_params(labelsize=self._text_font_size)


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

    def _drawLegend(self, ax: plt.Axes) -> None:
        # A discrete key: one line swatch per substituted base in its IGV color,
        # plus the gray used where the read matches the reference. Line2D proxies
        # mirror how reads are drawn (colored line segments).
        handles = [
            Line2D([], [], color=color, linewidth=3, label=base)
            for base, color in self.BASE_COLORS.items()
        ]
        handles.append(
            Line2D([], [], color=self.LINE_COLOR, linewidth=3, label='match'))
        ax.legend(
            handles=handles,
            title='Base',
            loc='center left',
            bbox_to_anchor=(1.0, 0.5),
            frameon=False,
            fontsize=self._text_font_size,
            title_fontsize=self._text_font_size,
        )
