from collections.abc import Iterable
from enum import Enum
from pathlib import Path

from matplotlib import pyplot as plt
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


def _plot_line_segment(
        ax: plt.Axes,
        line_segment: PolarLineSegment,
        n_points: int = 500,
        *args,
        **kwargs
) -> None:
    thetas = np.linspace(
        line_segment.start_coord.radians,
        line_segment.end_coord.radians,
        n_points)
    # NOTE: Originally radii was derived from `scipy.interp1d`, which has been
    #       superseded by `np.interp`. However, current use-case does not
    #       warrant it, as `np.linspace` should work.
    radii = np.linspace(
        line_segment.start_coord.radius,
        line_segment.end_coord.radius,
        n_points)
    ax.plot(thetas, radii, *args, **kwargs)


def plot(
        line_segments: Iterable[PolarLineSegment],
        fig_size: float,
        line_spacing: float,
        line_width: float,
        circle_size: float,
        clip: bool,
        figure_file: Path,
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
            PolarCoordinate(360, basis_radius))
        _plot_line_segment(ax, basis_curve, linewidth=5)

        # TODO: Draw clipped reads

        for line_segment in line_segments:
            _plot_line_segment(ax, line_segment, color='grey', linewidth=line_width)

    # TODO: Flip image so it is cw instead of ccw?
    plt.savefig(figure_file, bbox_inches='tight')
