"""Unit tests for the pure coloring decisions in ``plot``."""
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # headless: no display needed for these tests
from matplotlib import pyplot as plt

from concatmap.plot import DefaultPlotter
from concatmap.plot import IGV_BASE_COLORS
from concatmap.plot import MismatchPlotter
from concatmap.plot import MulticolorLinePlotter


_PLOTTER_KWARGS = dict(
    reads=[],
    reference_length=100,
    fig_size=4.0,
    line_spacing=0.02,
    line_width=0.75,
    circle_size=0.45,
    include_clipped_reads=False,
    figure_file=Path('unused.png'),  # never written: we don't call plot()
)


def test_igv_colors_are_canonical():
    assert IGV_BASE_COLORS == {
        'A': '#00C800',
        'C': '#0000C8',
        'G': '#D17105',
        'T': '#FF0000',
    }


def test_known_bases_map_to_their_igv_color():
    colors = MismatchPlotter.BASE_COLORS
    for base, expected in IGV_BASE_COLORS.items():
        assert colors.get(base) == expected


def test_ambiguous_base_has_no_color():
    # A None lookup is the signal to leave the gray read body showing.
    assert MismatchPlotter.BASE_COLORS.get('N') is None


def test_mismatch_plotter_avoids_red_for_basis_and_clips():
    # Red is the T substitution color, so the basis circle and clip extensions
    # must not also be red (see design critique).
    assert MismatchPlotter.BASIS_COLOR != 'red'
    assert MismatchPlotter.CLIPPED_COLOR != 'red'
    assert MismatchPlotter.CLIPPED_COLOR not in IGV_BASE_COLORS.values()


def test_list_values_retain_real_depth_range():
    # The colorbar labels itself with the raw min/max depth, so the range must
    # survive normalization (which otherwise discards it).
    plotter = MulticolorLinePlotter(values=[3, 1, 4, 1, 5, 9, 2, 6],
                                    **_PLOTTER_KWARGS)
    assert plotter._value_range == (1, 9)


def test_callable_values_have_no_depth_range():
    # A pre-built interpolator carries no absolute domain to label.
    plotter = MulticolorLinePlotter(values=lambda angles: angles,
                                    **_PLOTTER_KWARGS)
    assert plotter._value_range is None


def test_depth_legend_adds_a_colorbar_axes():
    plotter = MulticolorLinePlotter(values=[1, 2, 3, 4], **_PLOTTER_KWARGS)
    fig = plt.figure()
    ax = fig.add_subplot(111, polar=True)
    plotter._drawLegend(ax)
    assert len(fig.axes) == 2  # polar axes + colorbar
    plt.close(fig)


def test_depth_legend_is_noop_without_a_depth_range():
    plotter = MulticolorLinePlotter(values=lambda angles: angles,
                                    **_PLOTTER_KWARGS)
    fig = plt.figure()
    ax = fig.add_subplot(111, polar=True)
    plotter._drawLegend(ax)
    assert len(fig.axes) == 1  # nothing added
    plt.close(fig)


def test_by_base_legend_draws_a_base_color_key():
    plotter = MismatchPlotter(**_PLOTTER_KWARGS)
    fig = plt.figure()
    ax = fig.add_subplot(111, polar=True)
    plotter._drawLegend(ax)
    legend = ax.get_legend()
    assert legend is not None
    labels = [t.get_text() for t in legend.get_texts()]
    assert labels == ['A', 'C', 'G', 'T', 'match']
    plt.close(fig)


def _plotter_for_length(reference_length: int) -> DefaultPlotter:
    return DefaultPlotter(**{
        **_PLOTTER_KWARGS,
        'reference_length': reference_length,
    })


def test_tick_positions_are_round_and_exclude_the_wraparound():
    positions = _plotter_for_length(5000)._tickPositions()
    assert positions[0] == 0
    assert all(p < 5000 for p in positions)
    steps = {b - a for a, b in zip(positions, positions[1:])}
    assert steps == {500}


def test_tick_positions_adapt_to_reference_length():
    positions = _plotter_for_length(16569)._tickPositions()  # human mtDNA
    assert positions[:3] == [0, 2000, 4000]
    assert len(positions) <= DefaultPlotter._TICK_TARGET_COUNT + 1


def test_tick_label_alignment_faces_away_from_the_circle():
    align = DefaultPlotter._tickLabelAlignment
    assert align(0.0) == ('center', 'bottom')         # 12 o'clock
    assert align(math.pi / 2) == ('left', 'center')   # 3 o'clock (clockwise)
    assert align(math.pi) == ('center', 'top')        # 6 o'clock
    assert align(3 * math.pi / 2) == ('right', 'center')


def test_draw_ticks_adds_a_tick_and_label_per_position_without_rescaling():
    plotter = _plotter_for_length(1000)
    fig = plt.figure()
    ax = fig.add_subplot(111, polar=True)
    ax.set_rmax(1.0)
    plotter._drawTicks(ax)
    n = len(plotter._tickPositions())
    assert len(ax.lines) == n
    assert [t.get_text() for t in ax.texts][:3] == ['0', '100', '200']
    assert ax.get_rmax() == 1.0
    plt.close(fig)


def test_ticks_do_not_rescale_an_autoscaled_radial_axis():
    # Default and by-base plots leave rmax to lazy autoscale; ticks (which sit
    # past the reads) must not be folded into it.
    plotter = _plotter_for_length(1000)
    rmax = []
    for with_ticks in (False, True):
        fig = plt.figure()
        ax = fig.add_subplot(111, polar=True)
        ax.plot([0, 1], [0.1, 0.3])
        if with_ticks:
            plotter._drawTicks(ax)
        rmax.append(ax.get_rmax())
        plt.close(fig)
    assert rmax[0] == rmax[1]


def test_legend_clears_tick_labels():
    plotter = MismatchPlotter(**{**_PLOTTER_KWARGS, 'reference_length': 1000})
    fig = plt.figure()
    ax = fig.add_subplot(111, polar=True)
    ax.set_rmax(plotter.circle_size)  # as with reads: the stack fills the axes
    assert plotter._contentRight(ax) == 1.0  # no labels: legend stays put
    plotter._drawTicks(ax)
    right = plotter._contentRight(ax)
    assert right > 1.0
    plotter._drawLegend(ax)
    anchor_x = ax.get_legend().get_bbox_to_anchor().transformed(
        ax.transAxes.inverted()).x0
    assert anchor_x == right
    plt.close(fig)


def test_ticks_are_off_by_default():
    assert DefaultPlotter(**_PLOTTER_KWARGS).ticks is False


def test_default_plotter_legend_hook_is_noop():
    # The base hook draws nothing, so --legend is harmless in non-depth modes.
    plotter = DefaultPlotter(**_PLOTTER_KWARGS)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    assert plotter._drawLegend(ax) is None
    assert len(fig.axes) == 1
    plt.close(fig)
