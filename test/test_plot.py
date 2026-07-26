"""Unit tests for the pure coloring decisions in ``plot``."""
from concatmap.plot import IGV_BASE_COLORS
from concatmap.plot import MismatchPlotter


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
