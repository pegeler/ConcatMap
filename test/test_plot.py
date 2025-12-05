import math
from pathlib import Path
import matplotlib.pyplot as plt
from concatmap.struct import PolarCoordinate, PolarLineSegment
from concatmap.plot import plot, _plot_line_segment, _CoverageInterpolator, OutputFormat
import numpy as np
import pytest


def test_plot_creates_file(tmp_path):
    # Create a simple line segment
    seg = PolarLineSegment(
        PolarCoordinate(0, 0.5),
        PolarCoordinate(math.pi/4, 0.5)
    )
    outfile = tmp_path / "test_plot.png"
    plot(
        line_segments=[seg],
        fig_size=4,
        line_spacing=0.05,
        line_width=1.0,
        circle_size=0.5,
        include_clipped_reads=False,
        figure_file=outfile,
        coverage=None
    )
    assert outfile.exists()
    assert outfile.stat().st_size > 0  # file is not empty

def test_coverage_interpolator():
    coverage = [0, 10, 20, 30]
    interp = _CoverageInterpolator(coverage)
    #Note: need to pass in a numpy array to use the _CoverageInterpolate() function
    angles = np.array([0.0, math.pi/2, math.pi, 3*math.pi/2])
    values = interp(angles)
    assert len(values) == 4
    assert min(values) >= 0
    assert max(values) <= 30

def test_outputformat_enum_values():
    # check all expected output formats
    formats = [f.value for f in OutputFormat]
    assert '.png' in formats
    assert '.pdf' in formats
    assert '.svg' in formats


def test_plot_line_segment_without_coverage(tmp_path):
    fig, ax = plt.subplots(subplot_kw={'polar': True})
    seg = PolarLineSegment(
        PolarCoordinate(0, 1.0),
        PolarCoordinate(math.pi/2, 1.0)
    )
    _plot_line_segment(ax, seg, color='red', linewidth=2)
    # Save to file to ensure something was drawn
    outfile = tmp_path / "segment.png"
    fig.savefig(outfile)
    assert outfile.exists()
    assert outfile.stat().st_size > 0


def test_plot_line_segment_with_coverage(tmp_path):
    fig, ax = plt.subplots(subplot_kw={'polar': True})
    seg = PolarLineSegment(
        PolarCoordinate(0, 1.0),
        PolarCoordinate(math.pi/2, 1.0)
    )
    coverage = [0, 10, 20, 30]
    interp = _CoverageInterpolator(coverage)
    _plot_line_segment(ax, seg, color='blue', linewidth=2, coverage_interpolator=interp)
    outfile = tmp_path / "segment_cov.png"
    fig.savefig(outfile)
    assert outfile.exists()
    assert outfile.stat().st_size > 0


def test_plot_creates_file_with_and_without_coverage(tmp_path):
    seg = PolarLineSegment(
        PolarCoordinate(0, 0.5),
        PolarCoordinate(math.pi/4, 0.5)
    )
    outfile1 = tmp_path / "plot_grey.png"
    plot(
        line_segments=[seg],
        fig_size=4,
        line_spacing=0.05,
        line_width=1.0,
        circle_size=0.5,
        include_clipped_reads=False,
        figure_file=outfile1,
        coverage=None
    )
    assert outfile1.exists()
    assert outfile1.stat().st_size > 0

    outfile2 = tmp_path / "plot_cov.png"
    coverage = [0, 5, 10, 15, 20]
    plot(
        line_segments=[seg],
        fig_size=4,
        line_spacing=0.05,
        line_width=1.0,
        circle_size=0.5,
        include_clipped_reads=False,
        figure_file=outfile2,
        coverage=coverage
    )
    assert outfile2.exists()
    assert outfile2.stat().st_size > 0


@pytest.mark.parametrize("angles", [
    [0.0, math.pi/2, math.pi, 3*math.pi/2],  # list input
    np.array([0.0, math.pi/2, math.pi, 3*math.pi/2])  # numpy array input
])
def test_coverage_interpolator_accepts_list_and_array(angles):
    coverage = [0, 10, 20, 30]
    interp = _CoverageInterpolator(coverage)
    values = interp(angles)
    assert len(values) == 4
    assert min(values) >= 0
    assert max(values) <= 30