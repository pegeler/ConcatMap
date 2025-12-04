import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

from concatmap.struct import PolarCoordinate, PolarLineSegment
from concatmap import plot

# Make a few fake line segments
segments = [
    PolarLineSegment(PolarCoordinate(0, 0.5), PolarCoordinate(np.pi/4, 0.5)),
    PolarLineSegment(PolarCoordinate(np.pi/4, 0.6), PolarCoordinate(np.pi/2, 0.6)),
    PolarLineSegment(PolarCoordinate(np.pi/2, 0.7), PolarCoordinate(3*np.pi/4, 0.7)),
]

# Fake coverage values for testing
colored_segments = [
    [(segments[0], 10.0), (segments[1], 50.0), (segments[2], 90.0)]
]

# Test grey mode
plot.plot(
    line_segments=segments,
    fig_size=5,
    line_spacing=0.05,
    line_width=2,
    circle_size=0.5,
    include_clipped_reads=False,
    figure_file=Path("test_grey.png"),
)

# Test coverage mode
plot.plot(
    line_segments=colored_segments,
    fig_size=5,
    line_spacing=0.05,
    line_width=2,
    circle_size=0.5,
    include_clipped_reads=False,
    figure_file=Path("test_color.png"),
    coverage=[0, 10, 50, 90]  # dummy coverage array
)

print("Generated test_grey.png and test_color.png")
