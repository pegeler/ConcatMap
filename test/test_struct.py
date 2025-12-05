import math
import pytest
from concatmap.struct import PolarCoordinate, PolarLineSegment, SamFileRead

def test_polarcoordinate_degrees_and_get_pair():
    coord = PolarCoordinate(radians=math.pi, radius=1.0)
    assert pytest.approx(coord.degrees, rel=1e-6) == 180.0
    assert coord.get_pair(radians=True) == (math.pi, 1.0)
    assert coord.get_pair(radians=False)[0] == pytest.approx(180.0)

def test_polarline_segment_properties():
    start = PolarCoordinate(0, 1.0)
    end = PolarCoordinate(math.pi/2, 2.0)
    seg = PolarLineSegment(start, end)
    assert seg.thetas == (0, math.pi/2)
    assert seg.radii == (1.0, 2.0)

def test_samfileread_namedtuple():
    read = SamFileRead(10, 50, 12, 48)
    assert read.reference_start == 10
    assert read.reference_end == 50
    assert read.clipped_start == 12
    assert read.clipped_end == 48
