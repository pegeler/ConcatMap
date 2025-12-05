import math
import os
from pathlib import Path
from concatmap.utils import PositionToAngleConverter, chdir_temporarily
import pytest


def test_position_to_angle_converter():
    conv = PositionToAngleConverter(reference_length=360)
    # Position 1 should map to angle 0
    assert conv(1) == pytest.approx(0.0)
    # Position 91 should map to ~pi/2
    assert conv(91) == pytest.approx(math.pi/2, rel=1e-6)

def test_chdir_temporarily(tmp_path):
    original = Path.cwd()
    with chdir_temporarily(tmp_path):
        assert Path.cwd() == tmp_path
    assert Path.cwd() == original
