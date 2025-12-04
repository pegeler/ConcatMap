"""
Utility functions and classes.
"""
import math
import os
from contextlib import contextmanager
from pathlib import Path


@contextmanager
def chdir_temporarily(path: Path):
    cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(cwd)


class PositionToAngleConverter:
    """
    Callable class that will convert a position to the angle component of
    polar coordinates relative to the reference sequence length.
    """

    def __init__(self, reference_length: int) -> None:
        self.reference_length = reference_length
        self.rad_per_base = (2 * math.pi) / reference_length

    def __call__(self, pos: int) -> float:
        return (pos - 1) * self.rad_per_base
