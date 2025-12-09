"""
Utility functions and classes.
"""
import functools
import math
import os
from contextlib import contextmanager
from pathlib import Path

import numpy as np

from concatmap.typing import Array1D


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
        self.rad_per_base = math.tau / reference_length

    def __call__(self, pos: int) -> float:
        return pos * self.rad_per_base


class CoverageInterpolator:

    def __init__(self, coverage: list[float], normalize: bool):
        self.coverage = self._normalize(coverage) if normalize else coverage

        conv = PositionToAngleConverter(len(coverage))
        self.angles = [conv(i) for i, _ in enumerate(coverage, 1)]

    def __call__(self, angles: Array1D) -> Array1D:
        return np.interp(angles % math.tau, self.angles, self.coverage)

    def _normalize(self, xs: list[float]) -> list[float]:
        min_, max_ = self._minmax(xs)
        range_ = max_ - min_
        return [(x - min_) / range_ for x in xs] if range_ else [1.] * len(xs)

    @staticmethod
    def _minmax[T](xs: list[T]) -> tuple[T, T]:
        def f(minmax, x):
            min_, max_ = minmax
            return min(min_, x), max(max_, x)
        return functools.reduce(f, xs, (xs[0],) * 2)
