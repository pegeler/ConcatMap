"""
Utility functions and classes.
"""
import functools
import math
import os
from contextlib import contextmanager
from pathlib import Path
from typing import Self

import numpy as np

from concatmap.typing import Array1D


class Debug:
    """
    Singleton indicating whether this is a debug run.
    """
    _instance: Self | None = None
    is_debug: bool = False

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance


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

    The result is intentionally not wrapped to ``[0, tau)``: positions in the
    second concatenated copy, and negative clip projections that fall before the
    origin, return angles outside that range. Matplotlib's polar axes wrap them
    at render time, which is what lets origin-spanning reads draw as continuous
    arcs. Do not add a modulo here --- the clip-span guard in the plotter relies
    on the raw, unwrapped span.
    """

    def __init__(self, reference_length: int) -> None:
        self.reference_length = reference_length
        self.rad_per_base = math.tau / reference_length

    def __call__(self, pos: int) -> float:
        return pos * self.rad_per_base


class AngularCoordinatesInterpolator:
    """
    Callable class that will take a list of values and map each value to an
    angle based on its position within the list. Given an array of new angles,
    it will interpolate the value at each angle.
    """

    def __init__(self, values: list[float]):
        self.values = values[:]
        conv = PositionToAngleConverter(len(values))
        self.angles = [conv(i) for i, _ in enumerate(values)]
        # Close the loop with a sample at tau equal to position 0, so query
        # angles past the last position interpolate across the origin instead of
        # clamping to the last value (which left a one-bin color seam at North).
        self.angles.append(math.tau)
        self.values.append(self.values[0])

    def __call__(self, angles: Array1D) -> Array1D:
        return np.interp(angles % math.tau, self.angles, self.values)


def normalize(xs: list[float]) -> list[float]:
    """
    Scale the values of a list onto an interval [0, 1].

    :param xs: The list to be normalized.
    :return: A list the same length as the input.
    """
    min_, max_ = minmax(xs)
    range_ = max_ - min_
    return [(x - min_) / range_ for x in xs] if range_ else [1.] * len(xs)


def minmax[T](xs: list[T]) -> tuple[T, T]:
    """
    Get the min and max of a list in one pass through the data.

    :param xs: The list to be inspected.
    :return: A tuple containing the minimum and maximum values in the list.
    """
    def f(mm, x):
        min_, max_ = mm
        return min(min_, x), max(max_, x)
    return functools.reduce(f, xs, (xs[0],) * 2)
