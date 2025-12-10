import math
from dataclasses import dataclass
from typing import NamedTuple


@dataclass(slots=True, frozen=True)
class PolarCoordinate:
    """
    Store polar coordinate pairs.

    :ivar radians: The point's angular coordinate in radians.
    :ivar radius: The point's distance from the pole.
    """
    radians: float
    radius: float

    @property
    def degrees(self) -> float:
        """
        :return: Converted angular component in degrees.
        """
        return self.radians * 360 / math.tau

    def get_pair(self, radians: bool = False) -> tuple[float, float]:
        """
        Get the ordered pair of angle and radius.

        :param radians: Whether to report the angle in radians. Degrees is default.
        :return: An ordered pair of (angle, radius).
        """
        return self.radians if radians else self.degrees, self.radius


@dataclass(slots=True, frozen=True)
class PolarLineSegment:
    """
    Store a polar coordinate line segment.

    :ivar start_coord: The start coordinate.
    :ivar end_coord: The end coordinate.
    """
    start_coord: PolarCoordinate
    end_coord: PolarCoordinate

    @property
    def thetas(self) -> tuple[float, float]:
        return self.start_coord.radians, self.end_coord.radians

    @property
    def radii(self) -> tuple[float, float]:
        return self.start_coord.radius, self.end_coord.radius


class SamFileRead(NamedTuple):
    """
    A container to hold start and end positions from the samfile.

    Reference start and end follow Python convention of being zero-based and
    half-open.
    """
    reference_start: int
    reference_end: int
    clipped_start: int
    clipped_end: int
