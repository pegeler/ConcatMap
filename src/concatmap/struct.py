import math
from dataclasses import dataclass
from typing import NamedTuple


@dataclass(slots=True, frozen=True)
class PolarCoordinate:
    """
    Store polar coordinate pairs.

    :ivar angle: The point's angular coordinate in degrees.
    :ivar radius: The point's distance from the pole.
    """
    angle: float
    radius: float

    @property
    def radians(self) -> float:
        """
        :return: Converted angular component in radians.
        """
        return self.angle * 2 * math.pi / 360

    def get_pair(self, radians: bool = False) -> tuple[float, float]:
        """
        Get the ordered pair of angle and radius.

        :param radians: Whether to report the angle in radians. Degrees is default.
        :return: An ordered pair of (angle, radius).
        """
        return self.radians if radians else self.angle, self.radius


@dataclass(slots=True, frozen=True)
class PolarLineSegment:
    """
    Store a polar coordinate line segment.

    :ivar start_coord: The start coordinate.
    :ivar end_coord: The end coordinate.
    """
    start_coord: PolarCoordinate
    end_coord: PolarCoordinate


class SamFileRead(NamedTuple):
    """
    A container to hold start and end positions from the samfile.
    """
    reference_start: int
    reference_end: int
    clipped_start: int
    clipped_end: int
