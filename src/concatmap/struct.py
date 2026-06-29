import math
from collections.abc import Iterator
from dataclasses import dataclass
from enum import auto
from enum import StrEnum
from typing import NamedTuple


type ReadSegment = tuple[int, int]
"""An inclusive (start, end) pair of reference positions."""


class ReadSegmentType(StrEnum):
    MAPPED = auto()
    CLIPPED = auto()


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

    def getMappedSegment(self) -> ReadSegment:
        """
        The aligned portion of the read as one inclusive segment.

        :return: The inclusive (start, end) of the mapped reference span.
        """
        return self.reference_start, self.reference_end - 1

    def getClippedSegments(self) -> list[ReadSegment]:
        """
        The non-overlapping clip extensions as inclusive segments.

        :return: 0–2 segments: leading extension (if any), then trailing
            (if any).
        """
        segments = []
        if self.clipped_start < self.reference_start:
            segments.append((self.clipped_start, self.reference_start - 1))
        if self.clipped_end > self.reference_end:
            segments.append((self.reference_end, self.clipped_end - 1))
        return segments

    def getSegments(self, segment_type: ReadSegmentType) -> Iterator[ReadSegment]:
        """
        The segments for the requested kind of read sub-region.

        :param segment_type: Which sub-region to yield: the mapped portion or
            the non-overlapping clip extensions.
        :return: An iterator of inclusive segments.
        """
        match segment_type:
            case ReadSegmentType.MAPPED:
                yield self.getMappedSegment()
            case ReadSegmentType.CLIPPED:
                yield from self.getClippedSegments()
