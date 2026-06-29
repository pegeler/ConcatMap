import math
from collections.abc import Iterator
from dataclasses import dataclass
from enum import auto
from enum import StrEnum
from typing import NamedTuple


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

    def getEndpoints(self) -> tuple[int, int]:
        """
        Inclusive reference start/end of the aligned portion of the read.

        :return: A tuple representing start and end, inclusive.
        """
        return self.reference_start, self.reference_end - 1

    def clipSegments(self) -> list[tuple[int, int]]:
        """
        Inclusive start/end pairs for each non-overlapping clip extension.

        :return: 0–2 pairs: leading extension (if any), then trailing (if any).
        """
        segments = []
        if self.clipped_start < self.reference_start:
            segments.append((self.clipped_start, self.reference_start - 1))
        if self.clipped_end > self.reference_end:
            segments.append((self.reference_end, self.clipped_end - 1))
        return segments

    def segments(
            self,
            segment_type: ReadSegmentType,
    ) -> Iterator[tuple[int, int]]:
        """
        Inclusive start/end pairs for the requested kind of read sub-region.

        :param segment_type: Which sub-region to yield: the mapped portion or
            the non-overlapping clip extensions.
        :return: An iterator of inclusive (start, end) pairs.
        """
        match segment_type:
            case ReadSegmentType.MAPPED:
                yield self.getEndpoints()
            case ReadSegmentType.CLIPPED:
                yield from self.clipSegments()
