# ADR 0005: Implementation of `SamFileRead.getEndpoints()` method

**Date:** 2026-06-19
**Status:** Accepted

## Context

Three candidate implementations were considered for dispatching on `include_clipped`
to return the correct start/end pair from `SamFileRead`.

### Using `getattr()`

Construct attribute names from `StrEnum`s and access fields dynamically.

```python
import enum
from typing import NamedTuple


class EndpointType(enum.StrEnum):
    MAPPED = 'mapped'
    CLIPPED = 'clipped'


class TerminusType(enum.StrEnum):
    START = 'start'
    END = 'end'


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

    def getEndpoints(self, include_clipped: bool = False) -> tuple[int, int]:
        """
        Inclusive reference start/end of the read.

         :param include_clipped: Whether to include just mapped bases or
             also clipped bases.
        :return: A tuple representing start and end, inclusive.
        """
        endpoint_type = EndpointType.CLIPPED if include_clipped else EndpointType.MAPPED
        return (
            self._getEndpoint(endpoint_type, TerminusType.START),
            self._getEndpoint(endpoint_type, TerminusType.END) - 1
        )

    def _getEndpoint(
            self,
            endpoint_type: EndpointType,
            terminus_type: TerminusType,
    ) -> int:
        return getattr(self, f'{endpoint_type}_{terminus_type}')
```

### Using `self` with slices

Exploit the `NamedTuple` sequence interface, selecting fields by positional slice.

```python
from typing import NamedTuple


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

    def getEndpoints(self, include_clipped: bool = False) -> tuple[int, int]:
        """
        Inclusive reference start/end of the read.

         :param include_clipped: Whether to include just mapped bases or
             also clipped bases.
        :return: A tuple representing start and end, inclusive.
        """
        slice_ = slice(2, 4) if include_clipped else slice(2)
        start, end = self[slice_]
        return start, end - 1
```

### Using explicit `if/else` with named attributes

Branch on the flag and access fields by name directly.

```python
from typing import NamedTuple


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

    def getEndpoints(self, include_clipped: bool = False) -> tuple[int, int]:
        """
        Inclusive reference start/end of the read.

         :param include_clipped: Whether to include just mapped bases or
             also clipped bases.
        :return: A tuple representing start and end, inclusive.
        """
        if include_clipped:
            return self.clipped_start, self.clipped_end - 1
        return self.reference_start, self.reference_end - 1
```

## Decision

Using explicit `if/else` with named attributes.

## Rationale

The `getattr()` approach is the most extensible --- it would scale gracefully if
`SamFileRead` grew additional read types --- but it pays for that generality with
significant machinery (two `StrEnum`s, a private helper) that buys nothing given
the current two-pair structure. Rejected as over-engineered.

The slice approach is compact, but it silently couples to field declaration order.
Adding or reordering fields breaks the method with no type error or warning, and
the intent of `slice(2, 4)` is not obvious to a reader. Rejected as fragile.

The explicit `if/else` is slightly more verbose than the slice approach but names
its fields directly, making it immune to reordering and immediately self-documenting.
It is the simplest implementation that is also correct under future structural
changes to `SamFileRead`.
