# ADR 0005: Implementation of `SamFileRead.getEndpoints()` method

**Date:** 2026-06-19
**Status:** Accepted

## Context

I could either use slices into `self` or construct the attribute from
`StrEnum`s and then use `getattr` to access the start/end position for
each read type. Here are the example implementations:

### Using `getattr()`

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

## Using `self` with slices

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

## Decision

Using `self` with slices.

## Rationale

The number and types of fields we are likely to store in `SameFileRead` are
limited right now. The `self` with slice option is simpler, but less explicit.
That is fine as long as `SamFileRead` is a simple container for two pairs of
values. If the data we collect becomes richer, we can revisit implementation.
