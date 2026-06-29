# ADR 0006: Read sub-region dispatch via `SamFileRead.getSegments()`

**Date:** 2026-06-28
**Status:** Accepted

Supersedes [ADR 0005](0005-SamFileRead-getEndpoints-method-impl.md).

## Context

The plotter draws two kinds of read sub-region: the mapped portion (always)
and, with `-x`, the clipped extensions. Originally both were selected through a
single `getEndpoints(include_clipped: bool)` method that returned one
`(start, end)` pair (ADR 0005).

Two pressures broke that design:

1. **The clip case became 1:N.** To eliminate a moiré where the clip arc was
   overpainted by the mapped arc, the clip is now drawn as only its
   non-overlapping extensions --- a leading stub and/or a trailing stub. A read
   yields 0, 1, or 2 clip segments, but always exactly 1 mapped segment. A
   method returning a single pair can no longer serve both cases.

2. **Selection logic kept landing in the plotter.** Intermediate iterations
   put the mapped-vs-clipped dispatch in `plot.py` --- first as a
   `Callable` strategy argument, then as a private `_getSegments` staticmethod
   that took a `SamFileRead` and ignored `self`. Both are read logic wearing a
   plotter's clothes.

## Decision

Dispatch on a `ReadSegmentType` enum (`MAPPED` / `CLIPPED`) inside a new
`SamFileRead.getSegments(segment_type)` generator, and place **both** the enum
and the method in `struct.py`:

```python
type ReadSegment = tuple[int, int]


class ReadSegmentType(StrEnum):
    MAPPED = auto()
    CLIPPED = auto()


class SamFileRead(NamedTuple):
    ...
    def getSegments(
            self,
            segment_type: ReadSegmentType,
    ) -> Iterator[ReadSegment]:
        match segment_type:
            case ReadSegmentType.MAPPED:
                yield self.getMappedSegment()
            case ReadSegmentType.CLIPPED:
                yield from self.getClippedSegments()
```

`getMappedSegment()` and `getClippedSegments()` remain as named, individually
testable helpers; `getSegments()` is the type-directed dispatcher over them. The
plotter calls `read.getSegments(segment_type)` and never inspects which fields a
given segment type pulls. All three methods share one `verbObject` shape and a
single noun (`ReadSegment`), and the helper names mirror the enum members
(`MAPPED` -> `getMappedSegment`, `CLIPPED` -> `getClippedSegments`).
`ReadSegment` is a type alias for an inclusive `(start, end)` reference
interval --- the return shape of all three methods.

## Rationale

**Enum over `bool` flag.** A boolean parameter conflates two distinct
behaviors at the call site (`include_clipped=True` reads as a tweak, not a mode
switch). `ReadSegmentType.CLIPPED` is self-documenting and extensible if a
third sub-region is ever needed. The `match/case` is the project's idiom for
type-directed behavior over `if/elif` chains.

**Enum over `Callable` strategy.** A callable extractor
(`lambda r: [r.getMappedSegment()]`) was rejected: it pushes awkward machinery onto
every call site for the common default and hides what is really a small, closed
set of modes.

**Method on the data, enum in the data module.** `plot.py` imports from
`struct.py`, never the reverse; `struct.py` depends only on the stdlib. Read
sub-region selection is a property of a read, not of how it is rendered, so
`ReadSegmentType` is a domain concept that belongs in `struct.py`. Defining it
there lets `getSegments()` dispatch locally without inverting the layering, and it
keeps read logic out of the plotter --- behavior lives with the data it
describes.

## Consequences

- `getEndpoints(include_clipped)` is replaced by `getMappedSegment()`, a plain
  accessor for the mapped inclusive segment with no flag.
- The plotter's `_convertReadsToLineSegments` takes a `ReadSegmentType`
  (defaulting to `MAPPED`) and iterates `read.getSegments(segment_type)`,
  handling the 1:1 and 1:N cases through one loop.
- The over-length guard (ADR 0003) now runs per clip extension rather than over
  the whole clipped span.
- Adding a new read sub-region means adding a `ReadSegmentType` member and a
  `match` arm in `struct.py`; the plotter needs no change.
