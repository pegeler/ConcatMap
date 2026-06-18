# ADR 0001: Read-depth calculation via `samtools depth` subprocess

**Date:** 2026-06-18
**Status:** Accepted

## Context

The `--depth` mode colors each read arc by the sequencing depth at each
position along the reference. Depth must be computed per base across the
concatenated reference, then folded back to the single-copy reference length.

Two approaches were considered:

1. **Call `samtools depth` as a subprocess** --- write the sorted SAM, invoke
   `pysam.samtools.depth`, parse the TSV output.
2. **Reimplement depth in Python** --- walk the same `SamFileRead` records
   yielded by `read_samfile` and accumulate coverage in a Python array.

## Decision

Use `samtools depth` (approach 1).

## Rationale

`samtools depth` is battle-tested for this exact calculation and available
wherever pysam is (pysam bundles samtools). Reimplementing it would require
correctly handling CIGAR operations (deletions, introns, clips) to expand
each alignment record into reference-coordinate spans, which is non-trivial
and a source of new bugs.

## Known shortfall

`samtools depth -l` filters on **read length** (total query bases, including
clipped bases), while `read_samfile`'s `min_length` guard filters on
**reference span** (`r.reference_length`, the aligned portion only). For
heavily clipped reads these two metrics diverge: such a read may pass one
filter and fail the other, causing the depth colormap and the plotted lines
to reflect slightly different read sets.

On typical datasets (organellar/plasmid, long-read, moderate clip rates) the
effect is negligible. The residual is documented in `mapper.py` at the
`samtools.depth(...)` call site. Closing it fully would require approach 2.

## Consequences

- Depth calculation is robust and maintainable.
- A sorted SAM file must exist before `samtools depth` is invoked; the
  `--depth` and `--unsorted` flags are mutually exclusive to enforce this.
- Minor read-set divergence between depth and line views remains for
  heavily clipped reads.
