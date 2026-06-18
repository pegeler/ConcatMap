# ADR 0002: Drop secondary alignments, keep supplementary

**Date:** 2026-06-18
**Status:** Accepted

## Context

minimap2 with `--sam-hit-only` can emit three alignment classes for a single
read:

- **Primary** (flag 0x0): the best placement.
- **Secondary** (flag 0x100): an alternative placement of the same read at a
  different locus. Marks a read that maps to multiple locations.
- **Supplementary** (flag 0x800): a split-read segment --- a different portion
  of the same read's bases mapped to a separate location (chimeric alignment).

v1.x applied no flag filter, pulling all three classes with
`fetch(until_eof=True)`.

## Decision

Drop secondary alignments; keep supplementary alignments.

## Rationale

**Secondary alignments** represent the same molecule redrawn at a different
locus. Including them inflates apparent read density (2054 of 4838 records in
the Del1 test dataset were secondary, ~42%) without adding biological
information. Standard coverage tools (`samtools depth`, `samtools flagstat`)
drop secondary by default for the same reason.

**Supplementary alignments** carry the complementary segment of a split read.
They are essential for visualizing origin-spanning reads (the primary use-case
for the concatenated-reference design) and deletion-junction reads. Dropping
them would make the plot blind to the very events ConcatMap is built to detect.

This policy matches `samtools depth`'s default, keeping the depth colormap
and the plotted lines coherent with respect to flag filtering (see ADR 0001
for the remaining divergence).

## Consequences

- Read counts drop substantially on datasets with many alternative alignments.
- Supplementary alignments require correct CIGAR-based clip projection to
  render accurately (hard clips are common in split-read records); see the
  CIGAR clip fix in `mapper.py`.
- Hard break from v1.x behavior, noted in a comment at the filter site.
