# ADR 0003: Skip clip arcs that exceed one full turn

**Date:** 2026-06-18
**Status:** Accepted

## Context

ConcatMap projects clipped (unaligned) bases linearly onto the circular
reference: one clipped base = one reference position upstream or downstream of
the alignment. Each read contributes up to two non-overlapping clip extensions
(a leading stub before the alignment and a trailing stub after it), drawn as
separate arcs. When an extension is longer than the reference, its projected
arc sweeps past 2pi and overlaps itself, producing a misleading ring that
looks like uniform full coverage.

This occurs in practice with nanopore concatemers, chimeric reads, and long
reads that align only locally (observed in 458 of 2155 records, ~21%, on the
Del1 nanopore dataset).

Three options were considered:

1. **Skip** --- omit the clip arc entirely; still draw the mapped portion.
2. **Cap** --- truncate the clip projection at the reference boundary (wrap to
   at most one turn).
3. **Raise** --- error out or warn per-read; treat as unsupported input.

## Decision

Skip (option 1): for each clip extension, if
`abs(thetas[-1] - thetas[0]) > math.tau`, do not draw that extension's arc.
The guard is applied per extension, so a read with one in-scope extension and
one over-length extension still draws the in-scope one. The mapped read arc is
always plotted.

## Rationale

**Cap** was rejected because it would draw the clip as reaching exactly to
the origin, which is geometrically arbitrary and biologically misleading ---
the clip does not actually end there.

**Raise** was rejected because the affected reads are a minority and their
mapped portions are still informative. Stopping the run or flooding the log
for ~20% of reads is disproportionate.

**Skip** is honest: it omits only what cannot be drawn without fabricating
information. The user sees the aligned arc and can interpret it normally; the
absence of a clip arc signals "clip is out of scope" rather than "no clip
exists." The README "Scope and Limitations" section documents the case and the
affected read types.

The guard is a single `if` before `ax.plot`, so complexity cost is negligible.

## Consequences

- Over-length clip extensions for concatemers/chimeras are silently omitted
  (no warning emitted per-read to avoid log spam). A summary count could be
  added later. Because the guard is per extension, a read keeps any extension
  that is individually in scope --- earlier the whole clipped span was measured
  as one arc, so such a read was dropped wholesale.
- The tool is not appropriate for datasets where over-length clips are the
  primary signal of interest.
- `PositionToAngleConverter.__call__` intentionally returns unwrapped angles
  (outside `[0, tau)`) so the span check works on the raw difference; the
  modulo is applied only inside `AngularCoordinatesInterpolator.__call__`.
