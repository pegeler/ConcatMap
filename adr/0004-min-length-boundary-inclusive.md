# ADR 0004: `min_length` filter uses `>=` (inclusive lower bound)

**Date:** 2026-06-18
**Status:** Accepted

## Context

`read_samfile` filters out reads whose aligned reference span is shorter than
`min_length`. The boundary condition --- whether a read exactly equal to
`min_length` is kept --- differs between v1.x and v2:

- **v1.x:** `reference_length > min_length` (strict; a read equal to the
  minimum is excluded).
- **v2:** `reference_length >= min_length` (inclusive; a read equal to the
  minimum is kept).

## Decision

Keep `>=` in v2; do not revert to match v1.x.

## Rationale

`>=` is the more intuitive reading of "minimum length to plot": if the user
passes `-m 1000`, they expect reads of at least 1000 bases to appear. The
strict `>` behavior of v1.x is a one-read edge case with no biological
justification. This is a deliberate and acceptable hard break from v1.x.

## Consequences

- Reads whose reference span equals `min_length` exactly are now included
  where v1.x excluded them. The difference is at most one read per unique
  span length, negligible in practice.
- Migration Guide in the README notes the flag rename but not this boundary
  change, as it is unlikely to affect any real workflow.
