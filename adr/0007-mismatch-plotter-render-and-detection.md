# ADR 0007: `MismatchPlotter` render model and mismatch detection

**Date:** 2026-07-10
**Status:** Accepted

## Context

`MismatchPlotter` renders each read IGV/MSA alignment-track style: gray where
the read matches the reference, colored by nucleotide at each substituted base.
This required two decisions with non-obvious alternatives: how to *detect*
per-base substitutions, and how to *render* them.

The existing coloring plotter, `MulticolorLinePlotter` (used by `--depth`),
draws a read as ~200 linearized sample points and colors each pairwise segment
by interpolating an external per-position `values` array through
`AngularCoordinatesInterpolator` (`np.interp`) and a `LineCollection`. Reusing
that machinery for mismatches was the obvious first move, and the initial draft
did exactly that.

## Decision

### Detection: direct CIGAR walk, not MD tags

Detect substitutions by walking each read's CIGAR
(`mapper._find_mismatches`), comparing consumed query bases against the
in-memory concatenated reference string at each aligned (M/=/X) position.
Insertions/soft clips consume query only; deletions/ref-skips consume reference
only; hard clips/pads consume neither. Scope is substitutions only, consistent
with the project's one-base-per-reference-position model (README "Scope and
Limitations").

Rejected: MD tags / `get_aligned_pairs`. MD tags require minimap2 to emit them
(not guaranteed), and the reference is already loaded in `mapper.concatmap`, so
a direct comparison has no extra dependency and stays in this codebase's
coordinate space with no wraparound math.

The walk lives in `mapper.py` because it depends on `pysam` (CIGAR op
constants, `AlignedSegment`); the *representation* (`Mismatch`, and
`SamFileRead.mismatches`) lives in `struct.py`, which stays stdlib-only. This is
the same layering rule as ADR 0006: behavior lives with the data it describes,
and `plot.py`/`struct.py` never import `pysam`.

Detection is performed only when `--by_base` is requested:
`read_samfile(reference_sequence=...)` walks the CIGAR only when handed the
sequence, so default and depth runs pay nothing for an O(read-length)
comparison they never render. (Clip extents, by contrast, are O(1) per read and
are always computed — the cost asymmetry is why mismatches are gated and clips
are not.)

### Render: sparse per-mismatch overlays over a single gray body

Draw the whole read as one gray arc (as `DefaultPlotter` does), then overpaint
one short colored arc per substitution, spanning exactly that base's angular
slot `[conv(pos), conv(pos + 1)]` at the read's radius. Work is proportional to
the number of substitutions; a read with none is exactly the `DefaultPlotter`
path.

Rejected: the `MulticolorLinePlotter`-style fixed-stride
`LineCollection`. Depth is a smooth signal, so sampling it at a fixed ~200
points and interpolating is correct. **Substitutions are sparse, discrete,
single-base point features**, and a fixed stride aliases them: on a read longer
than ~200 bases each segment spans several bases, so a lone SNP either falls
between two sample midpoints and vanishes, or paints several bases wide at a
fractional-base offset. Rounding sample midpoints to the nearest integer
position (as the rejected draft did) only patches a discretization the sampling
itself introduces, and it reconstructs per-sample positions the pipeline had
already computed and discarded.

The overlay model is exact at any read length (each SNP at its true position and
true one-base width), efficient (O(substitutions), not a flat 200 segments per
read), and simpler: no `LineCollection`, no `set_rmax`, no position
reconstruction, no nearest-position rounding, and no `position -> base`
reverse-lookup dict (iterating the sparse `read.mismatches` directly asks the
cheaper question). It reuses the existing `PolarLineSegment` / `_linearize`
machinery rather than a parallel path. `AngularCoordinatesInterpolator` is
deliberately not reused — `np.interp` is meaningless for categorical color data.

### Threading read identity to the draw hook

`AbstractPlotter._convertReadsToLineSegments` yields a named `PlacedRead`
(read + `PolarLineSegment`), and `_drawLineSegment` takes an explicit, typed
`read` parameter. This replaces an earlier `(read, start, end,
PolarLineSegment)` positional tuple and an opaque `**kwargs` channel — both of
which ADR 0005 argues against (positional coupling to declaration order; intent
not obvious to a reader). Every mapped segment genuinely corresponds to a read,
so the contract "here is a segment and the read it belongs to; draw it" is
honest; `DefaultPlotter`/`MulticolorLinePlotter` simply ignore `read`.

## Consequences

- `SamFileRead` gains a defaulted `mismatches: tuple[Mismatch, ...] = ()` field;
  all existing positional constructions and plotters are unaffected.
- `MismatchPlotter` overrides `BASIS_COLOR` (→ black) and `CLIPPED_COLOR`
  (→ a neutral that is none of the four base colors), because red is the T
  substitution color; inheriting the red defaults would make the basis circle
  and clip extensions ambiguous with a T mismatch. (`MulticolorLinePlotter`
  overrides the same two for the analogous clash with its colormap.)
- Ambiguous query bases (e.g. `N`) are reported by the walk but have no entry in
  the color table, so they render as the gray body — one `dict.get` fallback,
  no special-casing.
- `-b` joins the existing `-u`/`-d` argparse mutually-exclusive group, so
  choosing exactly one view mode is enforced by argparse (the group scales to
  further modes without manual pairwise checks) and the usage synopsis reads
  `[-u | -d | -b]`. The tradeoff is that `-u -b` (unsorted *and* by-base) is
  rejected even though by-base detection has no sortedness requirement of its
  own; the group treats `-u` as a view-mode selector rather than an orthogonal
  sorting toggle. This was a deliberate choice for a single "pick one view"
  UX. If a future mode genuinely needs to combine with `-u`, revisit by
  splitting sorting from coloring (e.g. a `--color-by {none,depth,base}`
  argument alongside a standalone `-u`).
