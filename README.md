ConcatMap
=========

This program takes in a _fasta_ reference file and a _fastq_ sequencing file and
concatenates the reference sequence (two copies of the reference sequence
repeated in tandem) and then uses [minimap2](https://github.com/lh3/minimap2)
to map reads to the concatenated reference, generating a _.sam_ file that is
used for subsequent plotting.

## Requirements

This program relies on [minimap2](https://github.com/lh3/minimap2), an external
tool for sequence alignment. You must install it separately before using
ConcatMap. Follow the instructions on the
[project _README.md_](https://github.com/lh3/minimap2?tab=readme-ov-file#installation)
to install it for your operating system.

## Install

You can install ConcatMap directly from the GitHub repository using pip.
It's recommended to do this within a virtual environment.

```bash
pip install git+https://github.com/darylgohl/ConcatMap.git
```

## Usage

```
usage: concatmap [-h] -q QUERY_FILE -r REFERENCE_FILE [-o OUTPUT_DIR]
                 [-n OUTPUT_NAME] [-m MIN_LENGTH] [-l LINE_SPACING]
                 [-w LINE_WIDTH] [-c CIRCLE_SIZE] [-s FIG_SIZE] [-x]
                 [-f {eps,jpeg,jpg,pdf,pgf,png,ps,raw,rgba,svg,svgz,tif,tiff}]
                 [-u | -d]

ConcatMap Command Line Interface
================================

Map and plot reads against a circular reference.

Details
-------

This program takes in a fasta reference file and a fastq sequencing file and
concatenates the reference sequence (two copies of the reference sequence
repeated in tandem) and then uses minimap2 to map reads to the concatenated
reference, generating a .sam file that can be used for subsequent plotting.

options:
  -h, --help            show this help message and exit
  -q QUERY_FILE, --query_file QUERY_FILE
                        Input path for fastq sequencing reads file.
  -r REFERENCE_FILE, --reference_file REFERENCE_FILE
                        Input path for fasta reference file.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory for sam file and plot. (default: same
                        folder as input fastq file)
  -n OUTPUT_NAME, --output_name OUTPUT_NAME
                        Output name for sam file and plot. (default: same name
                        as input fastq file)
  -m MIN_LENGTH, --min_length MIN_LENGTH
                        Minimum mapped read length to plot. (default: 100)
  -l LINE_SPACING, --line_spacing LINE_SPACING
                        Radial spacing of each read on plot. (default: 0.02)
  -w LINE_WIDTH, --line_width LINE_WIDTH
                        Line width of each read on plot. (default: 0.75)
  -c CIRCLE_SIZE, --circle_size CIRCLE_SIZE
                        Size of central circle. (default: 0.45)
  -s FIG_SIZE, --fig_size FIG_SIZE
                        Size of figure. (default: 10.00)
  -x, --include_clipped_reads
                        Plot clipped portion of reads.
  -f {eps,jpeg,jpg,pdf,pgf,png,ps,raw,rgba,svg,svgz,tif,tiff}, --figure_format {eps,jpeg,jpg,pdf,pgf,png,ps,raw,rgba,svg,svgz,tif,tiff}
                        Format of saved figure. (default: pdf)
  -u, --unsorted        Plot from unsorted sam file.
  -d, --depth           Plot read line segments colored by read depth at each
                        position.
```

### Examples

Every invocation needs a query (`-q`) and a reference (`-r`); the flags below
layer on different views of the same data.

**Basic plot.** Map reads and render the mapped portions as a PNG, dropping
reads shorter than 1 kb:

```bash
concatmap -q reads.fastq -r reference.fasta -m 1000 -f png
```

**Include clipped bases.** Add `-x` to also draw each read's soft-clipped
extensions beyond its aligned span (useful for spotting concatemers and
reads that wrap the origin):

```bash
concatmap -q reads.fastq -r reference.fasta -m 1000 -f png -x
```

**Color by read depth.** `-d` colors each read by per-position coverage instead
of drawing it in a flat color. This computes coverage with `samtools depth`, so
it requires a position-sorted sam file (the default):

```bash
concatmap -q reads.fastq -r reference.fasta -m 1000 -f png -d
```

**Skip sorting for a quick look.** `-u` plots from the unsorted sam file,
skipping the position sort. It is mutually exclusive with `-d`, which needs
sorted input:

```bash
concatmap -q reads.fastq -r reference.fasta -m 1000 -f png -u
```

**Tune the layout.** The geometry flags control the rendering: `-s` figure
size, `-c` central circle size, `-l` radial spacing between reads, `-w` line
width. Write the output somewhere specific with `-o`/`-n`:

```bash
concatmap -q reads.fastq -r reference.fasta -m 1000 -f png \
          -s 50 -c 0.5 -l 0.02 -w 1.5 \
          -o plots -n my_plasmid
```

## Gallery

[![DOI: 10.1093/g3journal/jkab423](https://cdn.ncbi.nlm.nih.gov/pmc/blobs/9816/9210306/ed0056d2704b/jkab423f2.jpg)](https://pubmed.ncbi.nlm.nih.gov/34897429/#&gid=article-figures&pid=figure-2-uid-1)

## Scope and Limitations

ConcatMap projects each read (and, with `-x`, its clipped bases) linearly onto a
circular reference, one read base per reference position. This assumption breaks
down in a few cases:

- **Clips longer than the reference.** The clipped portion of a read is
  unaligned, so projecting it along the reference is only meaningful when it is
  shorter than the reference. Reads whose clipped span exceeds the reference
  length --- for example nanopore concatemers, chimeric/fusion reads, or long
  reads that align only locally --- would sweep past a full turn and overlap
  themselves into a misleading ring. Their clip arc is omitted (the mapped
  portion is still plotted). This tool is not intended for such reads.
- **Large references.** Bacterial chromosomes and other large genomes produce a
  dense, hard-to-read spiral; ConcatMap targets small circular references
  (organellar genomes, plasmids, viral genomes).
- **Very high read counts.** Each read occupies its own concentric ring, so deep
  datasets eventually overplot. Subsample, or use `-m` to raise the minimum
  read length.

## Migration Guide

The command-line interface was reworked in version 2.0. The boolean string
toggles of the legacy interface (for example `-u T` / `-x f`) are now ordinary
flags: present means on, absent means off. All other options keep their short
names and meanings. The table below maps the legacy flags to the new ones.

| Legacy            | New        | Notes                                                |
| ----------------- | ---------- | ---------------------------------------------------- |
| `-u T` / `-u F`   | `-u` / _omit_ | `--unsorted` is now a flag; sorted is the default. |
| `-x t` / `-x f`   | `-x` / _omit_ | `--include_clipped_reads`; present means on.       |
| `-f "png"`        | `-f png`   | Now a validated format name; quotes optional.        |
| `-m -c -l -w -s`  | unchanged  | Same short names and meanings.                       |
| `-q -r -o -n`     | unchanged  | Same short names and meanings.                       |
| _(none)_          | `-d`       | New `--depth` mode; colors reads by per-position depth. Mutually exclusive with `-u`. |

For example, the legacy invocation

```bash
ConcatMap -u F -x t -f "png" -m 1000 -q reads.fastq -r ref.fasta -n out
```

becomes

```bash
concatmap -x -f png -m 1000 -q reads.fastq -r ref.fasta -n out
```

(`-u F` is dropped because sorted is the default, and `-x t` becomes a bare
`-x`).

## Reference

<!-- TODO: Link to paper and include some figures -->
