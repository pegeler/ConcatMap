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
                 [-n OUTPUT_NAME] [-m MIN_LENGTH] [-u] [-l LINE_SPACING]
                 [-w LINE_WIDTH] [-c CIRCLE_SIZE] [-s FIG_SIZE] [-x]
                 [-f {eps,jpeg,jpg,pdf,pgf,png,ps,raw,rgba,svg,svgz,tif,tiff}]

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
                        Output directory for sam file and plot (default: same
                        folder as input fastq file)
  -n OUTPUT_NAME, --output_name OUTPUT_NAME
                        Output name for sam file and plot (default: same name
                        as input fastq file)
  -m MIN_LENGTH, --min_length MIN_LENGTH
                        Minimum mapped read length to plot (default: 100)
  -u, --unsorted        Plot from unsorted sam file
  -l LINE_SPACING, --line_spacing LINE_SPACING
                        Radial spacing of each read on plot (default: 0.02)
  -w LINE_WIDTH, --line_width LINE_WIDTH
                        Line width of each read on plot (default: 0.75)
  -c CIRCLE_SIZE, --circle_size CIRCLE_SIZE
                        Size of central circle (default: 0.45)
  -s FIG_SIZE, --fig_size FIG_SIZE
                        Size of figure (default: 10.00)
  -x, --clip            Plot clipped portion of reads
  -f {eps,jpeg,jpg,pdf,pgf,png,ps,raw,rgba,svg,svgz,tif,tiff}, --figure_format {eps,jpeg,jpg,pdf,pgf,png,ps,raw,rgba,svg,svgz,tif,tiff}
                        Format of saved figure. (default: pdf)
```

### Example

```bash
concatmap -q <PathToFASTQFile/InputFileName> \
          -r <PathToFASTAReferenceFile/ReferenceFileName> \
          -m 10000 \
          -f "png"
```

## Gallery

[![DOI: 10.1093/g3journal/jkab423](https://cdn.ncbi.nlm.nih.gov/pmc/blobs/9816/9210306/ed0056d2704b/jkab423f2.jpg)](https://pubmed.ncbi.nlm.nih.gov/34897429/#&gid=article-figures&pid=figure-2-uid-1)

## Reference

<!-- TODO: Link to paper and include some figures -->
