"""
ConcatMap Command Line Interface
================================

Map and plot reads against a circular reference.

Details
-------

This program takes in a fasta reference file and a fastq sequencing file and
concatenates the reference sequence (two copies of the reference sequence
repeated in tandem) and then uses minimap2 to map reads to the concatenated
reference, generating a .sam file that can be used for subsequent plotting.

MIT License
-----------

Copyright (c) 2025 Daryl Gohl, Paul Egeler

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
import argparse

from pathlib import Path

from concatmap import mapper
from concatmap import utils


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        '-q', '--fastq_file',
        required=True,
        type=Path,
        help='Input path for fastq sequencing reads file.')
    p.add_argument(
        '-r', '--fasta_file',
        required=True,
        type=Path,
        help='Input path for fasta reference file.')
    p.add_argument(
        '-o', '--output_dir',
        type=Path,
        help='Output directory for sam file and plot '
             '(default: same folder as input fastq file)')
    p.add_argument(
        '-n', '--output_file',
        type=Path,
        help='Output name for sam file and plot '
             '(default: same name as input fastq file)')
    p.add_argument(
        '-m', '--min_length',
        type=int,
        default=100,
        help='Minimum mapped read length to plot (default: %(default)d)')
    p.add_argument(
        '-u', '--unsorted',
        action='store_true',
        help='Plot from unsorted sam file')
    p.add_argument(
        '-l', '--line_spacing',
        type=float,
        default=0.02,
        help='Radial spacing of each read on plot (default: %(default).2f)')
    p.add_argument(
        '-w', '--line_width',
        type=float,
        default=0.75,
        help='Line width of each read on plot (default: %(default).2f)')
    p.add_argument(
        '-c', '--circle_size',
        type=float,
        default=0.45,
        help='Size of central circle (default: %(default).2f)')
    p.add_argument(
        '-s', '--fig_size',
        type=float,
        default=10,
        help='Size of figure (default: %(default).2f)')
    p.add_argument(
        '-x', '--clip',
        action='store_true',
        help='Plot clipped portion of reads')
    p.add_argument(
        '-f', '--figure_format',
        default='pdf',  # using string for default in help entry
        type=mapper.OutputFormat.__getitem__,
        metavar=f'{{{",".join([ext.name for ext in mapper.OutputFormat])}}}',
        choices=list(mapper.OutputFormat),
        help='Format of saved figure. (default: %(default)s)')

    args = p.parse_args(argv)

    if not args.output_dir:
        args.output_dir = args.fastq_file.parent

    if not args.output_file:
        args.output_file = args.fastq_file
    args.output_file = args.output_file.with_suffix('.sam')

    return args


def main():
    args = parse_args()

    with utils.chdir_temporarily(args.output_dir):
        mapper.concatmap(args)


if __name__ == '__main__':
    main()
