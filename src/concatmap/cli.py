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
"""
import argparse
import logging
from pathlib import Path

from concatmap import mapper
from concatmap import plot
from concatmap import utils


def get_logger(debug: bool) -> logging.Logger:
    log_level = logging.DEBUG if debug else logging.INFO

    logger = logging.getLogger(__package__)
    logger.setLevel(log_level)

    ch = logging.StreamHandler()
    ch.setLevel(log_level)
    logger.addHandler(ch)

    return logger


class ArgumentParser(argparse.ArgumentParser):
    """
    Argument parser that explains the ``--unsorted``/``--depth`` conflict.

    argparse routes a mutually-exclusive group violation through ``error()``
    with a terse "not allowed with" message. We intercept it to append the
    rationale, since there is only one such group in this parser.
    """

    def error(self, message: str) -> None:
        if '--depth' in message and '--unsorted' in message:
            message += (
                ' (--depth computes coverage with `samtools depth`, which '
                'requires a position-sorted sam file)'
            )
        super().error(message)


def parse_args(argv=None) -> argparse.Namespace:
    p = ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        '-q', '--query_file',
        required=True,
        type=Path,
        help='Input path for fastq sequencing reads file.')
    p.add_argument(
        '-r', '--reference_file',
        required=True,
        type=Path,
        help='Input path for fasta reference file.')
    p.add_argument(
        '-o', '--output_dir',
        type=Path,
        help='Output directory for sam file and plot. '
             '(default: same folder as input fastq file)')
    p.add_argument(
        '-n', '--output_name',
        type=Path,
        help='Output name for sam file and plot. '
             '(default: same name as input fastq file)')
    p.add_argument(
        '-m', '--min_length',
        type=int,
        default=100,
        help='Minimum mapped read length to plot. (default: %(default)d)')
    p.add_argument(
        '-l', '--line_spacing',
        type=float,
        default=0.02,
        help='Radial spacing of each read on plot. (default: %(default).2f)')
    p.add_argument(
        '-w', '--line_width',
        type=float,
        default=0.75,
        help='Line width of each read on plot. (default: %(default).2f)')
    p.add_argument(
        '-c', '--circle_size',
        type=float,
        default=0.45,
        help='Size of central circle. (default: %(default).2f)')
    p.add_argument(
        '-s', '--fig_size',
        type=float,
        default=10,
        help='Size of figure. (default: %(default).2f)')
    p.add_argument(
        '-x', '--include_clipped_reads',
        action='store_true',
        help='Plot clipped portion of reads.')
    p.add_argument(
        '-f', '--figure_format',
        default='pdf',  # using string for default in help entry
        type=plot.OutputFormat.__getitem__,
        metavar=f'{{{",".join([ext.name for ext in plot.OutputFormat])}}}',
        choices=list(plot.OutputFormat),
        help='Format of saved figure. (default: %(default)s)')
    g = p.add_mutually_exclusive_group()
    g.add_argument(
        '-u', '--unsorted',
        action='store_true',
        help='Plot from unsorted sam file. Mutually exclusive with --depth.')
    g.add_argument(
        '-d', '--depth',
        action='store_true',
        help='Plot read line segments colored by read depth at each position. '
             'Requires a position-sorted sam file, so it cannot be combined '
             'with --unsorted.')
    p.add_argument(
        '--debug',
        action='store_true',
        help=argparse.SUPPRESS)

    args = p.parse_args(argv)

    if not args.output_dir:
        args.output_dir = args.query_file.parent

    if not args.output_name:
        args.output_name = args.query_file.stem
    args.output_file_stem = args.output_dir / args.output_name

    if args.debug:
        utils.Debug().is_debug = True

    return args


def main():
    args = parse_args()

    output_dir: Path = args.output_dir
    if output_dir.exists() and not output_dir.is_dir():
        raise RuntimeError(f'{output_dir} is not a directory')
    output_dir.mkdir(parents=True, exist_ok=True)

    mapper.concatmap(args, get_logger(args.debug))


if __name__ == '__main__':
    main()
