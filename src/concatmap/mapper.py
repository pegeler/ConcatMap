"""
Driver code for processing files and generating plots.
"""
import functools
import logging
import subprocess
from argparse import Namespace
from collections.abc import Iterator
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pysam

from concatmap import plot
from concatmap.struct import SamFileRead

# Leading/trailing CIGAR operations representing clipped (unaligned) bases.
CLIP_OPS = frozenset({pysam.CSOFT_CLIP, pysam.CHARD_CLIP})

# TODO: This module is starting to look like it would benefit from class encapsulation.


def concat_sequence(record: SeqRecord) -> SeqRecord:
    """
    Concatenates a ``SeqRecord`` to itself.

    :param record: The ``SeqRecord`` to be doubled.
    :return: A new ``SeqRecord`` that is a concatenation of the original to
            itself. The new record ID will have '_CONCAT' appended to it.
    """
    new_id = record.id + '_CONCAT'
    return SeqRecord(
        record.seq * 2,
        id=new_id,
        name=new_id,
        description='concatenated reference file for ConcatMap')


def run_minimap(query_file: Path, reference_file: Path, sam_file: Path) -> None:
    """
    Makes a call out to the ``minimap2`` CLI to create a *sam* file from a
    query and reference file.

    Note, expects ``minimap2`` command to be on the search path.

    :param query_file: The location of the query *fastq* file.
    :param reference_file: The location of the reference *fasta* file.
    :param sam_file: The output *sam* file location.
    :return: Nothing. This function is called for its side effect.
    """
    cmd = [
        'minimap2', '--sam-hit-only', '-a',
        reference_file,
        query_file,
        '-o', sam_file,
    ]
    subprocess.run(cmd, check=True)


def read_samfile(
        sam_filename: Path,
        unsorted: bool,
        min_length: int,
        logger: logging.Logger | None = None,
) -> Iterator[SamFileRead]:
    """
    Generator that reads a *sam* file and yields a ``SamFileRead`` that contains
    the reference start and end positions as well as the clipped base positions.

    :param sam_filename: The file path.
    :param unsorted: Whether to leave the *sam* file unsorted or run the
            ``samtools`` sort.
    :param min_length: Minimum read length to be retained.
    :param logger: Optional logger.
    :return: An iterator of ``SamFileRead``.
    """
    if unsorted:
        samfile = pysam.AlignmentFile(sam_filename, 'r')
    else:
        sorted_sam_filename = sam_filename.with_stem(sam_filename.stem + '_sorted')
        if logger:
            logger.info('Writing sorted sam file to %s', sorted_sam_filename)
        pysam.samtools.sort(
            '-o', str(sorted_sam_filename), str(sam_filename), catch_stdout=False)
        samfile = pysam.AlignmentFile(sorted_sam_filename, 'r')

    for r in samfile.fetch(until_eof=True):
        # Skip secondary alignments (alternative placements of the same read).
        # Supplementary (split-read) alignments are kept: they carry the other
        # half of origin-spanning and deletion-junction reads. NOTE: v1.x
        # filtered neither. ``samtools depth`` (used by --depth) applies the
        # same policy by default, so the depth and line views agree on which
        # reads contribute.
        if not r.is_mapped or r.is_secondary or r.reference_length < min_length:
            continue
        # Take clip extents from the leading/trailing CIGAR clip ops so soft (S)
        # and hard (H) clips are handled uniformly. NOTE: do not derive these
        # from query_alignment_start/end and infer_read_length() --- those mix
        # read-coordinate systems (H excluded vs. included) and misproject a
        # leading hard clip onto the downstream end.
        cigar = r.cigartuples
        lead = cigar[0][1] if cigar[0][0] in CLIP_OPS else 0
        trail = cigar[-1][1] if cigar[-1][0] in CLIP_OPS else 0
        # Clipped-base positions projected onto the reference up/downstream.
        clipped_start = r.reference_start - lead
        clipped_end = r.reference_end + trail
        yield SamFileRead(
            r.reference_start,
            r.reference_end,
            clipped_start,
            clipped_end,
        )


def get_depths_at_positions(depth_file: Path, reference_length: int) -> list[float]:
    """
    Reads a depth file produced by ``samtools depth`` and returns a list of
    counts for each position in the reference sequence. Concatenated reference
    sequences are handled correctly, as positions will wrap around.

    :param depth_file: The tab-separated file output by ``samtools depth``.
    :param reference_length: The length of the reference sequence.
    :return: A list of counts for each position in the reference sequence.
    """
    counts = [0.] * reference_length
    with open(depth_file, 'r') as fh:
        for line in fh:
            pos, count = map(int, line.strip().split('\t')[1:])
            counts[(pos - 1) % reference_length] += count
    return counts


def concatmap(args: Namespace, logger: logging.Logger) -> None:
    logger.info('Reading reference file %s', args.reference_file)
    reference_record = SeqIO.read(args.reference_file, 'fasta')

    base_path = args.output_file_stem

    concat_record = concat_sequence(reference_record)
    # Scope the concatenated fasta to this run's output stem. NOTE: v1.x named
    # it after the reference sequence id (ignoring -n), so repeated runs against
    # one reference overwrote a single shared file. The record id still carries
    # the '_CONCAT' suffix and becomes the reference name in the sam file.
    concat_filename = base_path.with_name(f'{base_path.stem}_concat.fasta')
    logger.info('Writing concatenated file %s', concat_filename)
    with open(concat_filename, 'w') as fh:
        SeqIO.write(concat_record, fh, 'fasta')

    sam_filename = base_path.with_suffix('.sam')

    logger.info('Running minimap2')
    run_minimap(args.query_file, concat_filename, sam_filename)

    reads = list(read_samfile(sam_filename, args.unsorted, args.min_length, logger))
    logger.info('Read %d records from sam file', len(reads))

    if args.depth:
        # TODO: This should be refactored into two functions. Or better yet,
        #       we should make a ConcatMapper class and these could be private
        #       methods within it.
        sorted_sam_filename = sam_filename.with_stem(sam_filename.stem + '_sorted')
        depth_filename = base_path.with_name(base_path.stem + '_depth.tsv')
        logger.info('Writing depth file to %s', depth_filename)
        # samtools depth's default flag filtering matches read_samfile (drops
        # secondary, keeps supplementary). NOTE: -l filters on read length while
        # the plotted lines filter on reference span, so the depth and line
        # views can disagree slightly for heavily clipped reads.
        pysam.samtools.depth(
            '-aa',
            str(sorted_sam_filename),
            '-l', str(args.min_length),
            '-o', str(depth_filename),
        )
        depths = get_depths_at_positions(depth_filename, len(reference_record))
        plotter_class = functools.partial(plot.MulticolorLinePlotter, values=depths)
    else:
        plotter_class = plot.DefaultPlotter

    figure_file = base_path.with_suffix(args.figure_format.value)
    logger.info('Plotting output to %s', figure_file)
    plotter = plotter_class(
        reads=reads,
        reference_length=len(reference_record),
        fig_size=args.fig_size,
        line_spacing=args.line_spacing,
        line_width=args.line_width,
        circle_size=args.circle_size,
        include_clipped_reads=args.include_clipped_reads,
        figure_file=figure_file,
    )
    plotter.plot()
