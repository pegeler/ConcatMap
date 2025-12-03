"""
Driver code for processing files and generating plots.
"""
import logging
import subprocess
from argparse import Namespace
from collections.abc import Iterable
from collections.abc import Iterator
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pysam

from concatmap import plot
from concatmap.struct import PolarCoordinate
from concatmap.struct import PolarLineSegment
from concatmap.struct import SamFileRead
from concatmap.utils import PositionToAngleConverter


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
) -> Iterator[SamFileRead]:
    """
    Generator that reads a *sam* file and yields a ``SamFileRead`` that contains
    the reference start and end positions as well as the clipped base positions.

    :param sam_filename: The file path.
    :param unsorted: Whether to leave the *sam* file unsorted or run the
            ``samtools`` sort.
    :param min_length: Minimum read length to be retained.
    :return: An iterator of ``SamFileRead``.
    """
    if unsorted:
        samfile = pysam.AlignmentFile(sam_filename, 'r')
    else:
        sorted_sam_filename = sam_filename.with_stem(sam_filename.stem + '_sorted')
        pysam.samtools.sort(
            '-o', str(sorted_sam_filename), str(sam_filename), catch_stdout=False)
        samfile = pysam.AlignmentFile(sorted_sam_filename, 'r')

    for r in samfile.fetch(until_eof=True):
        if r.reference_length < min_length:
            continue
        # Start position of clipped bases relative to the reference upstream
        qs0 = r.reference_start - r.query_alignment_start
        # End position of clipped bases relative to the reference downstream
        qe1 = r.reference_end + r.infer_read_length() - r.query_alignment_end
        yield SamFileRead(r.reference_start, r.reference_end, qs0, qe1)


def convert_reads_to_line_segments(
        reads: Iterable[SamFileRead],
        reference_length: int,
        line_spacing: float,
        basis_radius: float,
) -> Iterator[PolarLineSegment]:
    pos_to_angle_converter = PositionToAngleConverter(reference_length)
    for i, read in enumerate(reads):
        line_segment = PolarLineSegment(
            PolarCoordinate(
                pos_to_angle_converter(read.reference_start),
                basis_radius + line_spacing * i),
            PolarCoordinate(
                pos_to_angle_converter(read.reference_end),
                basis_radius + line_spacing * i))
        yield line_segment


def concatmap(args: Namespace, logger: logging.Logger) -> None:
    logger.info('Reading reference file %s', args.reference_file)
    reference_record = SeqIO.read(args.reference_file, 'fasta')

    concat_record = concat_sequence(reference_record)
    concat_filename = args.output_file_stem.with_name(f'{concat_record.id}.fasta')
    logger.info('Writing concatenated file %s', concat_filename)
    with open(concat_filename, 'w') as fh:
        SeqIO.write(concat_record, fh, 'fasta')

    sam_filename = args.output_file_stem.with_suffix('.sam')

    logger.info('Running minimap2')
    run_minimap(args.query_file, concat_filename, sam_filename)

    reads = read_samfile(sam_filename, args.unsorted, args.min_length)
    line_segments = convert_reads_to_line_segments(
        reads,
        len(reference_record),
        args.line_spacing,
        args.circle_size,
    )

    figure_file = args.output_file_stem.with_suffix(args.figure_format.value)
    logger.info('Plotting output to %s', figure_file)
    plot.plot(
        line_segments,
        args.fig_size,
        args.line_spacing,
        args.line_width,
        args.circle_size,
        args.clip,
        figure_file,
    )
