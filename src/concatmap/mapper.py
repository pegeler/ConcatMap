"""
Driver code for processing files and generating plots.
"""
import math
import subprocess
from argparse import Namespace
from collections.abc import Iterable
from collections.abc import Iterator
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import NamedTuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from matplotlib import pyplot as plt
import numpy as np
import pysam


class OutputFormat(Enum):  # pylint: disable=invalid-name
    eps = '.eps'
    jpeg = '.jpeg'
    jpg = '.jpg'
    pdf = '.pdf'
    pgf = '.pgf'
    png = '.png'
    ps = '.ps'
    raw = '.raw'
    rgba = '.rgba'
    svg = '.svg'
    svgz = '.svgz'
    tif = '.tif'
    tiff = '.tiff'


class PositionToAngleConverter:
    """
    Callable class that will convert a position to the angle component of
    polar coordinates relative to the reference sequence length.
    """

    def __init__(self, reference_length: int) -> None:
        self.reference_length = reference_length
        self.deg_per_base = 360 / reference_length

    def __call__(self, pos: int) -> float:
        return pos * self.deg_per_base


@dataclass(slots=True, frozen=True)
class PolarCoordinate:
    """
    Store polar coordinate pairs.

    :ivar angle: The point's angular coordinate in degrees.
    :ivar radius: The point's distance from the pole.
    """
    angle: float
    radius: float

    @property
    def radians(self) -> float:
        """
        :return: Converted angular component in radians.
        """
        return self.angle * 2 * math.pi / 360

    def get_pair(self, radians: bool = False) -> tuple[float, float]:
        """
        Get the ordered pair of angle and radius.

        :param radians: Whether to report the angle in radians. Degrees is default.
        :return: An ordered pair of (angle, radius).
        """
        return self.radians if radians else self.angle, self.radius


@dataclass(slots=True, frozen=True)
class PolarLineSegment:
    """
    Store a polar coordinate line segment.

    :ivar start_coord: The start coordinate.
    :ivar end_coord: The end coordinate.
    """
    start_coord: PolarCoordinate
    end_coord: PolarCoordinate


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


class SamFileRead(NamedTuple):
    """
    A container to hold start and end positions from the samfile.
    """
    reference_start: int
    reference_end: int
    clipped_start: int
    clipped_end: int


def read_samfile(
        sam_filename: Path,
        unsorted: bool,
        min_length: int,
) -> Iterator[SamFileRead]:
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


def _plot_line_segment(
        ax: plt.Axes,
        line_segment: PolarLineSegment,
        n_points: int = 500,
        *args,
        **kwargs
) -> None:
    thetas = np.linspace(
        line_segment.start_coord.radians,
        line_segment.end_coord.radians,
        n_points)
    # NOTE: Originally radii was derived from `scipy.interp1d`, which has been
    #       superseded by `np.interp`. However, current use-case does not
    #       warrant it, as `np.linspace` should work.
    radii = np.linspace(
        line_segment.start_coord.radius,
        line_segment.end_coord.radius,
        n_points)
    ax.plot(thetas, radii, *args, **kwargs)


def plot(
        line_segments: Iterable[PolarLineSegment],
        fig_size: float,
        line_spacing: float,
        line_width: float,
        circle_size: float,
        clip: bool,
        figure_file: Path,
) -> None:
    with plt.style.context('ggplot'):
        fig = plt.figure(figsize=(fig_size, ) * 2)
        ax = fig.add_subplot(111, polar=True)
        ax.grid(False)
        ax.set_rticks([])
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_theta_zero_location('N')
        ax.set_facecolor('white')
        ax.axis('off')

        basis_radius = circle_size - 2 * line_spacing
        basis_curve = PolarLineSegment(
            PolarCoordinate(0, basis_radius),
            PolarCoordinate(360, basis_radius))
        _plot_line_segment(ax, basis_curve, linewidth=5)

        # TODO: Draw clipped reads

        for line_segment in line_segments:
            _plot_line_segment(ax, line_segment, color='grey', linewidth=line_width)

    # TODO: Flip image so it is cw instead of ccw?
    plt.savefig(figure_file, bbox_inches='tight')


def concatmap(args: Namespace) -> None:
    reference_record = SeqIO.read(args.reference_file, 'fasta')

    concat_record = concat_sequence(reference_record)
    with open(f'{concat_record.id}.fasta', 'w') as fh:
        SeqIO.write(concat_record, fh, 'fasta')

    sam_filename = args.output_file.with_suffix('.sam')

    run_minimap(args.query_file, Path(concat_record.id), sam_filename)

    reads = read_samfile(sam_filename, args.unsorted, args.min_length)
    line_segments = convert_reads_to_line_segments(
        reads,
        len(reference_record),
        args.line_spacing,
        args.circle_size,
    )

    plot(
        line_segments,
        args.fig_size,
        args.line_spacing,
        args.line_width,
        args.circle_size,
        args.clip,
        args.output_file.with_suffix(args.figure_format.value),
    )
