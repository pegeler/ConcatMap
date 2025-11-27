"""
Driver code for processing files and generating plots.
"""
import subprocess
from argparse import Namespace
from enum import Enum
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq  # Currently Unused
from Bio.SeqRecord import SeqRecord
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d  # Currently Unused
import numpy as np  # Currently Unused
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


def run_minimap(query_file: Path, reference_id_file: Path, sam_file: Path) -> None:
    cmd = [
        'minimap2', '--sam-hit-only', '-a',
        reference_id_file,
        query_file,
        '-o', sam_file,
    ]
    subprocess.run(cmd, check=True)


def concat_sequence(record: SeqRecord) -> SeqRecord:
    new_id = record.id + '_CONCAT'
    return SeqRecord(
        record.seq * 2,
        id=new_id,
        name=new_id,
        description='concatenated reference file for ConcatMap')


def read_samfile(sam_filename: Path, unsorted: bool, min_length: int) -> ...:
    # TODO: Currently considering the Array of Structs (AoS) vs Struct of Arrays
    #       (SoA) paradigm. This is related to Data Oriented Programming (DOP),
    #       which is popular in high performance applications, such as gaming.
    #       However, this will not be a performance bottleneck for us so I need
    #       to consider the best abstraction from a reader's standpoint.
    #       Note also to include this discussion in the manuscript.
    if unsorted:
        samfile = pysam.AlignmentFile(sam_filename, 'r')
    else:
        sorted_sam_filename = sam_filename.with_stem(sam_filename.stem + '_sorted')
        pysam.samtools.sort('-o', sorted_sam_filename, sam_filename, catch_stdout=False)
        samfile = pysam.AlignmentFile(sorted_sam_filename, 'r')

    for read in samfile.fetch(until_eof=True):
        if read.reference_length < min_length:
            continue
        ...

    return ...


class PolarCoordinatesConverter:
    """
    Callable class that will convert a position to polar coordinates relative
    to the reference sequence length.

    Note, coordinates will be converted to intuitive clockwise orientation
    rather than default counterclockwise.
    """

    def __init__(self, reference_length: int) -> None:
        self.reference_length = reference_length
        self.deg_per_base = 360 / reference_length

    def __call__(self, pos: int) -> float:
        return 360 - pos * self.deg_per_base


def concatmap(args: Namespace) -> None:
    reference_record = SeqIO.read(args.reference_file, 'fasta')

    concat_record = concat_sequence(reference_record)
    with open(concat_record.id, 'w') as fh:
        SeqIO.write(concat_record, fh, 'fasta')

    sam_filename = args.output_file.with_suffix('.sam')

    run_minimap(args.query_file, Path(concat_record.id), sam_filename)

    reads = read_samfile(sam_filename, args.unsorted, args.min_length)
    dedup_reference()  # TODO
    polar_coord_converter = PolarCoordinatesConverter(len(reference_record))
    convert_to_polar()  # TODO

    # Save figure
    figure_file = args.output_file.with_suffix(args.figure_format.value)
    plt.savefig(figure_file, bbox_inches='tight')
