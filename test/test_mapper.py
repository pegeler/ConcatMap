"""Unit tests for the pure CIGAR-walking mismatch detection in ``mapper``."""
from dataclasses import dataclass

import pysam

from concatmap.mapper import _find_mismatches
from concatmap.struct import Mismatch


@dataclass
class FakeSegment:
    """
    A minimal stand-in for ``pysam.AlignedSegment`` exposing only the three
    attributes ``_find_mismatches`` reads. Keeps the CIGAR walk testable
    without materializing a real alignment file.
    """
    query_sequence: str | None
    cigartuples: list[tuple[int, int]]
    reference_start: int


def find(query, cigartuples, reference, reference_start=0):
    segment = FakeSegment(query, cigartuples, reference_start)
    return list(_find_mismatches(segment, reference))


def test_all_match_yields_nothing():
    assert find('ACGT', [(pysam.CMATCH, 4)], 'ACGT') == []


def test_single_substitution():
    # read 'ACGT' vs ref 'AGGT': position 1 substituted C->G... read base is C.
    assert find('ACGT', [(pysam.CMATCH, 4)], 'AGGT') == [Mismatch(1, 'C')]


def test_run_of_substitutions():
    assert find('AAAA', [(pysam.CMATCH, 4)], 'TTAA') == [
        Mismatch(0, 'A'), Mismatch(1, 'A')]


def test_insertion_does_not_desync_downstream_positions():
    # 2M 2I 2M: the inserted 'TT' consumes query only; the trailing 'GG' must
    # still compare against reference positions 2,3 (not shifted by the insert).
    query = 'AA' 'TT' 'GG'
    reference = 'AAGG'
    assert find(query, [(pysam.CMATCH, 2), (pysam.CINS, 2), (pysam.CMATCH, 2)],
                reference) == []


def test_insertion_never_reported_as_mismatch():
    query = 'AA' 'TT' 'CC'  # trailing CC mismatches ref GG
    reference = 'AAGG'
    assert find(query, [(pysam.CMATCH, 2), (pysam.CINS, 2), (pysam.CMATCH, 2)],
                reference) == [Mismatch(2, 'C'), Mismatch(3, 'C')]


def test_deletion_advances_reference_only():
    # 2M 2D 2M: the deletion consumes reference positions 2,3; trailing query
    # 'GG' compares against reference positions 4,5.
    query = 'AAGG'
    reference = 'AAxxGG'
    assert find(query, [(pysam.CMATCH, 2), (pysam.CDEL, 2), (pysam.CMATCH, 2)],
                reference) == []


def test_leading_soft_clip_consumes_query_only():
    # 2S 4M: clipped 'NN' is skipped; alignment starts at reference_start.
    query = 'NN' 'ACGT'
    assert find(query, [(pysam.CSOFT_CLIP, 2), (pysam.CMATCH, 4)], 'ACGT') == []


def test_trailing_soft_clip_consumes_query_only():
    query = 'ACGT' 'NN'
    assert find(query, [(pysam.CMATCH, 4), (pysam.CSOFT_CLIP, 2)], 'ACGT') == []


def test_reference_start_offsets_positions():
    # Alignment starts at reference position 5; the mismatch is reported in
    # reference (concatenated) coordinates, not read coordinates.
    reference = 'xxxxxACGT'
    assert find('AGGT', [(pysam.CMATCH, 4)], reference, reference_start=5) == [
        Mismatch(6, 'G')]


def test_equal_op_skips_comparison():
    # '=' declares a match, so no comparison happens even if bases differ
    # (they shouldn't in valid output, but the op must be trusted, not walked).
    assert find('ACGT', [(pysam.CEQUAL, 4)], 'ZZZZ') == []


def test_diff_op_rederives_substituted_base():
    # 'X' declares a mismatch; we still read the actual query base to color it.
    assert find('ACGT', [(pysam.CDIFF, 4)], 'AAAA') == [
        Mismatch(1, 'C'), Mismatch(2, 'G'), Mismatch(3, 'T')]


def test_ambiguous_base_reported_and_left_for_plotter():
    # An 'N' query base at a matched position is a substitution against the ref
    # base; it is reported here and rendered gray by the plotter (no N color).
    assert find('ANGT', [(pysam.CMATCH, 4)], 'AAGT') == [Mismatch(1, 'N')]


def test_lowercase_query_is_upper_cased():
    # read_samfile upper-cases the reference; _find_mismatches upper-cases the
    # query. A lowercase soft-masked base must compare and report upper-cased.
    assert find('acgt', [(pysam.CMATCH, 4)], 'AAGT') == [Mismatch(1, 'C')]


def test_none_query_yields_nothing():
    assert find(None, [(pysam.CMATCH, 4)], 'ACGT') == []
