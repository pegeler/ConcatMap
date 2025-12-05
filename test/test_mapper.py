import logging
import subprocess
from pathlib import Path
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam

from concatmap.mapper import (
    concat_sequence,
    run_minimap,
    read_samfile,
    compute_coverage,
    convert_reads_to_line_segments,
)
from concatmap.struct import SamFileRead, PolarLineSegment


def test_concat_sequence_doubles_sequence():
    record = SeqRecord(Seq("ATGC"), id="ref")
    concat = concat_sequence(record)
    assert concat.id.endswith("_CONCAT")
    assert str(concat.seq) == "ATGCATGC"


def test_run_minimap_invokes_subprocess(monkeypatch, tmp_path):
    # Mock subprocess.run so we don't actually call minimap2
    called = {}
    def fake_run(cmd, check):
        called["cmd"] = cmd
        return 0
    monkeypatch.setattr(subprocess, "run", fake_run)

    q = tmp_path / "reads.fastq"
    r = tmp_path / "ref.fasta"
    s = tmp_path / "out.sam"
    q.write_text("")
    r.write_text(">ref\nATGC\n")

    run_minimap(q, r, s)
    assert "minimap2" in called["cmd"]


def test_read_samfile_and_compute_coverage(tmp_path):
    # Create a tiny SAM file with one alignment using pysam
    sam_path = tmp_path / "test.sam"
    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 10, "SN": "ref"}]}
    with pysam.AlignmentFile(sam_path, "w", header=header) as outf:
        a = pysam.AlignedSegment()
        a.query_name = "read1"
        a.query_sequence = "ATGCATGCAT"
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 0
        a.mapping_quality = 60
        a.cigar = [(0, 10)]  # 10M
        a.query_qualities = pysam.qualitystring_to_array("FFFFFFFFFF")
        outf.write(a)

    reads = list(read_samfile(sam_path, unsorted=True, min_length=5))
    assert isinstance(reads[0], SamFileRead)
    assert reads[0].reference_start == 0
    assert reads[0].reference_end == 10

    coverage = compute_coverage(reads, reference_length=10)
    assert len(coverage) == 10
    assert max(coverage) > 0


def test_convert_reads_to_line_segments():
    reads = [SamFileRead(1, 5, 0, 0), SamFileRead(2, 6, 0, 0)]
    segments = list(convert_reads_to_line_segments(reads, reference_length=10,
                                                   line_spacing=0.1,
                                                   basis_radius=0.5))
    assert all(isinstance(seg, PolarLineSegment) for seg in segments)
    assert len(segments) == 2
