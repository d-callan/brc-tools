from pathlib import Path

import pytest

from pangenome_helpers.pansn import PansnError, rename_fasta

DATA = Path(__file__).parent / "data" / "pansn"


def test_rename_fasta(tmp_path):
    output = tmp_path / "out.fa"
    n = rename_fasta(DATA / "input.fa", output, "SAMPLE", haplotype=2, delimiter="#")
    assert n == 2
    assert output.read_text() == ">SAMPLE#2#seq1\nACGT\n>SAMPLE#2#seq2 desc\nNNNN\n"


def test_rename_fasta_empty_headers_raises(tmp_path):
    empty = tmp_path / "empty.fa"
    empty.write_text("ACGT\n")
    with pytest.raises(PansnError):
        rename_fasta(empty, tmp_path / "out.fa", "sample")


def test_validate_sample_rules(tmp_path):
    out = tmp_path / "out.fa"
    with pytest.raises(PansnError):
        rename_fasta(DATA / "input.fa", out, "bad#sample")
    with pytest.raises(PansnError):
        rename_fasta(DATA / "input.fa", out, "bad sample")
