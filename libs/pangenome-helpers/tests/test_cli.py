"""Tests for pangenome-helpers CLI."""

from pathlib import Path
from unittest.mock import patch

import pytest

from pangenome_helpers.cli import main

DATA = Path(__file__).parent / "data"


def test_cli_help(capsys):
    """Test that --help works."""
    with pytest.raises(SystemExit) as exc_info:
        main(["-h"])
    assert exc_info.value.code == 0
    captured = capsys.readouterr()
    assert "pangenome-helpers" in captured.out
    assert "triage" in captured.out


def test_cli_no_args(capsys):
    """Test that no args prints help."""
    result = main([])
    assert result == 0
    captured = capsys.readouterr()
    assert "pangenome-helpers" in captured.out


def test_cli_triage(tmp_path):
    """Test triage command."""
    output_dir = tmp_path / "triage_out"
    result = main([
        "triage",
        str(DATA / "triage" / "liftoff.gff3"),
        str(DATA / "triage" / "query.fa"),
        str(DATA / "triage" / "reference.bed"),
        str(output_dir),
        "test_query",
        "--family-list", str(DATA / "triage" / "family.tsv"),
    ])
    assert result == 0
    assert (output_dir / "triage.tsv").exists()
    assert (output_dir / "summary.json").exists()


def test_cli_anchor_prep(tmp_path):
    """Test anchor-prep command."""
    output_bed = tmp_path / "filtered.bed"
    output_iso = tmp_path / "iso.tsv"
    result = main([
        "anchor-prep",
        str(DATA / "anchors" / "sample.gff3"),
        str(DATA / "anchors" / "raw.bed"),
        str(output_bed),
        str(output_iso),
    ])
    assert result == 0
    assert output_bed.exists()
    assert output_iso.exists()


def test_cli_pansn_rename(tmp_path):
    """Test pansn-rename command."""
    output = tmp_path / "out.fa"
    result = main([
        "pansn-rename",
        str(DATA / "pansn" / "input.fa"),
        str(output),
        "SAMPLE",
        "--haplotype", "2",
    ])
    assert result == 0
    assert output.exists()
    content = output.read_text()
    assert "SAMPLE#2#" in content


def test_cli_process_maf(tmp_path):
    """Test process-maf command."""
    output = tmp_path / "out.maf"
    result = main([
        "process-maf",
        "REF",
        str(DATA / "maf" / "sample.maf"),
        str(output),
    ])
    assert result == 0
    assert output.exists()


def test_cli_multiz_order(tmp_path):
    """Test multiz-order command."""
    output = tmp_path / "order.txt"
    result = main([
        "multiz-order",
        str(DATA / "multiz" / "compare.csv"),
        "hinge",
        "query1", "query2",
        str(output),
    ])
    assert result == 0
    assert output.exists()
    lines = output.read_text().strip().split("\n")
    assert len(lines) == 2


def test_cli_triage_invalid_file(tmp_path, capsys):
    """Test triage with invalid file."""
    result = main([
        "triage",
        "/nonexistent/file.gff3",
        str(DATA / "triage" / "query.fa"),
        str(DATA / "triage" / "reference.bed"),
        str(tmp_path),
        "test",
    ])
    assert result == 1
    captured = capsys.readouterr()
    assert "failed" in captured.err.lower()
