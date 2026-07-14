"""PanSN orchestration helpers (file IO + validation)."""

from __future__ import annotations

import gzip
from pathlib import Path

from genome_io.io import open_maybe_gz
from genome_io.pansn import rename_headers


class PansnError(ValueError):
    """Raised when inputs make PanSN header renaming impossible."""


def validate_sample(sample: str, delimiter: str) -> None:
    if delimiter in sample:
        raise PansnError(f"sample name {sample!r} contains delimiter {delimiter!r}")
    if any(c.isspace() for c in sample):
        raise PansnError(f"sample name {sample!r} contains whitespace")
    if not sample:
        raise PansnError("sample name is empty")


def rename_fasta(
    input_path: str | Path,
    output_path: str | Path,
    sample: str,
    *,
    haplotype: int = 1,
    delimiter: str = "#",
    gzip_output: bool = False,
) -> int:
    """Rename FASTA headers using PanSN SAMPLE/HAP/CONTIG prefixes.

    Returns the number of headers rewritten. Raises :class:`PansnError` if no
    headers are found.
    """

    validate_sample(sample, delimiter)
    in_fh = open_maybe_gz(input_path, "rt")
    try:
        if gzip_output:
            out_fh = gzip.open(output_path, "wt")
        else:
            out_fh = open(output_path, "wt")
        try:
            n = rename_headers(in_fh, out_fh, sample, haplotype, delimiter)
        finally:
            out_fh.close()
    finally:
        in_fh.close()
    if n == 0:
        raise PansnError("no FASTA headers found in input")
    return n
