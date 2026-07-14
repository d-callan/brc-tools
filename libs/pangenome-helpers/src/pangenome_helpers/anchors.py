"""Anchor-prep orchestration helpers."""

from __future__ import annotations

from pathlib import Path

from genome_io.gff import (
    build_isoforms,
    collect_protein_coding_genes,
    filter_bed12,
)


class AnchorPrepResult(tuple):
    """Named tuple summarizing anchor-prep outputs."""

    __slots__ = ()

    @property
    def bed_total(self) -> int:  # pragma: no cover - property wrappers are trivial
        return self[0]

    @property
    def bed_kept(self) -> int:  # pragma: no cover
        return self[1]

    @property
    def isoforms(self) -> int:  # pragma: no cover
        return self[2]


def prepare_anchor_inputs(
    gff_path: str | Path,
    raw_bed_path: str | Path,
    out_bed_path: str | Path,
    out_isoforms_path: str | Path,
) -> AnchorPrepResult:
    """Filter the raw gffread BED12 + emit isoforms TSV for TOGA2 anchors."""

    pc_genes = collect_protein_coding_genes(gff_path)
    n_in, n_kept = filter_bed12(raw_bed_path, out_bed_path, pc_genes)
    n_iso = build_isoforms(gff_path, out_isoforms_path)
    return AnchorPrepResult((n_in, n_kept, n_iso))
