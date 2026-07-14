"""Orthogroup filtering helpers shared across pangenome stages."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Sequence

from genome_io.orthology import load_ortho_table


@dataclass(frozen=True)
class Orthogroup:
    """Simple container for reference-centric orthogroup metadata."""

    orthogroup_id: str
    reference_gene_id: str


def load_intact_orthogroups(
    ortho_table: str | Path,
    ref_strain: str,
    strains: Sequence[str],
    *,
    min_intact: int = 2,
    ref_genes: Iterable[str] | None = None,
    normalize: Callable[[str], str] | None = None,
) -> list[Orthogroup]:
    """Return orthogroups that pass the intactness filter."""

    entries = load_ortho_table(
        ortho_table,
        list(strains),
        min_intact,
        ref_strain,
        normalize_fn=normalize,
    )
    allowed = set(ref_genes) if ref_genes is not None else None
    orthogroups: list[Orthogroup] = []
    for og_id, ref_gene in entries:
        if allowed is not None and ref_gene not in allowed:
            continue
        orthogroups.append(Orthogroup(og_id, ref_gene))
    return orthogroups
