"""Graph co-membership edge orchestration (Phase E graph edges)."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
from typing import Iterable

from genome_io.bed import load_bed_genes_by_source
from genome_io.pansn import load_graph_paths


def compute_graph_edges(
    paths_tsv: str | Path,
    annotations_pattern: str,
    strains: Iterable[str],
) -> list[dict]:
    """Return path co-membership edges limited to the provided strains."""

    allowed = set(strains)
    genes_by_strain = load_bed_genes_by_source(annotations_pattern)
    paths = load_graph_paths(paths_tsv)
    edges: list[dict] = []

    for path_id, strain_set in paths.items():
        strain_list = sorted(s for s in strain_set if s in allowed)
        if len(strain_list) < 2:
            continue
        for sa, sb in combinations(strain_list, 2):
            genes_a = [g for g in genes_by_strain.get(sa, []) if g[0] == path_id]
            genes_b = [g for g in genes_by_strain.get(sb, []) if g[0] == path_id]
            for _, _, _, gid_a in genes_a:
                for _, _, _, gid_b in genes_b:
                    edges.append(
                        {
                            "strain_a": sa,
                            "gene_a": gid_a,
                            "strain_b": sb,
                            "gene_b": gid_b,
                            "path_id": path_id,
                            "overlap": "1.000",
                        }
                    )
    return edges
