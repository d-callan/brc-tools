"""Reciprocal-best overlap orchestration (Phase E rbest edges)."""

from __future__ import annotations

import glob
from collections import defaultdict
from pathlib import Path
from typing import Iterable

from genome_io.bed import load_bed_genes_by_source
from genome_io.chains import iter_chains, project_gene
from genome_io.intervals import best_query_gene, index_by_chrom


class RbestEdge(dict):
    """Dictionary subclass for type clarity in tests."""


def compute_rbest_edges(
    chains_pattern: str,
    annotations_pattern: str,
    *,
    min_overlap: float = 0.90,
) -> list[RbestEdge]:
    """Compute reciprocal-best edges from chain projections.

    Parameters
    ----------
    chains_pattern : str
        Glob pattern for ``*.rbest.chain`` files.
    annotations_pattern : str
        Glob pattern for per-strain BED annotations.
    min_overlap : float, optional
        Minimum fractional coverage required on both target/query projections.
    """

    genes_by_strain = load_bed_genes_by_source(annotations_pattern)
    indexed_queries = {
        strain: index_by_chrom(entries) for strain, entries in genes_by_strain.items()
    }
    edges: list[RbestEdge] = []
    seen: set[tuple[str, str, str, str]] = set()

    for chain_path in sorted(glob.glob(chains_pattern)):
        stem = Path(chain_path).stem.replace(".rbest", "")
        parts = stem.split(".")
        if len(parts) < 2:
            continue
        strain_a, strain_b = parts[0], parts[1]
        q_index = indexed_queries.get(strain_b, {})
        genes_a_by_chrom = defaultdict(list)
        for chrom, start, end, gene_id in genes_by_strain.get(strain_a, []):
            genes_a_by_chrom[chrom].append((start, end, gene_id))

        for header, blocks in iter_chains(chain_path):
            if not blocks:
                continue
            chrom = header["tName"]
            q_chrom = header["qName"]
            for start, end, gene_id in genes_a_by_chrom.get(chrom, []):
                if end <= header["tStart"] or start >= header["tEnd"]:
                    continue
                aligned, qmin, qmax = project_gene(start, end, blocks)
                target_frac = aligned / max(1, end - start)
                if target_frac < min_overlap or qmin is None or qmax is None:
                    continue
                best = best_query_gene(qmin, qmax, q_index.get(q_chrom))
                if not best or best[1] < min_overlap:
                    continue
                key = (strain_a, gene_id, strain_b, best[0])
                if key in seen:
                    continue
                seen.add(key)
                edges.append(
                    RbestEdge(
                        strain_a=strain_a,
                        gene_a=gene_id,
                        strain_b=strain_b,
                        gene_b=best[0],
                        overlap_a=f"{target_frac:.3f}",
                        overlap_b=f"{best[1]:.3f}",
                    )
                )
    return edges
