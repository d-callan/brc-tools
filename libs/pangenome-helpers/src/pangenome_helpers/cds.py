"""Orthogroup CDS grouping orchestration (group_cds_by_og)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator

from genome_io.gff import normalize_gene_id, parse_gff_cds
from genome_io.orthology import load_ortho_table
from genome_io.sequence import (
    extract_cds,
    load_fasta_as_dict,
    strip_internal_stops,
    translate,
)

from .manifest import ensure_matching_collections, ensure_reference_not_in_queries, load_manifest_map


@dataclass(frozen=True)
class OrthogroupSequences:
    """CDS/protein sequences keyed by strain for a single orthogroup."""

    orthogroup_id: str
    reference_gene_id: str
    cds: Dict[str, str]
    proteins: Dict[str, str]


class OrthogroupFilterError(ValueError):
    """Raised when orthogroup inputs are inconsistent (e.g., missing reference gene)."""


def _load_fasta_cache(manifest: Dict[str, str]) -> Dict[str, Dict[str, str]]:
    return {strain: load_fasta_as_dict(path) for strain, path in manifest.items()}


def _extract_cds_from_segments(segments, fasta_dict) -> str:
    return extract_cds(segments, fasta_dict) if segments else ""


def iter_orthogroup_sequences(
    ortho_table: str | Path,
    ref_strain: str,
    ref_gff: str | Path,
    ref_fasta: str | Path,
    query_gff_manifest: str | Path,
    query_fasta_manifest: str | Path,
    *,
    min_intact: int = 2,
) -> Iterator[OrthogroupSequences]:
    """Yield orthogroup CDS/protein maps, mirroring group_cds_by_og logic."""

    query_gff_map = load_manifest_map(query_gff_manifest)
    query_fasta_map = load_manifest_map(query_fasta_manifest)
    ensure_reference_not_in_queries(ref_strain, query_gff_map)
    ensure_matching_collections(query_gff_manifest, query_fasta_manifest)

    strains = [ref_strain] + sorted(query_gff_map)
    ref_segments = parse_gff_cds(ref_gff)
    ref_fasta_dict = load_fasta_as_dict(ref_fasta)
    ortho_pairs = load_ortho_table(
        ortho_table,
        strains,
        min_intact,
        ref_strain,
        normalize_fn=normalize_gene_id,
    )

    target_genes = {gene_id for _, gene_id in ortho_pairs if gene_id in ref_segments}
    query_cds_maps: Dict[str, dict] = {}
    for strain, gff_path in query_gff_map.items():
        query_cds_maps[strain] = parse_gff_cds(gff_path, target_genes)

    query_fastas = _load_fasta_cache(query_fasta_map)

    for og_id, ref_gene_id in ortho_pairs:
        segments = ref_segments.get(ref_gene_id)
        if not segments:
            continue
        ref_cds = _extract_cds_from_segments(segments, ref_fasta_dict)
        if not ref_cds or len(ref_cds) % 3 != 0:
            continue
        ref_prot = translate(ref_cds).rstrip("*")
        if "*" in ref_prot:
            continue

        cds_map = {ref_strain: ref_cds}
        prot_map = {ref_strain: ref_prot}

        for strain in query_gff_map:
            segments = query_cds_maps[strain].get(ref_gene_id, [])
            cds_seq = _extract_cds_from_segments(segments, query_fastas[strain])
            if not cds_seq:
                continue
            cds_seq = cds_seq[: (len(cds_seq) // 3) * 3]
            cds_seq = strip_internal_stops(cds_seq)
            if not cds_seq:
                continue
            prot_seq = translate(cds_seq).rstrip("*")
            cds_map[strain] = cds_seq
            prot_map[strain] = prot_seq

        if len(cds_map) - 1 < min_intact:
            continue

        yield OrthogroupSequences(og_id, ref_gene_id, cds_map, prot_map)
