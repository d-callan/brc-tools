"""Consensus orthogroup orchestration (Phase E)."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Sequence

from genome_io.gff import normalize_gene_id
from genome_io.orthology import (
    UnionFind,
    collapse_positions,
    edge_weight,
    reciprocal_overlap,
)


class ConsensusResult(dict):
    """Individual orthogroup row produced by :func:`build_consensus_table`."""


def _iter_classification_rows(liftoff_dir: Path, anchors: Sequence[str], strains: Sequence[str]):
    base = Path(liftoff_dir)
    for anchor in anchors:
        sub = base / f"{anchor}-as-ref"
        if not sub.exists():
            continue
        for strain in strains:
            if strain == anchor:
                continue
            cls_path = sub / f"{strain}.classification.tsv"
            if not cls_path.exists():
                continue
            with open(cls_path) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    yield anchor, strain, row


def _seed_edges(
    liftoff_dir: str | Path,
    anchors: Sequence[str],
    strains: Sequence[str],
    uf: UnionFind,
    node_positions: dict[str, tuple[str, int, int]],
):
    pos_records: dict[tuple[str, str], list[tuple[int, int, str]]] = defaultdict(list)
    liftoff_dir = Path(liftoff_dir)

    for anchor, strain, row in _iter_classification_rows(liftoff_dir, anchors, strains):
        ref_gene = normalize_gene_id(row.get("reference_gene_id", ""))
        query_gene = normalize_gene_id(row.get("query_gene_id", ""))
        src = row.get("source", "")
        intact = row.get("intactness", "")
        if not query_gene or query_gene == "None":
            continue
        w = edge_weight(src, intact)
        if w == 0:
            continue
        anchor_node = f"{anchor}#{ref_gene}"
        query_node = f"{strain}#{query_gene}"
        uf.union(anchor_node, query_node)
        chrom = row.get("query_chrom", "")
        start_s = row.get("query_start", "")
        end_s = row.get("query_end", "")
        if chrom and start_s and end_s:
            try:
                start_i = int(start_s)
                end_i = int(end_s)
            except ValueError:
                continue
            pos_records[(strain, chrom)].append((start_i, end_i, query_node))
            node_positions[query_node] = (chrom, start_i, end_i)

    aliases = 0
    for (_strain, _chrom), recs in pos_records.items():
        recs.sort()
        for i in range(len(recs)):
            si, ei, ni = recs[i]
            for j in range(i + 1, len(recs)):
                sj, ej, nj = recs[j]
                if sj > ei:
                    break
                if uf.find(ni) == uf.find(nj):
                    continue
                if reciprocal_overlap((si, ei), (sj, ej)) >= 0.9:
                    uf.union(ni, nj)
                    aliases += 1
    return aliases


def _ingest_extra_edges(uf: UnionFind, edges: Iterable[dict]):
    count = 0
    for row in edges or []:
        sa = row.get("strain_a", "")
        sb = row.get("strain_b", "")
        ga = normalize_gene_id(row.get("gene_a", ""))
        gb = normalize_gene_id(row.get("gene_b", ""))
        if sa and sb and ga and gb:
            uf.union(f"{sa}#{ga}", f"{sb}#{gb}")
            count += 1
    return count


def build_consensus_table(
    liftoff_dir: str | Path,
    anchors: Sequence[str],
    strains: Sequence[str],
    ref_strain: str,
    *,
    rbest_edges: Iterable[dict] | None = None,
    graph_edges: Iterable[dict] | None = None,
) -> list[ConsensusResult]:
    """Return consensus orthogroup rows (orthogroup_id + strain columns)."""

    anchors = list(anchors)
    strains = list(strains)
    uf = UnionFind()
    node_positions: dict[str, tuple[str, int, int]] = {}

    aliases = _seed_edges(liftoff_dir, anchors, strains, uf, node_positions)
    _ingest_extra_edges(uf, rbest_edges)
    _ingest_extra_edges(uf, graph_edges)

    components: dict[str, set[str]] = defaultdict(set)
    for node in list(uf.parent.keys()):
        components[uf.find(node)].add(node)

    rows: list[ConsensusResult] = []
    total_strains = len(strains)
    for comp_nodes in components.values():
        per_strain: dict[str, list[str]] = defaultdict(list)
        for node in comp_nodes:
            parts = node.split("#", 1)
            if len(parts) == 2:
                per_strain[parts[0]].append(parts[1])
        if len(per_strain) < 2:
            continue
        strain_clusters = {
            strain: collapse_positions(genes, strain, node_positions)
            for strain, genes in per_strain.items()
        }
        present = [s for s in strains if s in strain_clusters]
        n_strains = len(present)
        max_copies = max(len(clusters) for clusters in strain_clusters.values())
        if n_strains == total_strains and max_copies == 1:
            label = "CORE-1:1"
        elif n_strains == total_strains and max_copies >= 2:
            label = "CORE-VAR"
        elif max_copies >= 3:
            label = "FAMILY"
        elif n_strains <= 2:
            label = "LINEAGE-SPECIFIC"
        else:
            label = "PARTIAL"
        row: ConsensusResult = ConsensusResult(
            orthogroup_id=f"OG{len(rows)+1:06d}",
            label=label,
            n_strains=n_strains,
            max_copies=max_copies,
        )
        for strain in strains:
            if strain in strain_clusters:
                row[strain] = ",".join("|".join(cluster) for cluster in strain_clusters[strain])
            else:
                row[strain] = "-"
        rows.append(row)

    rows.sort(key=lambda r: r["orthogroup_id"])
    return rows


def summarize_labels(rows: Sequence[dict]) -> Counter:
    """Return counts per label for reporting/logging."""

    return Counter(row["label"] for row in rows)
