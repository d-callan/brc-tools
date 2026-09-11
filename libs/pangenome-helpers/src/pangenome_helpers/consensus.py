"""Consensus orthogroup orchestration (Phase E).

Nodes are NATIVE genes. A projected annotation keeps the reference gene's
id -- Liftoff writes PVPAM_130008300 into PvC01's GFF -- so a projection is
resolved back onto the native gene at the same locus (via gene_beds) and
contributes an EDGE. One that resolves to nothing is dropped, because adding
it as a node would mean a second node for a locus that already has one.

Without the resolution step each physical gene enters the graph twice,
max_copies is 2 in 83% of groups, and since the labels key on max_copies
only 21 of 5,817 orthogroups come out CORE-1:1 (0.4% for eight conspecific
P. vivax strains). With the fix the same inputs give 3,990 of 5,731 (69.6%),
matching what the rbest edges alone produce (68.6%).
"""

from __future__ import annotations

import csv
import glob as _glob
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


def load_gene_beds(spec: str | Path) -> dict[str, dict[str, list[tuple[int, int, str]]]]:
    """Load per-strain native gene BEDs into ``{strain: {chrom: [(start, end, gene_id)]}}``.

    The BED filename stem is the strain name, matching phase_e_rbest_overlap.
    """
    idx: dict[str, dict[str, list[tuple[int, int, str]]]] = {}
    paths = _glob.glob(str(spec)) if any(c in str(spec) for c in "*?[") else \
        [str(p) for p in Path(spec).iterdir()] if Path(spec).is_dir() else [str(spec)]
    for path in sorted(paths):
        strain = Path(path).stem
        by_chrom: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
        with open(path) as fh:
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                if len(f) >= 4:
                    try:
                        by_chrom[f[0]].append((int(f[1]), int(f[2]), f[3]))
                    except ValueError:
                        continue
        for c in by_chrom:
            by_chrom[c].sort()
        if by_chrom:
            idx[strain] = dict(by_chrom)
    return idx


def resolve_to_native(
    gene_index: dict[str, dict[str, list[tuple[int, int, str]]]],
    strain: str,
    chrom: str,
    start: int,
    end: int,
    min_overlap: float,
) -> str | None:
    """Best-overlapping native gene id at these coordinates, or None.

    A projected annotation keeps the REFERENCE gene's id, so the projected id
    never equals the native one. Without this lookup each physical gene enters
    the graph twice, every group shows two copies per strain, and the labels
    (which key on max_copies) become artifacts.
    """
    best_id, best_ov = None, 0.0
    for ns, ne, gid in gene_index.get(strain, {}).get(chrom, ()):
        if ne <= start:
            continue
        if ns >= end:
            break
        ov = min(end, ne) - max(start, ns)
        if ov <= 0:
            continue
        r = min(ov / max(1, end - start), ov / max(1, ne - ns))
        if r > best_ov:
            best_ov, best_id = r, gid
    return best_id if best_ov >= min_overlap else None


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
    used_edges: set[frozenset[str]],
    gene_index: dict[str, dict[str, list[tuple[int, int, str]]]] | None = None,
    alias_overlap: float = 0.5,
    keep_unresolved_projections: bool = False,
    split_many2many: bool = False,
) -> dict[str, int]:
    """Seed union-find edges from Liftoff classification rows.

    Returns a dict of stats: aliases, aliased_to_native, unresolved_projections,
    many2many_skipped.
    """
    pos_records: dict[tuple[str, str], list[tuple[int, int, str]]] = defaultdict(list)
    liftoff_dir = Path(liftoff_dir)

    aliased_to_native = 0
    unresolved_projections = 0
    many2many_skipped = 0

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

        chrom = row.get("query_chrom", "")
        start_s = row.get("query_start", "")
        end_s = row.get("query_end", "")
        coords: tuple[str, int, int] | None = None
        if chrom and start_s and end_s:
            try:
                coords = (chrom, int(start_s), int(end_s))
            except ValueError:
                coords = None

        # Resolve projected gene to native gene by position
        if gene_index and coords:
            native = resolve_to_native(gene_index, strain, coords[0], coords[1],
                                       coords[2], alias_overlap)
            if native:
                query_gene = normalize_gene_id(native)
                aliased_to_native += 1
            else:
                unresolved_projections += 1
                if not keep_unresolved_projections:
                    continue

        ocls = (row.get("orthology_class") or "").strip()
        if split_many2many and ocls == "many2many":
            many2many_skipped += 1
            continue

        anchor_node = f"{anchor}#{ref_gene}"
        query_node = f"{strain}#{query_gene}"
        uf.union(anchor_node, query_node)
        used_edges.add(frozenset((anchor_node, query_node)))
        if coords:
            pos_records[(strain, coords[0])].append((coords[1], coords[2], query_node))
            node_positions[query_node] = coords

    # Interval-based aliasing (90% reciprocal overlap -> same physical gene)
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

    return {
        "aliases": aliases,
        "aliased_to_native": aliased_to_native,
        "unresolved_projections": unresolved_projections,
        "many2many_skipped": many2many_skipped,
    }


def _ingest_extra_edges(
    uf: UnionFind,
    edges: Iterable[dict] | None,
    used_edges: set[frozenset[str]],
) -> int:
    count = 0
    for row in edges or []:
        sa = row.get("strain_a", "")
        sb = row.get("strain_b", "")
        ga = normalize_gene_id(row.get("gene_a", ""))
        gb = normalize_gene_id(row.get("gene_b", ""))
        if sa and sb and ga and gb:
            node_a = f"{sa}#{ga}"
            node_b = f"{sb}#{gb}"
            uf.union(node_a, node_b)
            used_edges.add(frozenset((node_a, node_b)))
            count += 1
    return count


def build_consensus_table(
    liftoff_dir: str | Path,
    anchors: Sequence[str],
    strains: Sequence[str],
    *,
    rbest_edges: Iterable[dict] | None = None,
    gene_beds: str | Path | None = None,
    alias_overlap: float = 0.5,
    keep_unresolved_projections: bool = False,
    split_many2many: bool = False,
) -> list[ConsensusResult]:
    """Return consensus orthogroup rows (orthogroup_id + strain columns).

    Parameters
    ----------
    gene_beds
        Glob pattern or directory of per-strain native gene BED files.
        When provided, projected genes are resolved to native genes by
        reciprocal overlap (default 0.5). Without this, projected genes
        keep their reference-derived ids and each physical gene enters
        the graph twice, inflating max_copies.
    alias_overlap
        Minimum reciprocal overlap for resolving a projected gene to a
        native gene (default: 0.5).
    keep_unresolved_projections
        If True, add a projection that resolves to no native gene as a node
        of its own. Default False: drop it, since it is a locus the
        annotation missed, not an extra copy.
    split_many2many
        If True, skip projections whose orthology_class is "many2many".
        TOGA2 emits that as a warning that the genes should probably not
        collapse into one group; union-find can never undo the merge.
    """
    anchors = list(anchors)
    strains = list(strains)
    uf = UnionFind()
    node_positions: dict[str, tuple[str, int, int]] = {}
    used_edges: set[frozenset[str]] = set()

    gene_index = load_gene_beds(gene_beds) if gene_beds else None

    _seed_edges(
        liftoff_dir, anchors, strains, uf, node_positions, used_edges,
        gene_index=gene_index,
        alias_overlap=alias_overlap,
        keep_unresolved_projections=keep_unresolved_projections,
        split_many2many=split_many2many,
    )
    _ingest_extra_edges(uf, rbest_edges, used_edges)

    # Clique completeness: count edges per component
    comp_edges: dict[str, int] = Counter()
    for e in used_edges:
        comp_edges[uf.find(next(iter(e)))] += 1

    # Connected components
    components: dict[str, set[str]] = defaultdict(set)
    for node in list(uf.parent.keys()):
        components[uf.find(node)].add(node)

    rows: list[ConsensusResult] = []
    total_strains = len(strains)
    for cid, comp_nodes in sorted(components.items(), key=lambda kv: (-len(kv[1]), min(kv[1]))):
        per_strain: dict[str, list[str]] = defaultdict(list)
        for node in sorted(comp_nodes):
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

        # Clique: score against the support 1:1 evidence can actually provide.
        # Sum over strain pairs of min(copies_a, copies_b). Single-copy groups
        # reduce to plain density; multi-copy groups stay comparable.
        counts = sorted(len(c) for c in strain_clusters.values())
        expected_edges = sum(
            min(counts[i], counts[j])
            for i in range(len(counts))
            for j in range(i + 1, len(counts))
        )
        observed = comp_edges.get(cid, 0)
        clique = round(min(observed / expected_edges, 1.0), 3) if expected_edges else 1.0

        # Density: raw edge fraction over all member pairs.
        n_nodes = len(comp_nodes)
        pairs = n_nodes * (n_nodes - 1) // 2
        density = round(observed / pairs, 4) if pairs else 1.0

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
            clique=clique,
            density=density,
        )
        for strain in strains:
            if strain in strain_clusters:
                row[strain] = ",".join("|".join(cluster) for cluster in strain_clusters[strain])
            else:
                row[strain] = "-"
        rows.append(row)

    return rows


def summarize_labels(rows: Sequence[dict]) -> Counter:
    """Return counts per label for reporting/logging."""

    return Counter(row["label"] for row in rows)
