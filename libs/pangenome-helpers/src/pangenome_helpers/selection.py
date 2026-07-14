"""Selection track helpers (BUSTED + orthogroup BED builders)."""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
import tarfile
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping


@dataclass(slots=True)
class OrthogroupMapping:
    """Reference gene -> orthogroup metadata."""

    gene_id: str
    orthogroup_id: str
    label: str
    n_strains: int


def load_sizes(path: str | Path) -> dict[str, int]:
    """Load chromosome sizes from a ``.fa.fai`` or 2-column TSV."""

    sizes: dict[str, int] = {}
    with open(path) as handle:
        for line in handle:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            sizes[parts[0]] = int(parts[1])
    return sizes


def load_bed12(path: str | Path) -> dict[str, list[str]]:
    """Return gene_id -> BED12 fields, keeping the first isoform per gene."""

    bed: dict[str, list[str]] = {}
    with _open_maybe_gzip(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            name = fields[3]
            gene = name.rsplit(".", 1)[0] if "." in name else name
            bed.setdefault(gene, fields[:12])
    return bed


def load_ortholog_table(
    path: str | Path,
    ref_column: str,
    gene_prefix: str,
) -> tuple[dict[str, tuple[str, str, int]], dict[str, tuple[str, int]]]:
    """Return (gene -> (og_id, label, n_strains), og_id -> (label, n_strains))."""

    og_map: dict[str, tuple[str, str, int]] = {}
    og_info: dict[str, tuple[str, int]] = {}
    with _open_maybe_gzip(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        if ref_column not in fieldnames:
            raise ValueError(
                f"reference column {ref_column!r} missing from ortholog table (has {fieldnames})"
            )
        for row in reader:
            og_id = row.get("orthogroup_id")
            if not og_id:
                continue
            label = row.get("label", "")
            n_strains = int(row.get("n_strains", "0") or 0)
            og_info[og_id] = (label, n_strains)
            ref_field = row.get(ref_column, "-") or "-"
            if ref_field == "-":
                continue
            for gene in ref_field.split("|"):
                gene = gene.strip()
                if not gene or not gene.startswith(gene_prefix):
                    continue
                og_map[gene] = (og_id, label, n_strains)
    return og_map, og_info


def extract_busted_pvalues(source: str | Path) -> dict[str, float]:
    """Return ``gene_id -> p-value`` from a busted.json directory/tarball."""

    source_path = Path(source)
    if source_path.is_dir():
        return _extract_busted_from_dir(source_path)
    return _extract_busted_from_tar(source_path)


def build_selection_bed_rows(
    busted_results: Mapping[str, float],
    qvals: Mapping[str, float],
    og_map: Mapping[str, tuple[str, str, int]],
    bed12: Mapping[str, list[str]],
    chrom_sizes: Mapping[str, int],
) -> list[str]:
    """Return BED12+5 rows sorted by (chrom, start)."""

    sortable_rows: list[tuple[str, int, str]] = []
    for gene_id, pval in busted_results.items():
        if gene_id not in og_map or gene_id not in bed12:
            continue
        og_id, label, n_strains = og_map[gene_id]
        fields = bed12[gene_id]
        chrom = fields[0]
        if chrom not in chrom_sizes:
            continue
        qval = qvals.get(gene_id, 1.0)
        rgb = qval_to_rgb(qval)
        score = qval_to_score(qval)
        rgb_int = rgb_to_int(rgb)
        gene_family = is_variant_antigen(gene_id, label)
        row = "\t".join(
            [
                fields[0],
                fields[1],
                fields[2],
                og_id,
                str(score),
                fields[5],
                fields[6],
                fields[7],
                str(rgb_int),
                fields[9],
                fields[10],
                fields[11],
                og_id,
                str(n_strains),
                f"{pval:.6g}",
                f"{qval:.6g}",
                gene_family,
            ]
        )
        sortable_rows.append((fields[0], int(fields[1]), row))
    sortable_rows.sort(key=lambda item: (item[0], item[1]))
    return [row for _, _, row in sortable_rows]


def build_orthogroup_bed_rows(
    og_map: Mapping[str, tuple[str, str, int]],
    bed12: Mapping[str, list[str]],
    chrom_sizes: Mapping[str, int],
) -> list[str]:
    """Return BED12 rows for orthogroup membership."""

    sortable_rows: list[tuple[str, int, str]] = []
    for gene_id, (og_id, _label, n_strains) in og_map.items():
        if gene_id not in bed12:
            continue
        fields = bed12[gene_id]
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        if chrom not in chrom_sizes or end > chrom_sizes[chrom] or start >= end:
            continue
        rgb_int = _n_strains_to_rgb_int(n_strains)
        score = int(n_strains * 125)
        row = "\t".join(
            [
                chrom,
                str(start),
                str(end),
                og_id,
                str(score),
                fields[5],
                fields[6],
                fields[7],
                str(rgb_int),
                fields[9],
                fields[10],
                fields[11],
            ]
        )
        sortable_rows.append((chrom, start, row))
    sortable_rows.sort(key=lambda item: (item[0], item[1]))
    return [row for _, _, row in sortable_rows]


def bh_fdr(pvals: Mapping[str, float]) -> dict[str, float]:
    ordered = sorted(pvals.items(), key=lambda item: item[1])
    n = len(ordered)
    qvals: dict[str, float] = {}
    prev = 1.0
    for idx in range(n - 1, -1, -1):
        gene, pval = ordered[idx]
        q = min(prev, pval * n / (idx + 1))
        qvals[gene] = min(1.0, q)
        prev = q
    return qvals


def qval_to_rgb(q: float) -> str:
    if q < 0.01:
        return "255,0,0"
    if q < 0.05:
        return "255,128,0"
    if q < 0.10:
        return "200,200,0"
    return "128,128,128"


def qval_to_score(q: float) -> int:
    if q <= 0:
        return 1000
    return max(0, min(1000, int(-math.log10(q) * 100)))


def rgb_to_int(rgb: str) -> int:
    r, g, b = map(int, rgb.split(","))
    return (r << 16) | (g << 8) | b


def is_variant_antigen(_gene_id: str, label: str) -> str:
    families = {"VIR", "PIR", "PHIST", "DBP", "RBP", "SURFIN", "SERA"}
    lower = label.lower()
    for fam in families:
        if fam.lower() in lower:
            return fam
    return "other"


def _n_strains_to_rgb_int(n_strains: int) -> int:
    ratio = max(0.0, min(1.0, (n_strains - 1) / 7))
    r = max(0, int(255 * (1 - ratio)))
    g = max(0, int(255 * ratio))
    return (r << 16) | g


def _open_maybe_gzip(path: str | Path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def _extract_busted_from_dir(directory: Path) -> dict[str, float]:
    results: dict[str, float] = {}
    for root, _dirs, files in os.walk(directory):
        if "busted.json" not in files:
            continue
        gene_id = Path(root).name
        try:
            with open(Path(root) / "busted.json") as handle:
                data = json.load(handle)
        except Exception:
            continue
        pval = _extract_pvalue(data)
        if pval is not None:
            results[gene_id] = pval
    return results


def _extract_busted_from_tar(archive: Path) -> dict[str, float]:
    results: dict[str, float] = {}
    if not archive.exists():
        return results
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar.getmembers():
            if not member.name.endswith("busted.json"):
                continue
            parts = member.name.split("/")
            if len(parts) < 2:
                continue
            gene_id = parts[-2]
            extracted = tar.extractfile(member)
            if extracted is None:
                continue
            try:
                data = json.load(extracted)
            except Exception:
                continue
            pval = _extract_pvalue(data)
            if pval is not None:
                results[gene_id] = pval
    return results


def _extract_pvalue(payload: dict | None) -> float | None:
    if not payload:
        return None
    try:
        return float(payload.get("test results", {}).get("p-value"))
    except (TypeError, ValueError):
        return None
