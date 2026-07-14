"""Phase C.4 merge helpers (Liftoff-clean + TOGA2)."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List

from genome_io.gff import parse_gff_attributes_to_dict

GENE_TYPES = {"gene", "protein_coding_gene", "ncRNA_gene", "pseudogene"}


@dataclass
class MergeOutputs:
    """Container for merge outputs (classification rows + merged GFF lines)."""

    classification_rows: List[dict]
    merged_gff_lines: List[str]


def load_liftoff_clean(path: str | Path) -> dict[str, list[str]]:
    genes: dict[str, list[str]] = {}
    current_ref = None
    path = Path(path)
    if not path.exists():
        return genes
    with path.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            if feature in GENE_TYPES:
                attrs = parse_gff_attributes_to_dict(fields[8])
                gid = attrs.get("ID", "")
                current_ref = _normalize_gene_id(gid)
                genes.setdefault(current_ref, []).append(line)
            elif current_ref is not None:
                genes[current_ref].append(line)
    return genes


def load_toga_loss_summary(path: str | Path) -> dict[str, str]:
    path = Path(path)
    if not path.exists():
        return {}
    status: dict[str, str] = {}
    with path.open() as handle:
        for row in handle:
            row = row.rstrip("\n")
            if not row or row.startswith("level"):
                continue
            parts = row.split("\t")
            if len(parts) >= 3 and parts[0] == "PROJECTION":
                status[parts[1]] = parts[2]
    return status


def load_toga_orthology(path: str | Path) -> dict[str, list[dict[str, str]]]:
    path = Path(path)
    if not path.exists():
        return {}
    out: dict[str, list[dict[str, str]]] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            ref_gene = row.get("t_gene")
            if not ref_gene:
                continue
            out.setdefault(ref_gene, []).append(
                {
                    "q_gene": row.get("q_gene", ""),
                    "q_tx": row.get("q_transcript", ""),
                    "class": row.get("orthology_class", ""),
                    "t_tx": row.get("t_transcript", ""),
                }
            )
    return out


def load_query_bed(path: str | Path) -> dict[str, str]:
    path = Path(path)
    if not path.exists():
        return {}
    out: dict[str, str] = {}
    with path.open() as handle:
        for row in handle:
            if not row.strip() or row.startswith("#"):
                continue
            fields = row.rstrip("\n").split("\t")
            if len(fields) >= 4:
                out[fields[3]] = row
    return out


def load_reference_genes(path: str | Path) -> set[str]:
    refs: set[str] = set()
    path = Path(path)
    if not path.exists():
        return refs
    with path.open() as handle:
        for row in handle:
            fields = row.rstrip("\n").split("\t")
            if len(fields) >= 4:
                refs.add(fields[3])
    return refs


def merge_annotations(
    query: str,
    liftoff_clean: dict[str, list[str]],
    loss_summary: dict[str, str],
    orthology: dict[str, list[dict[str, str]]],
    query_annotation: dict[str, str],
    query_genes: dict[str, str],
    reference_genes: Iterable[str],
) -> MergeOutputs:
    ref_set = set(reference_genes)
    rows: list[dict] = []
    seen_refs: set[str] = set()

    # Liftoff clean rows
    for ref_id, lines in liftoff_clean.items():
        gene_line = lines[0]
        fields = gene_line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        attrs = parse_gff_attributes_to_dict(fields[8])
        rows.append(
            {
                "reference_gene_id": ref_id,
                "query_gene_id": attrs.get("ID", ""),
                "source": "liftoff",
                "intactness": "I",
                "query_chrom": fields[0],
                "query_start": fields[3],
                "query_end": fields[4],
                "query_strand": fields[6],
                "orthology_class": "liftoff_clean",
            }
        )
        seen_refs.add(ref_id)

    # TOGA projections
    q_to_ref = {}
    for ref_id, projections in orthology.items():
        for proj in projections:
            q_gene = proj.get("q_gene") or ""
            q_to_ref[q_gene] = (
                ref_id,
                proj.get("class", ""),
                proj.get("t_tx", ""),
            )
            if q_gene in ("", "None"):
                rows.append(
                    {
                        "reference_gene_id": ref_id,
                        "query_gene_id": "",
                        "source": "cesar2",
                        "intactness": "L",
                        "query_chrom": "",
                        "query_start": "",
                        "query_end": "",
                        "query_strand": "",
                        "orthology_class": proj.get("class", "one2zero"),
                    }
                )
            else:
                bed_line = query_annotation.get(q_gene) or query_genes.get(q_gene)
                chrom = start = end = strand = ""
                if bed_line:
                    parts = bed_line.rstrip("\n").split("\t")
                    chrom = parts[0] if len(parts) > 0 else ""
                    start = parts[1] if len(parts) > 1 else ""
                    end = parts[2] if len(parts) > 2 else ""
                    strand = parts[5] if len(parts) > 5 else ""
                status = _status_for_projection(loss_summary, proj.get("t_tx", ""))
                rows.append(
                    {
                        "reference_gene_id": ref_id,
                        "query_gene_id": q_gene,
                        "source": "cesar2",
                        "intactness": status,
                        "query_chrom": chrom,
                        "query_start": start,
                        "query_end": end,
                        "query_strand": strand,
                        "orthology_class": proj.get("class", ""),
                    }
                )
            seen_refs.add(ref_id)

    # Missing genes
    for ref_id in sorted(ref_set - seen_refs):
        rows.append(
            {
                "reference_gene_id": ref_id,
                "query_gene_id": "",
                "source": "none",
                "intactness": "M",
                "query_chrom": "",
                "query_start": "",
                "query_end": "",
                "query_strand": "",
                "orthology_class": "unprojected",
            }
        )

    merged_gff = _build_merged_gff_lines(query, liftoff_clean, query_genes, q_to_ref, loss_summary)
    return MergeOutputs(rows, merged_gff)


def _build_merged_gff_lines(
    query: str,
    liftoff_clean: dict[str, list[str]],
    query_genes: dict[str, str],
    q_to_ref: dict[str, tuple[str, str, str]],
    loss_summary: dict[str, str],
) -> list[str]:
    lines = ["##gff-version 3", f"# Phase C.4 merged annotation for {query}"]
    for ref_id, entries in liftoff_clean.items():
        for line in entries:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 9:
                attrs = fields[8].rstrip(";") + ";source=liftoff;intactness=I"
                fields[8] = attrs
                lines.append("\t".join(fields))

    for q_gene, bed_line in query_genes.items():
        parts = bed_line.rstrip("\n").split("\t")
        if len(parts) < 6:
            continue
        chrom, start, end, name, _score, strand = parts[:6]
        ref_gene, orth_class, t_tx = q_to_ref.get(q_gene, ("", "", ""))
        intactness = _status_for_projection(loss_summary, t_tx)
        attrs = (
            f"ID={name};reference_gene_id={ref_gene};source=cesar2;"
            f"intactness={intactness};orthology_class={orth_class}"
        )
        start_1 = str(int(start) + 1) if start.isdigit() else start
        lines.append(
            "\t".join([chrom, "TOGA2", "protein_coding_gene", start_1, end, ".", strand, ".", attrs])
        )
    return lines


def _status_for_projection(loss_summary: dict[str, str], transcript: str) -> str:
    if not transcript:
        return "?"
    prefix = f"{transcript}#"
    for key, value in loss_summary.items():
        if key == transcript or key.startswith(prefix):
            return value
    return "?"


def _normalize_gene_id(gene_id: str) -> str:
    if "_" not in gene_id:
        return gene_id
    core, suffix = gene_id.rsplit("_", 1)
    if suffix.isdigit() and len(suffix) <= 2 and not core.endswith("_"):
        return core
    return gene_id
