"""Phase C.2 triage orchestration helpers."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence

from genome_io.gff import parse_gff_attributes_to_dict


STOP_CODONS = {"TAA", "TAG", "TGA"}
CANONICAL_SPLICE = {("GT", "AG"), ("AT", "AC")}
GENE_TYPES = {"gene", "protein_coding_gene", "ncRNA_gene", "pseudogene"}
TRANSCRIPT_TYPES = {"mRNA", "transcript", "pseudogenic_transcript", "tRNA", "rRNA", "ncRNA", "snoRNA", "snRNA"}


@dataclass
class GeneRecord:
    gene_id: str
    reference_id: str
    chrom: str
    start: int
    end: int
    strand: str
    attrs: dict
    transcripts: list


@dataclass
class TriageSettings:
    core_identity_min: float = 0.95
    core_coverage_min: float = 0.90
    family_identity_min: float = 0.85
    subtelomere_bp: int = 100_000


@dataclass
class TriageResult:
    triage_rows: List[dict]
    flagged_reference_ids: set[str]
    clean_gene_ids: set[str]
    needs_cesar2_bed_lines: List[str]
    clean_gff_lines: List[str]
    summary: dict


def parse_liftoff_gff(path: str | Path) -> list[GeneRecord]:
    genes = {}
    tx_to_gene = {}
    tx_features: dict = defaultdict(lambda: {"exon": [], "CDS": []})
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _src, feature, start, end, _score, strand, _phase, attr_str = parts[:9]
            attrs = parse_gff_attributes_to_dict(attr_str)
            if feature in GENE_TYPES:
                gid = attrs.get("ID")
                if gid:
                    genes[gid] = GeneRecord(gid, _normalize_gene_id(gid), chrom, int(start), int(end), strand, attrs, [])
            elif feature in TRANSCRIPT_TYPES:
                tx_id = attrs.get("ID")
                parent = attrs.get("Parent")
                if tx_id and parent:
                    tx_to_gene[tx_id] = parent
            elif feature in ("exon", "CDS"):
                parent = attrs.get("Parent")
                if parent:
                    tx_features[parent][feature].append((int(start), int(end)))
    for tx, gene in tx_to_gene.items():
        if gene in genes:
            exons = sorted(tx_features[tx]["exon"])
            cdss = sorted(tx_features[tx]["CDS"])
            genes[gene].transcripts.append((tx, exons, cdss))
    return list(genes.values())


def _normalize_gene_id(gid: str) -> str:
    if not gid or gid == "None" or "_" not in gid:
        return gid
    core, suffix = gid.rsplit("_", 1)
    if suffix.isdigit() and len(suffix) <= 2 and not core.endswith("_"):
        return core
    return gid


def read_reference_bed(lines: Iterable[str]) -> dict[str, str]:
    mapping = {}
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 4:
            mapping[parts[3]] = line
    return mapping


def read_family_list(lines: Iterable[str]) -> dict[str, str]:
    families = {}
    for ln in lines:
        if not ln.strip() or ln.startswith("#"):
            continue
        parts = ln.rstrip("\n").split("\t")
        if len(parts) >= 2:
            families[parts[0]] = parts[1]
    return families


def run_triage(
    genes: Sequence[GeneRecord],
    fasta_sequences: Mapping[str, str],
    reference_bed_map: Mapping[str, str],
    family_map: Mapping[str, str] | None,
    settings: TriageSettings,
    *,
    liftoff_gff_lines: Iterable[str] | None = None,
) -> TriageResult:
    chrom_sizes = {chrom: len(seq) for chrom, seq in fasta_sequences.items()}
    triage_rows: list[dict] = []
    flagged_ref_ids: set[str] = set()
    clean_gene_ids: set[str] = set()
    rule_counter: dict[str, int] = defaultdict(int)

    for gene in genes:
        triggers, is_family = triage_gene(gene, fasta_sequences, chrom_sizes, family_map or {}, settings)
        if triggers:
            flagged_ref_ids.add(gene.reference_id)
            for t in triggers:
                rule_counter[t.split("_", 1)[0]] += 1
        else:
            clean_gene_ids.add(gene.gene_id)
        triage_rows.append(
            {
                "gene_id": gene.gene_id,
                "reference_id": gene.reference_id,
                "chrom": gene.chrom,
                "start": gene.start,
                "end": gene.end,
                "strand": gene.strand,
                "is_family": is_family,
                "sequence_ID": gene.attrs.get("sequence_ID", ""),
                "coverage": gene.attrs.get("coverage", ""),
                "extra_copy_number": gene.attrs.get("extra_copy_number", "0"),
                "valid_ORFs": gene.attrs.get("valid_ORFs", gene.attrs.get("valid_ORF", "")),
                "decision": "CESAR2_FALLBACK" if triggers else "LIFTOFF_OK",
                "rules_triggered": ",".join(triggers),
            }
        )

    bed_lines = [reference_bed_map[ref] for ref in sorted(flagged_ref_ids) if ref in reference_bed_map]
    clean_gff_lines = list(_filter_clean_gff(liftoff_gff_lines or [], clean_gene_ids))
    summary = {
        "total_genes": len(genes),
        "liftoff_clean": len(clean_gene_ids),
        "needs_cesar2": len(flagged_ref_ids),
        "fallback_rate": len(flagged_ref_ids) / len(genes) if genes else 0,
        "rule_counts": dict(rule_counter),
        "thresholds": settings.__dict__,
    }
    return TriageResult(triage_rows, flagged_ref_ids, clean_gene_ids, bed_lines, clean_gff_lines, summary)


def _filter_clean_gff(lines: Iterable[str], clean_gene_ids: set[str]) -> Iterable[str]:
    current_gene = None
    for line in lines:
        if line.startswith("#"):
            yield line
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        attrs = parse_gff_attributes_to_dict(parts[8])
        if parts[2] in GENE_TYPES:
            current_gene = attrs.get("ID")
        if current_gene in clean_gene_ids:
            yield line


def triage_gene(
    gene: GeneRecord,
    fasta_sequences: Mapping[str, str],
    chrom_sizes: Mapping[str, int],
    family_map: Mapping[str, str],
    settings: TriageSettings,
) -> tuple[list[str], bool]:
    triggers: list[str] = []
    family_membership = family_map.get(gene.reference_id) or family_map.get(gene.gene_id)
    is_family = family_membership is not None
    if is_family:
        triggers.append("R8_family")
    valid_orfs_attr = gene.attrs.get("valid_ORFs", gene.attrs.get("valid_ORF", ""))
    if valid_orfs_attr and valid_orfs_attr.lower() in {"0", "false"}:
        triggers.append("R1a_valid_ORF")
    frame_bad = False
    internal_stop = False
    for _tx_id, _exons, cdss in gene.transcripts:
        if not cdss:
            continue
        cds_nt = _extract_cds_sequence(cdss, gene.chrom, gene.strand, fasta_sequences)
        if len(cds_nt) % 3 != 0:
            frame_bad = True
        if _has_internal_stop(cds_nt):
            internal_stop = True
    if frame_bad:
        triggers.append("R1b_cds_length")
    if internal_stop:
        triggers.append("R1c_internal_stop")
    seq_id = float(gene.attrs.get("sequence_ID", "1.0") or 1.0)
    id_min = settings.family_identity_min if is_family else settings.core_identity_min
    if seq_id < id_min:
        triggers.append("R2_identity")
    coverage = float(gene.attrs.get("coverage", "1.0") or 1.0)
    if coverage < settings.core_coverage_min:
        triggers.append("R3_coverage")
    extra = int(gene.attrs.get("extra_copy_number", "0") or 0)
    if extra > 0 and not is_family:
        triggers.append("R4_extra_copies")
    if gene.attrs.get("partial_mapping", "").lower() == "true":
        triggers.append("R5_partial")
    if _has_splice_issue(gene, fasta_sequences):
        triggers.append("R6_splice")
    if _is_subtelomeric(gene.chrom, gene.start, gene.end, chrom_sizes, settings.subtelomere_bp):
        triggers.append("R7_subtelomeric")
    return triggers, is_family


def _extract_cds_sequence(segments: Sequence[tuple[int, int]], chrom: str, strand: str, fasta_sequences: Mapping[str, str]) -> str:
    seq = fasta_sequences.get(chrom, "")
    if not seq:
        return ""
    parts: list[str] = []
    ordered = sorted(segments, reverse=(strand == "-"))
    for start, end in ordered:
        start_idx = max(0, start - 1)
        end_idx = min(len(seq), end)
        parts.append(seq[start_idx:end_idx])
    cds = "".join(parts)
    return _revcomp(cds) if strand == "-" else cds


def _has_internal_stop(cds_nt: str) -> bool:
    if len(cds_nt) < 6 or len(cds_nt) % 3 != 0:
        return False
    for i in range(0, len(cds_nt) - 3, 3):
        if cds_nt[i:i + 3].upper() in STOP_CODONS:
            return True
    return False


def _has_splice_issue(gene: GeneRecord, fasta_sequences: Mapping[str, str]) -> bool:
    for _tx_id, exons, _cds in gene.transcripts:
        if len(exons) < 2:
            continue
        exons_sorted = sorted(exons, reverse=(gene.strand == "-"))
        for i in range(len(exons_sorted) - 1):
            donor_site = _slice(fasta_sequences, gene.chrom, exons_sorted[i][1] + 1, exons_sorted[i][1] + 2)
            acceptor_site = _slice(fasta_sequences, gene.chrom, exons_sorted[i + 1][0] - 1, exons_sorted[i + 1][0])
            if (donor_site, acceptor_site) not in CANONICAL_SPLICE:
                return True
    return False


def _is_subtelomeric(chrom: str, start: int, end: int, chrom_sizes: Mapping[str, int], flank_bp: int) -> bool:
    size = chrom_sizes.get(chrom)
    if not size:
        return False
    return start < flank_bp or end > size - flank_bp


def _slice(fasta_sequences: Mapping[str, str], chrom: str, start: int, end: int) -> str:
    seq = fasta_sequences.get(chrom, "")
    if not seq:
        return ""
    start_idx = max(0, start - 1)
    end_idx = min(len(seq), end)
    if start_idx >= end_idx:
        return ""
    return seq[start_idx:end_idx].upper()


def _revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]
