"""Phase C.2 triage orchestration and reporting."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

from .triage import TriageResult, TriageSettings, parse_liftoff_gff, read_family_list, read_reference_bed, run_triage


@dataclass(slots=True)
class Phase2Output:
    """Paths to Phase C.2 output files."""

    triage_tsv: Path
    needs_cesar2_bed: Path
    liftoff_clean_gff: Path
    summary_json: Path


def write_phase_c2_outputs(
    result: TriageResult,
    output_dir: str | Path,
    query_name: str,
    settings: TriageSettings,
) -> Phase2Output:
    """Write triage results to TSV, BED, GFF, and summary JSON files."""

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    triage_tsv = output_dir / "triage.tsv"
    needs_cesar2_bed = output_dir / "needs_cesar2.bed"
    liftoff_clean_gff = output_dir / "liftoff_clean.gff3"
    summary_json = output_dir / "summary.json"

    _write_triage_tsv(triage_tsv, result.triage_rows)
    _write_needs_cesar2_bed(needs_cesar2_bed, result.needs_cesar2_bed_lines)
    _write_liftoff_clean_gff(liftoff_clean_gff, result.clean_gff_lines)
    _write_summary_json(summary_json, result, query_name, settings)

    return Phase2Output(triage_tsv, needs_cesar2_bed, liftoff_clean_gff, summary_json)


def orchestrate_phase_c2(
    liftoff_gff_path: str | Path,
    query_fasta_path: str | Path,
    reference_bed_path: str | Path,
    output_dir: str | Path,
    query_name: str,
    family_list_path: str | Path | None = None,
    settings: TriageSettings | None = None,
) -> Phase2Output:
    """Full Phase C.2 pipeline: parse, triage, and write outputs."""

    if settings is None:
        settings = TriageSettings()

    genes = parse_liftoff_gff(liftoff_gff_path)
    fasta_sequences = _load_fasta_sequences(query_fasta_path)
    reference_bed_map = read_reference_bed(_read_lines(reference_bed_path))
    family_map = read_family_list(_read_lines(family_list_path)) if family_list_path else {}
    liftoff_lines = _read_lines(liftoff_gff_path)

    result = run_triage(
        genes,
        fasta_sequences,
        reference_bed_map,
        family_map,
        settings,
        liftoff_gff_lines=liftoff_lines,
    )

    return write_phase_c2_outputs(result, output_dir, query_name, settings)


def _write_triage_tsv(path: Path, rows: list[dict]) -> int:
    """Write triage decisions to TSV. Returns number of rows written."""

    if not rows:
        return 0
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def _write_needs_cesar2_bed(path: Path, lines: list[str]) -> int:
    """Write flagged genes BED file. Returns number of lines written."""

    with open(path, "w") as handle:
        for line in lines:
            handle.write(line)
    return len(lines)


def _write_liftoff_clean_gff(path: Path, lines: list[str]) -> int:
    """Write filtered GFF3 file. Returns number of lines written."""

    with open(path, "w") as handle:
        for line in lines:
            handle.write(line)
    return len(lines)


def _write_summary_json(
    path: Path,
    result: TriageResult,
    query_name: str,
    settings: TriageSettings,
) -> None:
    """Write summary statistics and thresholds to JSON."""

    summary = {
        "query": query_name,
        "total_genes": len(result.triage_rows),
        "liftoff_clean": len(result.clean_gene_ids),
        "needs_cesar2": len(result.flagged_reference_ids),
        "needs_cesar2_in_bed": len(result.needs_cesar2_bed_lines),
        "fallback_rate": len(result.flagged_reference_ids) / len(result.triage_rows) if result.triage_rows else 0.0,
        "rule_counts": result.summary.get("rule_counts", {}),
        "thresholds": asdict(settings),
    }
    with open(path, "w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)


def _load_fasta_sequences(path: str | Path) -> dict[str, str]:
    """Load FASTA sequences into memory as a dict: chrom -> sequence."""

    sequences: dict[str, str] = {}
    current_chrom = None
    current_seq = []

    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_chrom is not None:
                    sequences[current_chrom] = "".join(current_seq).upper()
                current_chrom = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

    if current_chrom is not None:
        sequences[current_chrom] = "".join(current_seq).upper()

    return sequences


def _read_lines(path: str | Path | None) -> Iterable[str]:
    """Read lines from a file, yielding each line."""

    if path is None:
        return
    with open(path) as handle:
        for line in handle:
            yield line
