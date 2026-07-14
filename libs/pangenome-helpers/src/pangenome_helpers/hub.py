"""Assembly hub builders (genomes.txt, trackDb, selection BEDs)."""

from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

DEFAULTPOS_RE = re.compile(r"^[^\s:]+:\d+-\d+$")
FIELD_ORDER = [
    "genome",
    "trackDb",
    "groups",
    "description",
    "twoBitPath",
    "organism",
    "defaultPos",
    "scientificName",
    "htmlPath",
]
REQUIRED_COLUMNS = ["accession", "defaultPos", "organism", "scientificName", "description"]


@dataclass(slots=True)
class GenomeRecord:
    """Information needed for a genomes.txt stanza."""

    accession: str
    trackDb: str
    groups: str
    description: str
    twoBitPath: str
    organism: str
    defaultPos: str
    scientificName: str
    htmlPath: str


def read_metadata_tsv(path: str | Path) -> list[dict[str, str]]:
    """Return metadata rows parsed from a TSV file."""

    with Path(path).open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [row for row in reader]


def build_genome_records(metadata_rows: Iterable[dict[str, str]]) -> list[GenomeRecord]:
    """Validate metadata rows and expand defaults for genomes.txt rendering."""

    records: list[GenomeRecord] = []
    for line_no, row in enumerate(metadata_rows, start=2):
        acc = (row.get("accession") or "").strip()
        if not acc:
            raise ValueError(f"row {line_no}: empty accession")
        default_pos = (row.get("defaultPos") or "").strip()
        if not DEFAULTPOS_RE.match(default_pos):
            raise ValueError(
                f"row {line_no} ({acc}): defaultPos must be chrom:start-end, got {default_pos!r}"
            )

        def col(name: str, default: str) -> str:
            val = (row.get(name) or "").strip()
            return val if val else default

        records.append(
            GenomeRecord(
                accession=acc,
                trackDb=col("trackDb", f"{acc}/trackDb.txt"),
                groups=col("groups", f"{acc}/groups.txt"),
                description=_require(row, "description", line_no, acc),
                twoBitPath=col("twoBitPath", f"{acc}/{acc}.2bit"),
                organism=_require(row, "organism", line_no, acc),
                defaultPos=default_pos,
                scientificName=_require(row, "scientificName", line_no, acc),
                htmlPath=col("htmlPath", f"{acc}/description.html"),
            )
        )
    return records


def render_genomes_txt(records: Sequence[GenomeRecord]) -> str:
    """Convert genome records into UCSC genomes.txt content."""

    if not records:
        raise ValueError("no genomes to render")
    stanzas = [
        "\n".join(f"{field} {_record_value(record, field)}" for field in FIELD_ORDER)
        for record in records
    ]
    return "\n\n".join(stanzas) + "\n"


@dataclass(slots=True)
class TrackDbConfig:
    """Configuration for generating a trackDb.txt section."""

    assembly: str
    strain: str
    species_panel: Sequence[str]
    anchor_strains: Sequence[str]
    maf_url: str
    chains_dir: str = "chains"
    group: str = "brc_pangenome"
    include_selection: bool = False
    html_multiz: str = "../shared/multiz.html"
    html_chains: str = "../shared/chains.html"
    html_annot: str = "../shared/annot.html"
    html_select: str = "../shared/selection.html"


def render_trackdb(config: TrackDbConfig) -> str:
    panel_pairs = [_split_pair(item) for item in config.species_panel]
    anchors = [_split_pair(item) for item in config.anchor_strains]
    targets = [acc for acc, _label in panel_pairs if acc != config.assembly]
    out: list[str] = []
    _emit_maf(out, config, panel_pairs)
    _emit_chains(out, config, targets, dict(panel_pairs))
    _emit_annot(out, config, anchors)
    if config.include_selection:
        _emit_select(out, config)
    return "\n".join(out).rstrip("\n") + "\n"


# ---------------------------------------------------------------------------
# trackDb helpers
# ---------------------------------------------------------------------------


def _split_pair(value: str) -> tuple[str, str]:
    if "=" in value:
        left, right = value.split("=", 1)
        return left.strip(), right.strip()
    return value, value


def _emit_maf(out: list[str], config: TrackDbConfig, pairs: Sequence[tuple[str, str]]):
    name = f"{config.strain.replace('-', '_').replace(' ', '_')}_multiz"
    out.extend(
        [
            f"track {name}",
            f"html {config.html_multiz}",
            f"shortLabel {config.strain} multiz",
            f"longLabel  {config.strain} multi-z alignment",
            "type bigMaf",
            f"bigDataUrl {config.maf_url}",
            f"group {config.group}",
            "visibility pack",
        ]
    )
    if pairs:
        out.append("speciesOrder " + " ".join(acc for acc, _ in pairs))
        labels = " ".join(f'{acc}="{label}"' for acc, label in pairs)
        out.append(f"speciesLabels {labels}")
    out.append("")


def _emit_chains(
    out: list[str], config: TrackDbConfig, targets: Sequence[str], label_map: dict[str, str]
):
    out.extend(
        [
            "track brc_pangenome_chains",
            f"html {config.html_chains}",
            "compositeTrack on",
            "shortLabel Pangenome chains",
            f"longLabel  Pairwise chain alignments ({len(targets)} targets)",
            "type bigChain",
            f"group {config.group}",
            "visibility hide",
            "",
        ]
    )
    for tgt in targets:
        label = label_map.get(tgt, tgt)
        out.extend(
            [
                f"    track chain_to_{tgt}",
                f"    html {config.html_chains}",
                "    parent brc_pangenome_chains off",
                f"    shortLabel chain to {label}",
                f"    longLabel  Chain alignment from {config.assembly} to {label}",
                f"    type bigChain {tgt}",
                f"    bigDataUrl {config.chains_dir}/{config.assembly}_to_{tgt}.bigChain.bb",
                f"    linkDataUrl {config.chains_dir}/{config.assembly}_to_{tgt}.bigChain.link.bb",
                "    visibility hide",
                "",
            ]
        )


def _emit_annot(out: list[str], config: TrackDbConfig, anchors: Sequence[tuple[str, str]]):
    out.extend(
        [
            "track brc_pangenome_annot",
            f"html {config.html_annot}",
            "compositeTrack on",
            "shortLabel Pangenome annot",
            "longLabel  Gene projections (Liftoff + TOGA2)",
            "type bigBed 12",
            f"group {config.group}",
            "visibility pack",
            "",
        ]
    )
    for _acc, label in anchors:
        out.extend(
            [
                f"    track annot_from_{label}",
                f"    html {config.html_annot}",
                "    parent brc_pangenome_annot off",
                f"    shortLabel annot from {label}",
                f"    longLabel  Genes projected from {label}",
                "    type bigBed 12",
                f"    bigDataUrl annot_from_{label}.bb",
                "    visibility dense",
                "",
            ]
        )


def _emit_select(out: list[str], config: TrackDbConfig):
    out.extend(
        [
            "track brc_pangenome_select",
            f"html {config.html_select}",
            "compositeTrack on",
            "shortLabel Pangenome select",
            "longLabel  BUSTED selection + orthogroup membership",
            "type bigBed 12",
            f"group {config.group}",
            "visibility hide",
            "",
        ]
    )
    subs = [
        ("selection_strict", "Selection (strict)", "BUSTED selection, strict core set"),
        ("selection_relaxed", "Selection (relaxed)", "BUSTED selection, relaxed core set"),
        ("orthogroup_membership", "Orthogroups", "Orthogroup membership per gene"),
    ]
    for name, short, long in subs:
        out.extend(
            [
                f"    track {name}",
                f"    html {config.html_select}",
                "    parent brc_pangenome_select off",
                f"    shortLabel {short}",
                f"    longLabel  {long}",
                "    type bigBed 12 +",
                f"    bigDataUrl {name}.bb",
                "    visibility dense",
                "",
            ]
        )


# ---------------------------------------------------------------------------
# Selection helpers shared with build_hub_bb
# ---------------------------------------------------------------------------


def bh_fdr(pvals: dict[str, float]) -> dict[str, float]:
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


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _require(row: dict[str, str], name: str, line_no: int, acc: str) -> str:
    val = (row.get(name) or "").strip()
    if not val:
        raise ValueError(f"row {line_no} ({acc}): missing column {name!r}")
    return val


def _record_value(record: GenomeRecord, field: str) -> str:
    if field == "genome":
        return record.accession
    return getattr(record, field)
