"""Manifest helpers specific to pangenome workflows."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Tuple

from genome_io.io import read_manifest


class ManifestError(ValueError):
    """Raised when collection manifests are malformed or inconsistent."""


def load_manifest_map(path: str | Path, *, delimiter: str = "\t") -> Dict[str, str]:
    """Return an identifier -> path mapping from a two-column manifest.

    Wraps :func:`genome_io.io.read_manifest` but enforces uniqueness so downstream
    workflows can rely on deterministic strain ordering.
    """

    entries = read_manifest(path, delimiter=delimiter)
    mapping: Dict[str, str] = {}
    for ident, value in entries:
        if ident in mapping:
            raise ManifestError(f"duplicate identifier {ident!r} in manifest {path}")
        mapping[ident] = value
    return mapping


def ensure_matching_collections(
    gff_manifest: str | Path,
    fasta_manifest: str | Path,
) -> Tuple[list[str], Dict[str, str], Dict[str, str]]:
    """Validate that paired GFF/FASTA manifests list the same set of strains.

    Returns (strain_list, gff_map, fasta_map) where ``strain_list`` is sorted for
    deterministic iteration order.
    """

    gff_map = load_manifest_map(gff_manifest)
    fasta_map = load_manifest_map(fasta_manifest)
    gff_strains = set(gff_map)
    fasta_strains = set(fasta_map)
    if gff_strains != fasta_strains:
        missing_gff = sorted(fasta_strains - gff_strains)
        missing_fasta = sorted(gff_strains - fasta_strains)
        raise ManifestError(
            "GFF/FASTA manifests list different strains: "
            f"missing_gff={missing_gff} missing_fasta={missing_fasta}"
        )
    strains = sorted(gff_map)
    return strains, gff_map, fasta_map


def ensure_reference_not_in_queries(ref: str, queries: Iterable[str]) -> None:
    """Guard against accidently listing the reference strain in query manifests."""

    if ref in set(queries):
        raise ManifestError(f"reference strain {ref!r} also present in queries")
