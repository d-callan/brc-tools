from pathlib import Path

import pytest

from pangenome_helpers.manifest import (
    ManifestError,
    ensure_matching_collections,
    ensure_reference_not_in_queries,
    load_manifest_map,
)

DATA = Path(__file__).parent / "data" / "manifest"


def test_load_manifest_map_deduplicates():
    mapping = load_manifest_map(DATA / "gff.tsv")
    assert mapping == {"ref": "ref.gff", "strainA": "A.gff", "strainB": "B.gff"}


def test_load_manifest_map_duplicate_raises():
    with pytest.raises(ManifestError):
        load_manifest_map(DATA / "dup.tsv")


def test_ensure_matching_collections_returns_sorted_strains():
    strains, gff_map, fasta_map = ensure_matching_collections(DATA / "gff.tsv", DATA / "fasta.tsv")
    assert strains == ["ref", "strainA", "strainB"]
    assert gff_map["strainA"] == "A.gff"
    assert fasta_map["strainB"] == "B.fa"


def test_ensure_matching_collections_mismatch_raises(tmp_path):
    bad_fa = tmp_path / "fa.tsv"
    bad_fa.write_text("strainA\tA.fa\n")
    with pytest.raises(ManifestError):
        ensure_matching_collections(DATA / "gff.tsv", bad_fa)


def test_ensure_reference_not_in_queries():
    ensure_reference_not_in_queries("ref", ["strainA", "strainB"])
    with pytest.raises(ManifestError):
        ensure_reference_not_in_queries("ref", ["ref", "strainA"])
