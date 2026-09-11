"""Tests for selection track helpers."""

from pangenome_helpers.selection import (
    _n_strains_to_rgb_int,
    bh_fdr,
    build_orthogroup_bed_rows,
    build_selection_bed_rows,
)


def test_build_selection_bed_rows_name_is_gene_id():
    """BED name field (col 4) should be gene_id, not og_id."""
    bed12 = {
        "gene1": ["chr1", "100", "200", "gene1", "0", "+", "100", "200", "0", "1", "100", "0"],
    }
    chrom_sizes = {"chr1": 1000}
    og_map = {"gene1": ("OG000001", "CORE-1:1", 4)}
    busted = {"gene1": 0.001}
    qvals = bh_fdr(busted)

    rows = build_selection_bed_rows(busted, qvals, og_map, bed12, chrom_sizes)
    assert len(rows) == 1
    fields = rows[0].split("\t")
    assert fields[3] == "gene1", f"name field should be gene_id, got {fields[3]}"
    assert fields[12] == "OG000001", f"extra field 13 should be og_id, got {fields[12]}"


def test_build_orthogroup_bed_rows_name_is_gene_id():
    """Orthogroup BED name field (col 4) should be gene_id, not og_id."""
    bed12 = {
        "gene1": ["chr1", "100", "200", "gene1", "0", "+", "100", "200", "0", "1", "100", "0"],
        "gene2": ["chr1", "300", "400", "gene2", "0", "+", "300", "400", "0", "1", "100", "0"],
    }
    chrom_sizes = {"chr1": 1000}
    og_map = {
        "gene1": ("OG000001", "CORE-1:1", 4),
        "gene2": ("OG000001", "CORE-1:1", 4),
    }

    rows = build_orthogroup_bed_rows(og_map, bed12, chrom_sizes)
    assert len(rows) == 2
    for row in rows:
        fields = row.split("\t")
        gene = fields[3]
        assert gene in ("gene1", "gene2"), f"name field should be gene_id, got {gene}"


def test_n_strains_to_rgb_int_does_not_saturate():
    """With max_strains > 8, colors for 8 and 12 strains should differ."""
    rgb_8 = _n_strains_to_rgb_int(8, max_strains=12)
    rgb_12 = _n_strains_to_rgb_int(12, max_strains=12)
    assert rgb_8 != rgb_12, "colors should differ when max_strains allows differentiation"


def test_n_strains_to_rgb_int_saturates_at_max():
    """At max_strains, the color should be fully green (0 red)."""
    rgb_max = _n_strains_to_rgb_int(12, max_strains=12)
    r = (rgb_max >> 16) & 0xFF
    assert r == 0, f"expected 0 red at max, got {r}"


def test_n_strains_to_rgb_int_backward_compat():
    """Default max_strains=8 preserves original behavior for n_strains <= 8."""
    assert _n_strains_to_rgb_int(1) == _n_strains_to_rgb_int(1, max_strains=8)
    assert _n_strains_to_rgb_int(8) == _n_strains_to_rgb_int(8, max_strains=8)
