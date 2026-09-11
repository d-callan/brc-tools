from pathlib import Path

from pangenome_helpers.consensus import build_consensus_table, summarize_labels

DATA = Path(__file__).parent / "data" / "consensus"


def test_build_consensus_table_core_group(tmp_path):
    liftoff_dir = DATA / "liftoff"
    anchors = ["anchorA"]
    strains = ["anchorA", "strainB"]
    rows = build_consensus_table(
        liftoff_dir,
        anchors,
        strains,
        rbest_edges=[{"strain_a": "anchorA", "gene_a": "anchorGene", "strain_b": "strainB", "gene_b": "queryGene"}],
    )
    assert len(rows) == 1
    row = rows[0]
    assert row["label"] == "CORE-1:1"
    assert row["anchorA"] != "-"
    assert row["strainB"].startswith("queryGene")
    counts = summarize_labels(rows)
    assert counts["CORE-1:1"] == 1


def test_build_consensus_table_with_gene_beds_resolves_projected_to_native():
    """Projected gene (queryGene) should resolve to native gene (nativeGeneB)
    when gene_beds is provided and the coordinates overlap."""
    liftoff_dir = DATA / "liftoff"
    anchors = ["anchorA"]
    strains = ["anchorA", "strainB"]
    rows = build_consensus_table(
        liftoff_dir,
        anchors,
        strains,
        rbest_edges=[{"strain_a": "anchorA", "gene_a": "anchorGene", "strain_b": "strainB", "gene_b": "nativeGeneB"}],
        gene_beds=str(DATA / "gene_beds" / "*.bed"),
    )
    assert len(rows) == 1
    row = rows[0]
    assert row["label"] == "CORE-1:1"
    assert "nativeGeneB" in row["strainB"]
    assert "queryGene" not in row["strainB"]


def test_build_consensus_table_clique_and_density_columns():
    """Clique and density columns should be present and numeric."""
    liftoff_dir = DATA / "liftoff"
    anchors = ["anchorA"]
    strains = ["anchorA", "strainB"]
    rows = build_consensus_table(
        liftoff_dir,
        anchors,
        strains,
        rbest_edges=[{"strain_a": "anchorA", "gene_a": "anchorGene", "strain_b": "strainB", "gene_b": "queryGene"}],
    )
    assert len(rows) == 1
    row = rows[0]
    assert "clique" in row
    assert "density" in row
    assert isinstance(row["clique"], float)
    assert isinstance(row["density"], float)
    assert 0.0 <= row["clique"] <= 1.0
    assert 0.0 <= row["density"] <= 1.0


def test_build_consensus_table_keep_unresolved_projections():
    """When keep_unresolved_projections=True, a projection that doesn't
    resolve to a native gene should still be added as a node."""
    liftoff_dir = DATA / "liftoff"
    anchors = ["anchorA"]
    strains = ["anchorA", "strainB"]
    rows = build_consensus_table(
        liftoff_dir,
        anchors,
        strains,
        gene_beds=str(DATA / "gene_beds" / "*.bed"),
        keep_unresolved_projections=True,
    )
    # queryGene at chr2:0-100 overlaps nativeGeneB at chr2:0-100,
    # so it should be resolved and not kept as a separate node
    assert len(rows) == 1
    row = rows[0]
    assert "nativeGeneB" in row["strainB"]


def test_build_consensus_table_alias_edges_counted_in_density(tmp_path):
    """Alias edges from interval-based aliasing should be counted in used_edges
    so they contribute to density.

    Setup: two projections onto strainB at chr1:0-100 and chr1:0-95 (95%
    reciprocal overlap).  This creates 4 nodes, 3 edges (2 classification + 1
    alias), 6 pairs -> density = 3/6 = 0.5.  Without the fix the alias edge
    is not counted and density = 2/6 = 0.333.
    """
    liftoff_dir = tmp_path / "liftoff"
    (liftoff_dir / "anchorA-as-ref").mkdir(parents=True)
    cls_path = liftoff_dir / "anchorA-as-ref" / "strainB.classification.tsv"
    cls_path.write_text(
        "reference_gene_id\tquery_gene_id\tsource\tintactness\tquery_chrom\tquery_start\tquery_end\n"
        "anchorGene1\tqueryGene1\tcesar2\tI\tchr1\t0\t100\n"
        "anchorGene2\tqueryGene2\tcesar2\tI\tchr1\t0\t95\n"
    )

    rows = build_consensus_table(
        liftoff_dir,
        anchors=["anchorA"],
        strains=["anchorA", "strainB"],
    )
    assert len(rows) == 1
    row = rows[0]
    assert row["density"] == 0.5, f"expected 0.5 (3 edges / 6 pairs), got {row['density']}"
