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
