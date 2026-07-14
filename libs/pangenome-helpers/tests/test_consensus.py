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
        ref_strain="anchorA",
        rbest_edges=[{"strain_a": "anchorA", "gene_a": "anchorGene", "strain_b": "strainB", "gene_b": "queryGene"}],
    )
    assert len(rows) == 1
    row = rows[0]
    assert row["label"] == "CORE-1:1"
    assert row["anchorA"] != "-"
    assert row["strainB"].startswith("queryGene")
    counts = summarize_labels(rows)
    assert counts["CORE-1:1"] == 1
