from pathlib import Path

from pangenome_helpers.graph_edges import compute_graph_edges

DATA = Path(__file__).parent / "data" / "graph"


def test_compute_graph_edges_emits_pairs():
    edges = compute_graph_edges(
        DATA / "graph_paths.tsv",
        str(DATA / "*.bed"),
        strains=["sampleA", "sampleB"],
    )
    assert edges
    pair = edges[0]
    assert pair["strain_a"] == "sampleA"
    assert pair["strain_b"] == "sampleB"
    assert pair["path_id"] == "contig1"
