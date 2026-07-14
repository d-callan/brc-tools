from pathlib import Path

from pangenome_helpers.overlap import compute_rbest_edges

DATA = Path(__file__).parent / "data" / "overlap"


def test_compute_rbest_edges_simple(tmp_path):
    edges = compute_rbest_edges(
        str(DATA / "*.chain"),
        str(DATA / "*.bed"),
        min_overlap=0.5,
    )
    assert len(edges) == 1
    edge = edges[0]
    assert edge["gene_a"] == "geneA"
    assert edge["gene_b"] == "geneB"
    assert edge["strain_a"] == "strainA"
    assert edge["strain_b"] == "strainB"
