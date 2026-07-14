from pathlib import Path

from pangenome_helpers.orthology import Orthogroup, load_intact_orthogroups

DATA = Path(__file__).parent / "data" / "orthology"


def test_load_intact_orthogroups_filters_by_min_strains():
    orthos = load_intact_orthogroups(
        DATA / "ortho.tsv",
        ref_strain="ref",
        strains=["ref", "A", "B"],
        min_intact=2,
    )
    assert orthos == [Orthogroup("OG1", "refGene")]


def test_load_intact_orthogroups_respects_reference_set():
    orthos = load_intact_orthogroups(
        DATA / "ortho.tsv",
        ref_strain="ref",
        strains=["ref", "A"],
        min_intact=1,
        ref_genes={"refGene"},
    )
    assert len(orthos) == 1
    assert orthos[0].reference_gene_id == "refGene"
