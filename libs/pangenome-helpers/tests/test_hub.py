from pathlib import Path

from pangenome_helpers.hub import (
    TrackDbConfig,
    build_genome_records,
    render_genomes_txt,
    render_trackdb,
)

DATA = Path(__file__).parent / "data" / "hub"


def test_render_genomes_txt():
    rows = [{k: v for k, v in zip(
        [
            "accession",
            "defaultPos",
            "organism",
            "scientificName",
            "description",
        ],
        DATA.joinpath("metadata.tsv").read_text().splitlines()[1].split("\t"),
    )}]
    records = build_genome_records(rows)
    text = render_genomes_txt(records)
    assert "genome ACC1" in text
    assert text.endswith("\n")


def test_render_trackdb():
    config = TrackDbConfig(
        assembly="ACC1",
        strain="Strain1",
        species_panel=["ACC1=Strain1", "ACC2=Strain2"],
        anchor_strains=["ACC1=Strain1"],
        maf_url="Strain1.multiz.maf.bb",
        include_selection=True,
    )
    text = render_trackdb(config)
    assert "track Strain1_multiz" in text
    assert "track brc_pangenome_select" in text
