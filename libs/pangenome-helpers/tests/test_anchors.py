from pathlib import Path

from pangenome_helpers.anchors import AnchorPrepResult, prepare_anchor_inputs

DATA = Path(__file__).parent / "data" / "anchors"


def test_prepare_anchor_inputs(tmp_path):
    out_bed = tmp_path / "filtered.bed"
    out_iso = tmp_path / "iso.tsv"
    result = prepare_anchor_inputs(
        DATA / "sample.gff3",
        DATA / "raw.bed",
        out_bed,
        out_iso,
    )
    assert isinstance(result, AnchorPrepResult)
    assert result.bed_total == 2
    assert result.bed_kept == 2
    assert result.isoforms == 2
    assert out_bed.read_text().count("\n") == 2
    assert "geneA" in out_iso.read_text()
