from pathlib import Path

from pangenome_helpers.merge import (
    MergeOutputs,
    load_liftoff_clean,
    load_query_bed,
    load_reference_genes,
    load_toga_loss_summary,
    load_toga_orthology,
    merge_annotations,
)

DATA = Path(__file__).parent / "data" / "merge"


def test_merge_annotations_combines_sources():
    liftoff = load_liftoff_clean(DATA / "liftoff_clean.gff3")
    loss = load_toga_loss_summary(DATA / "loss_summary.tsv")
    ortho = load_toga_orthology(DATA / "orthology.tsv")
    query_bed = load_query_bed(DATA / "query_annotation.bed")
    query_genes = load_query_bed(DATA / "query_genes.bed")
    reference = load_reference_genes(DATA / "ref.bed")

    outputs = merge_annotations(
        query="Q",
        liftoff_clean=liftoff,
        loss_summary=loss,
        orthology=ortho,
        query_annotation=query_bed,
        query_genes=query_genes,
        reference_genes=reference,
    )
    assert isinstance(outputs, MergeOutputs)
    assert len(outputs.classification_rows) == 2  # liftoff + TOGA record
    assert outputs.classification_rows[0]["source"] == "liftoff"
    merged_gff = "\n".join(outputs.merged_gff_lines)
    assert "source=cesar2" in merged_gff
