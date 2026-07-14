from pathlib import Path

from pangenome_helpers.triage import (
    TriageSettings,
    parse_liftoff_gff,
    read_family_list,
    read_reference_bed,
    run_triage,
)

DATA = Path(__file__).parent / "data" / "triage"


def test_run_triage_flags_family_gene():
    genes = parse_liftoff_gff(DATA / "liftoff.gff3")
    fasta = {"chr1": (DATA / "query.fa").read_text().splitlines()[1]}
    ref_bed = read_reference_bed((DATA / "reference.bed").read_text().splitlines())
    family = read_family_list((DATA / "family.tsv").read_text().splitlines())

    result = run_triage(genes, fasta, ref_bed, family, TriageSettings(subtelomere_bp=1))
    assert result.summary["total_genes"] == 2
    assert "gene1" in result.clean_gene_ids
    assert "gene2" in result.flagged_reference_ids
    assert any("gene2" in line for line in result.needs_cesar2_bed_lines)
