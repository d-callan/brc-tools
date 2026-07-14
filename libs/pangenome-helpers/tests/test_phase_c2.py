from pathlib import Path
from tempfile import TemporaryDirectory

from pangenome_helpers.phase_c2 import orchestrate_phase_c2, write_phase_c2_outputs
from pangenome_helpers.triage import TriageResult, TriageSettings


def test_orchestrate_phase_c2_writes_outputs():
    data_dir = Path(__file__).parent / "data" / "triage"
    with TemporaryDirectory() as tmpdir:
        output = orchestrate_phase_c2(
            liftoff_gff_path=data_dir / "liftoff.gff3",
            query_fasta_path=data_dir / "query.fa",
            reference_bed_path=data_dir / "reference.bed",
            output_dir=tmpdir,
            query_name="test_query",
            family_list_path=data_dir / "family.tsv",
        )
        assert output.triage_tsv.exists()
        assert output.needs_cesar2_bed.exists()
        assert output.liftoff_clean_gff.exists()
        assert output.summary_json.exists()
        summary = output.summary_json.read_text()
        assert "test_query" in summary
        assert "total_genes" in summary


def test_write_phase_c2_outputs():
    result = TriageResult(
        triage_rows=[
            {
                "gene_id": "gene1",
                "reference_id": "gene1",
                "chrom": "chr1",
                "start": 100,
                "end": 200,
                "strand": "+",
                "is_family": False,
                "sequence_ID": "0.95",
                "coverage": "0.90",
                "extra_copy_number": "0",
                "valid_ORFs": "1",
                "decision": "LIFTOFF_OK",
                "rules_triggered": "",
            }
        ],
        flagged_reference_ids=set(),
        clean_gene_ids={"gene1"},
        needs_cesar2_bed_lines=[],
        clean_gff_lines=["##gff-version 3\n"],
        summary={"rule_counts": {}, "thresholds": {}},
    )
    settings = TriageSettings()
    with TemporaryDirectory() as tmpdir:
        output = write_phase_c2_outputs(result, tmpdir, "test", settings)
        assert output.triage_tsv.exists()
        assert output.summary_json.exists()
        triage_content = output.triage_tsv.read_text()
        assert "gene1" in triage_content
