from pathlib import Path

from pangenome_helpers.cds import OrthogroupSequences, iter_orthogroup_sequences

DATA = Path(__file__).parent / "data" / "cds"


def test_iter_orthogroup_sequences_yields_sequences(tmp_path):
    gff_manifest = tmp_path / "gff_manifest.tsv"
    fasta_manifest = tmp_path / "fasta_manifest.tsv"
    with open(gff_manifest, "w") as fh:
        fh.write(f"strainA\t{DATA / 'strainA.gff3'}\n")
        fh.write(f"strainB\t{DATA / 'strainB.gff3'}\n")
    with open(fasta_manifest, "w") as fh:
        fh.write(f"strainA\t{DATA / 'strainA.fa'}\n")
        fh.write(f"strainB\t{DATA / 'strainB.fa'}\n")

    sequences = list(
        iter_orthogroup_sequences(
            DATA / "ortho.tsv",
            ref_strain="ref",
            ref_gff=DATA / "ref.gff3",
            ref_fasta=DATA / "ref.fa",
            query_gff_manifest=gff_manifest,
            query_fasta_manifest=fasta_manifest,
            min_intact=1,
        )
    )
    assert len(sequences) == 1
    og = sequences[0]
    assert isinstance(og, OrthogroupSequences)
    assert og.orthogroup_id == "OG000001"
    assert og.reference_gene_id == "gene1"
    assert set(og.cds) == {"ref", "strainA", "strainB"}
    assert set(og.proteins) == {"ref", "strainA", "strainB"}
