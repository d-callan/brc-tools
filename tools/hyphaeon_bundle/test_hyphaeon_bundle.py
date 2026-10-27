#!/usr/bin/env python3
"""Unit tests for hyphaeon_bundle.py (coordinate math + stats)."""

import importlib.util
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    'hyphaeon_bundle',
    Path(__file__).parent / 'hyphaeon_bundle.py')
hb = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(hb)


# --------------------------------------------------------------------------
# reference coordinate machinery
# --------------------------------------------------------------------------

def test_cds_genomic_positions_plus_single_exon():
    # gene chr1:1000-1180, CDS 1000-1180, + strand, one 180-bp block
    b = ['chr1', '1000', '1180', 'g', '0', '+', '1000', '1180', '0',
         '1', '180,', '0,']
    chrom, strand, pos = hb.cds_genomic_positions(b)
    assert chrom == 'chr1' and strand == '+'
    assert pos == list(range(1000, 1180))


def test_cds_genomic_positions_minus_multiexon():
    # minus-strand gene, 2 exons: chr1:100-160 and chr1:300-330
    # CDS = full exons. On '-' the first CDS codon is at the exon-2 3' end.
    b = ['chr1', '100', '330', 'g', '0', '-', '100', '330', '0',
         '2', '60,30,', '0,200,']
    chrom, strand, pos = hb.cds_genomic_positions(b)
    # '-' strand: transcription starts in the LAST exon (329..300),
    # then continues into the first exon (159..100)
    expected = list(range(329, 299, -1)) + list(range(159, 99, -1))
    assert pos == expected
    assert len(pos) == 90  # 30 codons


def test_cds_positions_respect_thick_bounds():
    # UTRs trimmed: thickStart/thickEnd clip the exon
    b = ['chr1', '1000', '1180', 'g', '0', '+', '1030', '1150', '0',
         '1', '180,', '0,']
    _, _, pos = hb.cds_genomic_positions(b)
    assert pos == list(range(1030, 1150))


def test_aln_to_ref_codon_map_with_gap():
    # ref has a gap codon at alignment position 3 -> maps to nearest flanking
    ref_seq = 'GCT' + 'AAA' + '---' + 'TTT'
    r_of = hb.aln_to_ref_codon_map(ref_seq, 4)
    assert r_of == {1: 1, 2: 2, 3: 2, 4: 3}


def test_codon_span_to_genome():
    # + strand, CDS at chr1:1000-1180 -> codon 5 occupies nt 1012..1014
    b = ['chr1', '1000', '1180', 'g', '0', '+', '1000', '1180', '0',
         '1', '180,', '0,']
    _, _, pos = hb.cds_genomic_positions(b)
    r_of = {i: i for i in range(1, 61)}
    assert hb.codon_span_to_genome(r_of, pos, 5, 5) == (1012, 1015)
    # span codons 10..20 -> 1027..1059+1
    assert hb.codon_span_to_genome(r_of, pos, 10, 20) == (1027, 1060)


def test_codon_span_to_genome_minus_strand():
    # '-' strand, CDS chr1:100-130 single block, 10 codons.
    # codon 1 = rightmost 3 nt (129,128,127) -> [127,130); codon 10 = leftmost.
    b = ['chr1', '100', '130', 'g', '0', '-', '100', '130', '0',
         '1', '30,', '0,']
    _, strand, pos = hb.cds_genomic_positions(b)
    r_of = {i: i for i in range(1, 11)}
    assert hb.codon_span_to_genome(r_of, pos, 1, 1) == (127, 130)
    assert hb.codon_span_to_genome(r_of, pos, 10, 10) == (100, 103)
    # whole-gene span covers the CDS
    assert hb.codon_span_to_genome(r_of, pos, 1, 10) == (100, 130)


def test_normalize_gene_id():
    # Liftoff extra-copy suffixes (_1, _2, ...) collapse; real locus ids
    # with longer numeric tails (or '_' prefixes) are untouched
    assert hb.normalize_gene_id('GENE_2') == 'GENE'
    assert hb.normalize_gene_id('GENE_12') == 'GENE'
    assert hb.normalize_gene_id('GENE_123') == 'GENE_123'
    assert hb.normalize_gene_id('GENE') == 'GENE'
    assert hb.normalize_gene_id('Q8_other') == 'Q8_other'


# --------------------------------------------------------------------------
# reference-gene resolution (set-membership, no prefix)
# --------------------------------------------------------------------------

def test_ref_gene_resolution_composite_cell(tmp_path):
    tsv = tmp_path / 'ortho.tsv'
    tsv.write_text(
        'orthogroup_id\tlabel\tn_strains\tREF\tQ1\n'
        'OG1\tCORE-1:1\t2\tQ9_alias|REF_geneA,Q8_other\tQ1g\n'
        'OG2\tFAMILY\t2\t-\tQ1h\n')
    rows, ref_gene = hb.load_ortholog_table(
        str(tsv), 'REF', {'REF_geneA', 'REF_geneB'})
    assert ref_gene['OG1'] == 'REF_geneA'   # member of ref BED12 names wins
    assert ref_gene['OG2'] == ''
    assert rows['OG1']['label'] == 'CORE-1:1'


def test_ref_gene_resolution_fallback_first_candidate(tmp_path):
    tsv = tmp_path / 'ortho.tsv'
    tsv.write_text(
        'orthogroup_id\tlabel\tn_strains\tREF\tQ1\n'
        'OG1\tCORE-1:1\t2\tNOVELa|NOVELb\tQ1g\n')
    _, ref_gene = hb.load_ortholog_table(str(tsv), 'REF', {'REF_geneA'})
    assert ref_gene['OG1'] == 'NOVELa'      # no member -> first candidate
