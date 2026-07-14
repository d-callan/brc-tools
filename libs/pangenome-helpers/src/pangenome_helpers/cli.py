"""Command-line interface for pangenome-helpers orchestration functions."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from . import (
    build_consensus_table,
    build_genome_records,
    build_orthogroup_bed_rows,
    build_selection_bed_rows,
    compute_graph_edges,
    compute_rbest_edges,
    derive_multiz_order,
    ensure_matching_collections,
    extract_busted_pvalues,
    iter_orthogroup_sequences,
    load_bed12,
    load_intact_orthogroups,
    load_liftoff_clean,
    load_manifest_map,
    load_ortholog_table,
    load_query_bed,
    load_reference_genes,
    load_sizes,
    load_toga_loss_summary,
    load_toga_orthology,
    merge_annotations,
    orchestrate_phase_c2,
    prepare_anchor_inputs,
    process_maf_file,
    read_family_list,
    read_reference_bed,
    render_genomes_txt,
    render_trackdb,
    rename_fasta,
    run_triage,
    summarize_labels,
    write_phase_c2_outputs,
    TriageSettings,
)


def main(argv: list[str] | None = None) -> int:
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="pangenome-helpers",
        description="Pangenome workflow orchestration helpers",
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Phase C.2: Triage
    triage_parser = subparsers.add_parser(
        "triage",
        help="Phase C.2 triage orchestration",
    )
    triage_parser.add_argument("liftoff_gff", help="Liftoff GFF3 file")
    triage_parser.add_argument("query_fasta", help="Query FASTA file")
    triage_parser.add_argument("reference_bed", help="Reference BED file")
    triage_parser.add_argument("output_dir", help="Output directory")
    triage_parser.add_argument("query_name", help="Query strain name")
    triage_parser.add_argument(
        "--family-list",
        help="Family list TSV (gene_id<TAB>family_name)",
    )
    triage_parser.add_argument(
        "--core-identity-min",
        type=float,
        default=0.95,
        help="Minimum identity for core genes (default: 0.95)",
    )
    triage_parser.add_argument(
        "--core-coverage-min",
        type=float,
        default=0.90,
        help="Minimum coverage for core genes (default: 0.90)",
    )
    triage_parser.add_argument(
        "--family-identity-min",
        type=float,
        default=0.85,
        help="Minimum identity for family genes (default: 0.85)",
    )
    triage_parser.add_argument(
        "--subtelomere-bp",
        type=int,
        default=100_000,
        help="Subtelomeric flank size in bp (default: 100000)",
    )
    triage_parser.set_defaults(func=cmd_triage)

    # Phase C.4: Merge
    merge_parser = subparsers.add_parser(
        "merge",
        help="Phase C.4 merge annotations",
    )
    merge_parser.add_argument("query", help="Query strain name")
    merge_parser.add_argument("liftoff_clean_gff", help="Liftoff clean GFF3")
    merge_parser.add_argument("loss_summary_tsv", help="TOGA loss summary TSV")
    merge_parser.add_argument("orthology_tsv", help="TOGA orthology TSV")
    merge_parser.add_argument("query_annotation_bed", help="Query annotation BED")
    merge_parser.add_argument("query_genes_bed", help="Query genes BED")
    merge_parser.add_argument("reference_bed", help="Reference genes BED")
    merge_parser.add_argument("output_dir", help="Output directory")
    merge_parser.set_defaults(func=cmd_merge)

    # Phase E: Consensus
    consensus_parser = subparsers.add_parser(
        "consensus",
        help="Phase E consensus orthogroup building",
    )
    consensus_parser.add_argument("liftoff_dir", help="Liftoff classification directory")
    consensus_parser.add_argument("anchors", nargs="+", help="Anchor strain names")
    consensus_parser.add_argument("--strains", nargs="+", required=True, help="All strain names")
    consensus_parser.add_argument("--ref-strain", required=True, help="Reference strain name")
    consensus_parser.add_argument(
        "--rbest-edges",
        help="Reciprocal-best edges JSON file",
    )
    consensus_parser.add_argument(
        "--graph-edges",
        help="Graph edges JSON file",
    )
    consensus_parser.add_argument("output_tsv", help="Output consensus TSV")
    consensus_parser.set_defaults(func=cmd_consensus)

    # Phase E: Reciprocal-best edges
    rbest_parser = subparsers.add_parser(
        "rbest-edges",
        help="Phase E reciprocal-best overlap edges",
    )
    rbest_parser.add_argument("chains_pattern", help="Glob pattern for *.chain files")
    rbest_parser.add_argument("annotations_pattern", help="Glob pattern for *.bed files")
    rbest_parser.add_argument(
        "--min-overlap",
        type=float,
        default=0.90,
        help="Minimum fractional overlap (default: 0.90)",
    )
    rbest_parser.add_argument("output_json", help="Output edges JSON")
    rbest_parser.set_defaults(func=cmd_rbest_edges)

    # Phase E: Graph edges
    graph_parser = subparsers.add_parser(
        "graph-edges",
        help="Phase E graph co-membership edges",
    )
    graph_parser.add_argument("paths_tsv", help="ODGI paths TSV")
    graph_parser.add_argument("annotations_pattern", help="Glob pattern for *.bed files")
    graph_parser.add_argument("--strains", nargs="+", required=True, help="Strain names to include")
    graph_parser.add_argument("output_json", help="Output edges JSON")
    graph_parser.set_defaults(func=cmd_graph_edges)

    # Hub building
    hub_parser = subparsers.add_parser(
        "hub",
        help="UCSC hub manifest and track building",
    )
    hub_parser.add_argument("genomes_metadata_tsv", help="Genomes metadata TSV")
    hub_parser.add_argument("output_dir", help="Output directory")
    hub_parser.add_argument(
        "--hub-name",
        help="Hub name (default: from metadata)",
    )
    hub_parser.add_argument(
        "--hub-email",
        help="Hub contact email (default: from metadata)",
    )
    hub_parser.set_defaults(func=cmd_hub)

    # PanSN FASTA renaming
    pansn_parser = subparsers.add_parser(
        "pansn-rename",
        help="Rename FASTA headers with PanSN prefixes",
    )
    pansn_parser.add_argument("input_fasta", help="Input FASTA file")
    pansn_parser.add_argument("output_fasta", help="Output FASTA file")
    pansn_parser.add_argument("sample", help="Sample name")
    pansn_parser.add_argument(
        "--haplotype",
        type=int,
        default=1,
        help="Haplotype number (default: 1)",
    )
    pansn_parser.add_argument(
        "--delimiter",
        default="#",
        help="PanSN delimiter (default: #)",
    )
    pansn_parser.add_argument(
        "--gzip",
        action="store_true",
        help="Gzip output",
    )
    pansn_parser.set_defaults(func=cmd_pansn_rename)

    # MAF processing
    maf_parser = subparsers.add_parser(
        "process-maf",
        help="Filter, reorder, and sort MAF file",
    )
    maf_parser.add_argument("ref_species", help="Reference species name")
    maf_parser.add_argument("input_maf", help="Input MAF file")
    maf_parser.add_argument("output_maf", help="Output MAF file")
    maf_parser.set_defaults(func=cmd_process_maf)

    # Multiz ordering
    multiz_parser = subparsers.add_parser(
        "multiz-order",
        help="Derive multiz query order from sourmash compare matrix",
    )
    multiz_parser.add_argument("compare_csv", help="Sourmash compare CSV")
    multiz_parser.add_argument("hinge", help="Hinge strain name")
    multiz_parser.add_argument("queries", nargs="+", help="Query strain names")
    multiz_parser.add_argument("output_txt", help="Output ordered strain names (one per line)")
    multiz_parser.set_defaults(func=cmd_multiz_order)

    # Anchor prep
    anchor_parser = subparsers.add_parser(
        "anchor-prep",
        help="Prepare anchor inputs (filter BED12, emit isoforms)",
    )
    anchor_parser.add_argument("gff", help="GFF3 file")
    anchor_parser.add_argument("raw_bed", help="Raw gffread BED12")
    anchor_parser.add_argument("output_bed", help="Output filtered BED")
    anchor_parser.add_argument("output_isoforms", help="Output isoforms TSV")
    anchor_parser.set_defaults(func=cmd_anchor_prep)

    # CDS grouping
    cds_parser = subparsers.add_parser(
        "group-cds",
        help="Group CDS by orthogroup",
    )
    cds_parser.add_argument("ortho_table", help="Ortholog table TSV")
    cds_parser.add_argument("ref_strain", help="Reference strain name")
    cds_parser.add_argument("ref_gff", help="Reference GFF3")
    cds_parser.add_argument("ref_fasta", help="Reference FASTA")
    cds_parser.add_argument("query_gff_manifest", help="Query GFF manifest")
    cds_parser.add_argument("query_fasta_manifest", help="Query FASTA manifest")
    cds_parser.add_argument("output_json", help="Output orthogroup sequences JSON")
    cds_parser.add_argument(
        "--min-intact",
        type=int,
        default=2,
        help="Minimum intact strains (default: 2)",
    )
    cds_parser.set_defaults(func=cmd_group_cds)

    args = parser.parse_args(argv)
    if not args.command:
        parser.print_help()
        return 0
    return args.func(args)


def cmd_triage(args: argparse.Namespace) -> int:
    """Phase C.2 triage orchestration."""
    try:
        settings = TriageSettings(
            core_identity_min=args.core_identity_min,
            core_coverage_min=args.core_coverage_min,
            family_identity_min=args.family_identity_min,
            subtelomere_bp=args.subtelomere_bp,
        )
        output = orchestrate_phase_c2(
            liftoff_gff_path=args.liftoff_gff,
            query_fasta_path=args.query_fasta,
            reference_bed_path=args.reference_bed,
            output_dir=args.output_dir,
            query_name=args.query_name,
            family_list_path=args.family_list,
            settings=settings,
        )
        print(f"✓ Triage complete")
        print(f"  TSV: {output.triage_tsv}")
        print(f"  BED: {output.needs_cesar2_bed}")
        print(f"  GFF: {output.liftoff_clean_gff}")
        print(f"  Summary: {output.summary_json}")
        return 0
    except Exception as e:
        print(f"✗ Triage failed: {e}", file=sys.stderr)
        return 1


def cmd_merge(args: argparse.Namespace) -> int:
    """Phase C.4 merge annotations."""
    try:
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        liftoff = load_liftoff_clean(args.liftoff_clean_gff)
        loss = load_toga_loss_summary(args.loss_summary_tsv)
        ortho = load_toga_orthology(args.orthology_tsv)
        query_bed = load_query_bed(args.query_annotation_bed)
        query_genes = load_query_bed(args.query_genes_bed)
        reference = load_reference_genes(args.reference_bed)

        outputs = merge_annotations(
            query=args.query,
            liftoff_clean=liftoff,
            loss_summary=loss,
            orthology=ortho,
            query_annotation=query_bed,
            query_genes=query_genes,
            reference_genes=reference,
        )

        classification_tsv = output_dir / "classification.tsv"
        merged_gff = output_dir / "merged.gff3"

        with open(classification_tsv, "w") as fh:
            if outputs.classification_rows:
                keys = list(outputs.classification_rows[0].keys())
                fh.write("\t".join(keys) + "\n")
                for row in outputs.classification_rows:
                    fh.write("\t".join(str(row.get(k, "")) for k in keys) + "\n")

        with open(merged_gff, "w") as fh:
            for line in outputs.merged_gff_lines:
                fh.write(line + "\n")

        print(f"✓ Merge complete")
        print(f"  Classification: {classification_tsv}")
        print(f"  Merged GFF: {merged_gff}")
        return 0
    except Exception as e:
        print(f"✗ Merge failed: {e}", file=sys.stderr)
        return 1


def cmd_consensus(args: argparse.Namespace) -> int:
    """Phase E consensus orthogroup building."""
    try:
        rbest_edges = None
        if args.rbest_edges:
            with open(args.rbest_edges) as fh:
                rbest_edges = json.load(fh)

        graph_edges = None
        if args.graph_edges:
            with open(args.graph_edges) as fh:
                graph_edges = json.load(fh)

        rows = build_consensus_table(
            liftoff_dir=args.liftoff_dir,
            anchors=args.anchors,
            strains=args.strains,
            ref_strain=args.ref_strain,
            rbest_edges=rbest_edges,
            graph_edges=graph_edges,
        )

        with open(args.output_tsv, "w") as fh:
            if rows:
                keys = list(rows[0].keys())
                fh.write("\t".join(keys) + "\n")
                for row in rows:
                    fh.write("\t".join(str(row.get(k, "")) for k in keys) + "\n")

        counts = summarize_labels(rows)
        print(f"✓ Consensus complete: {len(rows)} orthogroups")
        for label, count in sorted(counts.items()):
            print(f"  {label}: {count}")
        print(f"  Output: {args.output_tsv}")
        return 0
    except Exception as e:
        print(f"✗ Consensus failed: {e}", file=sys.stderr)
        return 1


def cmd_rbest_edges(args: argparse.Namespace) -> int:
    """Phase E reciprocal-best overlap edges."""
    try:
        edges = compute_rbest_edges(
            chains_pattern=args.chains_pattern,
            annotations_pattern=args.annotations_pattern,
            min_overlap=args.min_overlap,
        )
        with open(args.output_json, "w") as fh:
            json.dump(edges, fh, indent=2)
        print(f"✓ Reciprocal-best edges complete: {len(edges)} edges")
        print(f"  Output: {args.output_json}")
        return 0
    except Exception as e:
        print(f"✗ Reciprocal-best edges failed: {e}", file=sys.stderr)
        return 1


def cmd_graph_edges(args: argparse.Namespace) -> int:
    """Phase E graph co-membership edges."""
    try:
        edges = compute_graph_edges(
            paths_tsv=args.paths_tsv,
            annotations_pattern=args.annotations_pattern,
            strains=args.strains,
        )
        with open(args.output_json, "w") as fh:
            json.dump(edges, fh, indent=2)
        print(f"✓ Graph edges complete: {len(edges)} edges")
        print(f"  Output: {args.output_json}")
        return 0
    except Exception as e:
        print(f"✗ Graph edges failed: {e}", file=sys.stderr)
        return 1


def cmd_hub(args: argparse.Namespace) -> int:
    """UCSC hub manifest and track building."""
    try:
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        genomes = build_genome_records(args.genomes_metadata_tsv)
        hub_name = args.hub_name or "Pangenome Hub"
        hub_email = args.hub_email or "contact@example.com"

        genomes_txt = render_genomes_txt(genomes, hub_name, hub_email)
        with open(output_dir / "genomes.txt", "w") as fh:
            fh.write(genomes_txt)

        trackdb_txt = render_trackdb(genomes)
        with open(output_dir / "trackDb.txt", "w") as fh:
            fh.write(trackdb_txt)

        print(f"✓ Hub complete: {len(genomes)} genomes")
        print(f"  genomes.txt: {output_dir / 'genomes.txt'}")
        print(f"  trackDb.txt: {output_dir / 'trackDb.txt'}")
        return 0
    except Exception as e:
        print(f"✗ Hub failed: {e}", file=sys.stderr)
        return 1


def cmd_pansn_rename(args: argparse.Namespace) -> int:
    """Rename FASTA headers with PanSN prefixes."""
    try:
        n = rename_fasta(
            input_path=args.input_fasta,
            output_path=args.output_fasta,
            sample=args.sample,
            haplotype=args.haplotype,
            delimiter=args.delimiter,
            gzip_output=args.gzip,
        )
        print(f"✓ PanSN rename complete: {n} headers")
        print(f"  Output: {args.output_fasta}")
        return 0
    except Exception as e:
        print(f"✗ PanSN rename failed: {e}", file=sys.stderr)
        return 1


def cmd_process_maf(args: argparse.Namespace) -> int:
    """Filter, reorder, and sort MAF file."""
    try:
        result = process_maf_file(
            ref_species=args.ref_species,
            input_path=args.input_maf,
            output_path=args.output_maf,
        )
        print(f"✓ MAF processing complete")
        print(f"  Kept: {result.kept}")
        print(f"  Dropped: {result.dropped}")
        print(f"  Output: {args.output_maf}")
        return 0
    except Exception as e:
        print(f"✗ MAF processing failed: {e}", file=sys.stderr)
        return 1


def cmd_multiz_order(args: argparse.Namespace) -> int:
    """Derive multiz query order from sourmash compare matrix."""
    try:
        ordered = derive_multiz_order(
            compare_csv=args.compare_csv,
            hinge=args.hinge,
            queries=args.queries,
        )
        with open(args.output_txt, "w") as fh:
            for strain in ordered:
                fh.write(strain + "\n")
        print(f"✓ Multiz order complete: {len(ordered)} strains")
        for strain in ordered:
            print(f"  {strain}")
        return 0
    except Exception as e:
        print(f"✗ Multiz order failed: {e}", file=sys.stderr)
        return 1


def cmd_anchor_prep(args: argparse.Namespace) -> int:
    """Prepare anchor inputs."""
    try:
        result = prepare_anchor_inputs(
            gff_path=args.gff,
            raw_bed_path=args.raw_bed,
            out_bed_path=args.output_bed,
            out_isoforms_path=args.output_isoforms,
        )
        print(f"✓ Anchor prep complete")
        print(f"  BED total: {result.bed_total}, kept: {result.bed_kept}")
        print(f"  Isoforms: {result.isoforms}")
        print(f"  Output BED: {args.output_bed}")
        print(f"  Output isoforms: {args.output_isoforms}")
        return 0
    except Exception as e:
        print(f"✗ Anchor prep failed: {e}", file=sys.stderr)
        return 1


def cmd_group_cds(args: argparse.Namespace) -> int:
    """Group CDS by orthogroup."""
    try:
        sequences = list(
            iter_orthogroup_sequences(
                ortho_table=args.ortho_table,
                ref_strain=args.ref_strain,
                ref_gff=args.ref_gff,
                ref_fasta=args.ref_fasta,
                query_gff_manifest=args.query_gff_manifest,
                query_fasta_manifest=args.query_fasta_manifest,
                min_intact=args.min_intact,
            )
        )
        output_data = [
            {
                "orthogroup_id": seq.orthogroup_id,
                "reference_gene_id": seq.reference_gene_id,
                "cds": seq.cds,
                "proteins": seq.proteins,
            }
            for seq in sequences
        ]
        with open(args.output_json, "w") as fh:
            json.dump(output_data, fh, indent=2)
        print(f"✓ CDS grouping complete: {len(sequences)} orthogroups")
        print(f"  Output: {args.output_json}")
        return 0
    except Exception as e:
        print(f"✗ CDS grouping failed: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
