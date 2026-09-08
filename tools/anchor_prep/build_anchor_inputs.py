#!/usr/bin/env python3
"""Build per-anchor BED12 + isoforms TSV for TOGA2 input.

Ported from
/home/anton/git/Pv4-pangenome/v3/pipeline/impl/setup/build_anchor_inputs.sh
(the gffread BED12 -> protein-coding filter, plus the mRNA ID/Parent
isoforms extraction). The gffread --bed call itself is run by the Galaxy
command block; this script consumes the raw BED12 plus the source GFF3 and
emits the filtered BED12 and the isoforms TSV.
"""
import argparse
import pathlib
import re


def collect_protein_coding_genes(gff_path):
    """Gene IDs of protein_coding_gene / gene features (ID= from column 9)."""
    pc = set()
    with pathlib.Path(gff_path).open() as fin:
        for ln in fin:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] not in ("protein_coding_gene", "gene"):
                continue
            for attr in f[8].split(";"):
                if attr.startswith("ID="):
                    pc.add(attr[3:])
                    break
    return pc


def filter_bed12(raw_bed_path, out_path, pc):
    """Keep BED12 rows whose gene maps to a protein-coding gene.

    gffread --bed emits 12 BED columns plus extra attribute columns; the
    name column (column 4) holds the *transcript* id, while the gene id
    lives in the trailing attribute column as geneID=....

    phase_c4_merge keys its reference gene set off the BED name column
    (fields[3]) and reconciles it against the clean-GFF *gene* IDs, so the
    BED name MUST be the GENE id (not the transcript id) or every gene
    shows up as a spurious unprojected `.t1` row (Finding F3). We therefore
    rewrite column 4 to the gene id before emitting the 12-column BED.
    """
    n_in = n_kept = 0
    with pathlib.Path(raw_bed_path).open() as fin, pathlib.Path(out_path).open("w") as fout:
        for ln in fin:
            ln = ln.rstrip("\n")
            if not ln:
                continue
            f = ln.split("\t")
            n_in += 1
            if len(f) < 12:
                continue
            gid = None
            if len(f) >= 13:
                m = re.search(r"geneID=([^;]+)", f[12])
                if m:
                    gid = m.group(1)
            if gid is None:
                # No geneID= attribute; fall back to the existing BED name.
                gid = f[3]
            if not pc or gid in pc:
                out_fields = f[:12]
                out_fields[3] = gid  # gene id in the name column (F3)
                fout.write("\t".join(out_fields) + "\n")
                n_kept += 1
    return n_in, n_kept


def build_isoforms(gff_path, out_path):
    """Isoforms TSV: gene_id<TAB>transcript_id from mRNA features.

    transcript_id = ID= of the mRNA, gene_id = Parent= of the mRNA.
    """
    rows = []
    with pathlib.Path(gff_path).open() as fin:
        for ln in fin:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "mRNA":
                continue
            tx = gene = ""
            for attr in f[8].split(";"):
                if attr.startswith("ID="):
                    tx = attr[3:]
                elif attr.startswith("Parent="):
                    gene = attr[7:]
            if tx and gene:
                rows.append((gene, tx))
    rows.sort(key=lambda r: (r[0], r[1]))
    with pathlib.Path(out_path).open("w") as fout:
        for gene, tx in rows:
            fout.write(gene + "\t" + tx + "\n")
    return len(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gff", required=True, help="Source GFF3 annotation")
    ap.add_argument("--raw-bed", required=True, help="Raw BED12 from gffread --bed")
    ap.add_argument("--out-bed", required=True, help="Filtered protein-coding BED12 output")
    ap.add_argument("--out-isoforms", required=True, help="Isoforms TSV output")
    args = ap.parse_args()

    pc = collect_protein_coding_genes(args.gff)
    n_in, n_kept = filter_bed12(args.raw_bed, args.out_bed, pc)
    n_iso = build_isoforms(args.gff, args.out_isoforms)
    print(f"bed12: kept {n_kept}/{n_in}  isoforms: {n_iso}")


if __name__ == "__main__":
    main()
