#!/usr/bin/env python3
"""Reduce a protein FASTA to ONE protein per gene, keeping the longest.

⛔ WHY THIS EXISTS AT ALL: BUSCO's `duplicated` FRACTION IS AN ANNOTATION ARTIFACT UNLESS EVERY
MEMBER OF THE PANEL USES ONE RULE. BUSCO protein mode scores a proteome against single-copy
orthologs, so two isoforms of one gene both hit the same ortholog and the group is called
duplicated. `complete` and `missing` barely move; `duplicated` moves a lot.

⚠ AND `duplicated` IS A NUMBER SOME PANELS ACTUALLY READ. Where a panel holds haplotype-resolved
assemblies, duplication legitimately means the assembly has not collapsed its haplotypes --
workflows/workflow_descriptions.md says so. Mix isoform noise in and you cannot tell annotation
style from assembly biology, which destroys the one signal the number carries there.

⛔ THE RULE IS "LONGEST PER GENE" BECAUSE A PROVIDER'S OWN RULE CANNOT COVER A MIXED PANEL. The
`_protein_primary.faa.gz` files staged before this existed were built by taking RefSeq's DESIGNATED
primary transcript: their record counts equal NCBI's `protein_coding` gene totals exactly, and the
records they discard average LONGER than the ones they keep, so it is demonstrably not the longest
isoform. That designation is a RefSeq annotation attribute -- GenBank submitter annotations and
other publishers do not carry it -- so the rule stops working the moment a proteome arrives from
anywhere else. Longest-per-gene needs nothing but the FASTA, which is why it is the rule that can
cover a panel assembled from several sources.

⚠ HOW THE GENE IS IDENTIFIED IS PER-SOURCE, AND GUESSING IT IS THE WAY TO GET THIS SILENTLY WRONG.
A pattern that fails to group isoforms leaves the file unchanged and reports success, which looks
exactly like a proteome that was already one-per-gene. So `--gene-regex` is REQUIRED, and the
summary prints how many genes had more than one transcript: if that is 0 on a file you expected to
reduce, the pattern did not match, not the data.

    # transcript-suffixed ids, e.g. Xx_1.10g000010.m01.polypeptide -> gene Xx_1.10g000010
    python3 scripts/primary_proteome.py in.fasta --out out.faa.gz \\
        --gene-regex '^(?P<gene>[^ ]+?)\\.m[0-9]+\\.polypeptide'

    # NCBI style: XP_030478000.1 with the gene in the FASTA description as [gene=NAME]
    python3 scripts/primary_proteome.py in.faa.gz --out out.faa.gz \\
        --gene-regex 'gene=(?P<gene>[^]]+)'
"""

from __future__ import annotations

import argparse
import gzip
import pathlib
import re
import sys


def opener(path: pathlib.Path):
    """gzip or plain, decided by the bytes rather than by the suffix.

    ⚠ NOT `path.suffix == ".gz"`. The panel's files arrive from three sources with inconsistent
    naming -- GigaDB ships `.fasta` uncompressed, NCBI ships `.faa.gz` -- and a suffix test would
    hand gzip bytes to the text reader as mojibake rather than failing.
    """
    with path.open("rb") as fh:
        magic = fh.read(2)
    return gzip.open if magic == b"\x1f\x8b" else open


def records(path: pathlib.Path):
    """Yield (header_without_'>', sequence) for a FASTA, streaming."""
    fn = opener(path)
    name, seq = None, []
    with fn(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].rstrip("\n"), []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", type=pathlib.Path)
    ap.add_argument("--out", type=pathlib.Path, required=True,
                    help="written gzipped if the name ends .gz")
    ap.add_argument("--gene-regex", required=True,
                    help="regex with a named group `gene`, matched against the FASTA header. "
                         "⚠ REQUIRED: a pattern that matches nothing leaves the file unchanged "
                         "and looks like success.")
    a = ap.parse_args()

    try:
        pat = re.compile(a.gene_regex)
    except re.error as e:
        sys.exit(f"--gene-regex is not a valid regex: {e}")
    if "gene" not in (pat.groupindex or {}):
        sys.exit(f"--gene-regex {a.gene_regex!r} has no named group `gene`. Write it as "
                 f"`(?P<gene>...)` so this knows which part identifies the gene.")

    best: dict[str, tuple[int, str, str]] = {}
    total = unmatched = 0
    for header, seq in records(a.fasta):
        total += 1
        m = pat.search(header)
        if not m:
            unmatched += 1
            # ⛔ AN UNMATCHED HEADER IS ITS OWN GENE, NOT A DISCARD. Dropping it would silently
            # shorten the proteome, and a proteome that is quietly short is the failure this whole
            # exercise is about. It is counted and reported instead.
            gene = header.split()[0]
        else:
            gene = m.group("gene")
        prev = best.get(gene)
        if prev is None or len(seq) > prev[0]:
            best[gene] = (len(seq), header, seq)

    if unmatched:
        print(f"  ⚠ {unmatched} of {total} header(s) did not match --gene-regex and were each "
              f"treated as their own gene. If that is most of the file, the pattern is wrong.")

    groups = total - len(best)
    print(f"  {total} record(s) -> {len(best)} gene(s); {groups} record(s) dropped as shorter "
          f"isoforms")
    if groups == 0:
        # Not an error: plenty of published proteomes are natively one-per-gene. But it is
        # indistinguishable from a pattern that matched nothing, so say which it was.
        print(f"  ⚠ NOTHING WAS REDUCED. Either this proteome is already one protein per gene, or "
              f"--gene-regex grouped nothing ({total - unmatched} header(s) did match). Check "
              f"before assuming the former.")

    out = a.out
    fn = gzip.open if out.name.endswith(".gz") else open
    with fn(out, "wt") as fh:
        # Sorted by gene so the file is byte-reproducible across runs; a proteome that differs
        # between two runs of the same input cannot be checksummed into a provenance record.
        for gene in sorted(best):
            _, header, seq = best[gene]
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    print(f"  wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
