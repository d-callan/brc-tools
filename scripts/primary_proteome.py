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

⛔ THE RULE IS "LONGEST PER GENE" BECAUSE A PROVIDER'S OWN RULE CANNOT COVER A MIXED PANEL, NOT
BECAUSE THE PROVIDER'S RULE IS WORSE. RefSeq designates a primary transcript per gene; GenBank
submitter annotations and other publishers do not carry that attribute at all, so a pipeline built
on it stops working the moment a proteome arrives from anywhere else. Longest-per-gene needs
nothing but the FASTA and a gene assignment, which is why it is the rule that can cover a panel
assembled from several sources.

⚠ AND ON THIS PANEL IT SELECTS THE SAME PROTEIN REFSEQ ALREADY DID -- MEASURED, both RefSeq members,
all 54,043 genes, IDENTICAL ACCESSIONS AND SEQUENCES. Recording that because the opposite was
written down here first, on this evidence: the discarded records average LONGER than the kept ones
(cs10 513 aa vs 438; ASM2916894v1 525 vs 440). That is true and it is not evidence. It pools across
DIFFERENT GENE SETS -- a gene with several isoforms is a longer gene -- so the pooled order reverses
the within-gene order, and within every single gene the kept record is maximal. A pooled mean cannot
answer a per-group question. `--check` is the test that actually answers it, which is why it exists.

⚠ HOW THE GENE IS IDENTIFIED IS PER-SOURCE, AND GUESSING IT IS THE WAY TO GET THIS SILENTLY WRONG.
A pattern that fails to group isoforms leaves the file unchanged and reports success, which looks
exactly like a proteome that was already one-per-gene. So one of `--gene-regex` / `--gene-map-gff3`
is REQUIRED, and the summary prints how many records were dropped: if that is 0 on a file you
expected to reduce, the grouping did not work, not the data.

    # transcript-suffixed ids, e.g. Xx_1.10g000010.m01.polypeptide -> gene Xx_1.10g000010
    python3 scripts/primary_proteome.py in.fasta --out out.faa.gz \\
        --gene-regex '^(?P<gene>[^ ]+?)\\.m[0-9]+\\.polypeptide'

    # NCBI/RefSeq: the header carries NO gene, so the gene comes from the annotation
    python3 scripts/primary_proteome.py GCF_x_protein.faa.gz --out out.faa.gz \\
        --gene-map-gff3 GCF_x_genomic.gff.gz

    # audit a proteome someone else reduced -- writes nothing, exits non-zero if it is not
    # one-longest-protein-per-gene
    python3 scripts/primary_proteome.py GCF_x_protein.faa.gz \\
        --check GCF_x_protein_primary.faa.gz --gene-map-gff3 GCF_x_genomic.gff.gz

    python3 scripts/primary_proteome.py --self-test    # prove each guard fires, with its message

⛔ AND FOR AN NCBI PROTEOME THE REGEX ROUTE DOES NOT EXIST -- IT IS NOT A MATTER OF FINDING THE
RIGHT PATTERN. An FTP `*_protein.faa.gz` header is an accession plus a PRODUCT NAME:

    >XP_030477617.1 rust resistance kinase Lr10 isoform X1 [Cannabis sativa]

Zero of cs10's 33,674 headers and zero of ASM2916894v1's 39,959 contain `[gene=`; that bracketed
form belongs to `*_translated_cds.faa`, a different file. Grouping on the product name would merge
unrelated paralogs that share a name and split one gene whose isoforms are described differently,
so the gene has to come from the GFF3, where the parent chain states it. `--gene-map-gff3` walks
`CDS.protein_id -> Parent -> ... -> gene` using the GFF3's own `ID`/`Parent` links rather than any
naming convention, since RefSeq's `rna-`/`gene-` prefixes are a convention and submitter GFF3s do
not share it.

⚠ AN UNRESOLVED ACCESSION IS A REFUSAL IN THIS MODE, unlike the regex route where an unmatched
header is its own gene. A regex can legitimately not apply to some headers; a GFF3 that does not
mention a protein the FASTA contains means the two files are from DIFFERENT annotation releases,
and grouping the remainder would produce a proteome that looks complete and is silently mixed.
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


GENE_TYPES = frozenset({"gene", "pseudogene"})


def attrs_of(col9: str) -> dict[str, str]:
    """Split a GFF3 attribute column. Only the keys this needs, no URL-decoding."""
    out = {}
    for field in col9.rstrip().split(";"):
        key, sep, val = field.partition("=")
        if sep:
            out[key] = val
    return out


def gene_map_from_gff3(path: pathlib.Path) -> dict[str, str]:
    """protein accession -> gene id, from the GFF3's own ID/Parent chain.

    ⚠ WALKED, NOT PATTERN-MATCHED. A RefSeq CDS's parent is `rna-XM_...` whose parent is
    `gene-LOC...`, but those prefixes are a RefSeq convention and a submitter GFF3 uses its own
    ids, so the walk stops on a feature whose TYPE is gene-like rather than on a name.
    """
    parent_of: dict[str, str] = {}
    type_of: dict[str, str] = {}
    prot_parent: dict[str, str] = {}
    attr_gene: dict[str, str] = {}
    fn = opener(path)
    with fn(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            col = line.split("\t")
            if len(col) < 9:
                continue
            at = attrs_of(col[8])
            fid = at.get("ID")
            par = at.get("Parent", "").split(",")[0]
            if fid and fid not in type_of:
                # setdefault semantics: a multi-row CDS repeats its ID, first row wins.
                type_of[fid] = col[2]
                parent_of[fid] = par
            if col[2] == "CDS" and "protein_id" in at:
                prot_parent.setdefault(at["protein_id"], par)
                if "gene" in at:
                    attr_gene.setdefault(at["protein_id"], at["gene"])

    if not prot_parent:
        sys.exit(f"{path}: no CDS feature carries a `protein_id` attribute, so no protein can be "
                 f"mapped to a gene. Is this a GFF3 with annotation, or only an assembly report?")

    out: dict[str, str] = {}
    orphans: list[str] = []
    for prot, par in prot_parent.items():
        cur, seen = par, set()
        while cur and cur not in seen:
            seen.add(cur)
            if type_of.get(cur) in GENE_TYPES:
                out[prot] = cur
                break
            cur = parent_of.get(cur, "")
        else:
            orphans.append(prot)

    if orphans:
        # ⛔ NOT "group it by its transcript instead". That is precisely the mis-grouping this mode
        # exists to prevent: every isoform becomes its own gene and nothing is reduced.
        sys.exit(f"{path}: {len(orphans)} protein(s) have no gene feature anywhere up their "
                 f"Parent chain, e.g. {', '.join(sorted(orphans)[:5])}. Grouping them by "
                 f"transcript would silently keep every isoform.")

    genes = len(set(out.values()))
    print(f"  gene map: {len(out)} protein(s) -> {genes} gene(s) from {path.name}")
    _cross_check(path, out, attr_gene)
    return out


def _cross_check(path: pathlib.Path, chain: dict[str, str], attr: dict[str, str]) -> None:
    """Second, independent route to the same grouping: the CDS line's own `gene=` attribute.

    ⚠ THE COMPARISON IS DELIBERATELY ONE-SIDED, and the asymmetry is the whole content of it.
    `gene=` carries a SYMBOL, and two distinct gene features can legitimately share one -- paralogs
    annotated with the same name -- so the attribute route MERGING two chain-genes is explainable
    and is only reported. The attribute route SPLITTING one chain-gene is not explainable by
    anything in the format: the same gene feature cannot have two symbols. That direction means one
    of the two readings is wrong, and since they are the only two the file offers, there is nothing
    left to break the tie -- so it refuses rather than pick.
    """
    if not attr:
        return
    shared = [p for p in chain if p in attr]
    if not shared:
        return
    symbols: dict[str, set[str]] = {}
    for prot in shared:
        symbols.setdefault(chain[prot], set()).add(attr[prot])
    split = {g: v for g, v in symbols.items() if len(v) > 1}
    if split:
        ex = "; ".join(f"{g} -> {sorted(v)}" for g, v in sorted(split.items())[:3])
        sys.exit(f"{path}: {len(split)} gene feature(s) whose CDS rows disagree with the Parent "
                 f"chain about which gene they belong to, e.g. {ex}. The two routes this file "
                 f"offers contradict each other, so neither can be trusted for grouping.")
    per_symbol: dict[str, set[str]] = {}
    for gene, syms in symbols.items():
        per_symbol.setdefault(next(iter(syms)), set()).add(gene)
    merged = sum(len(v) - 1 for v in per_symbol.values() if len(v) > 1)
    note = f"; {merged} gene(s) share a symbol with another" if merged else ""
    print(f"  cross-check: the CDS `gene=` attribute agrees with the Parent chain on all "
          f"{len(shared)} protein(s){note}")


def check_reduced(source: pathlib.Path, reduced: pathlib.Path,
                  gmap: dict[str, str]) -> int:
    """Assert `reduced` holds exactly one protein per gene, and the longest one.

    ⛔ THIS IS THE TEST A POOLED MEAN CANNOT BE. "The discarded records are longer on average" is a
    statement about two different gene sets; "this gene's kept record is its longest" is a statement
    about one gene, asserted here for every gene separately. See the module docstring for the
    occasion on which the difference mattered.
    """
    src = {}
    for header, seq in records(source):
        src[header.split()[0]] = seq
    kept = {}
    for header, seq in records(reduced):
        kept[header.split()[0]] = seq

    problems: list[str] = []

    strangers = sorted(set(kept) - set(src))
    if strangers:
        problems.append(f"{len(strangers)} record(s) in {reduced.name} are not in the source, "
                        f"e.g. {', '.join(strangers[:5])}")
    altered = sorted(k for k in kept if k in src and kept[k] != src[k])
    if altered:
        problems.append(f"{len(altered)} record(s) have a different sequence from the source, "
                        f"e.g. {', '.join(altered[:5])} -- a reduction must not edit sequences")

    by_gene: dict[str, list[str]] = {}
    for prot, gene in gmap.items():
        if prot in src:
            by_gene.setdefault(gene, []).append(prot)

    wrong_count, not_longest = [], []
    for gene, prots in by_gene.items():
        here = [p for p in prots if p in kept]
        if len(here) != 1:
            wrong_count.append(f"{gene} ({len(here)})")
        elif len(src[here[0]]) != max(len(src[p]) for p in prots):
            not_longest.append(f"{gene} kept {here[0]} at {len(src[here[0]])} aa, longest is "
                               f"{max(len(src[p]) for p in prots)} aa")
    if wrong_count:
        problems.append(f"{len(wrong_count)} gene(s) do not contribute exactly one record, "
                        f"e.g. {'; '.join(sorted(wrong_count)[:5])}")
    if not_longest:
        problems.append(f"{len(not_longest)} gene(s) kept a shorter isoform, "
                        f"e.g. {'; '.join(sorted(not_longest)[:3])}")

    if not by_gene:
        # ⛔ 0 == 0 IS NOT A PASS. An empty gene map compares nothing and would report success.
        problems.append(f"the gene map named none of {source.name}'s {len(src)} protein(s), so "
                        f"nothing was actually checked")

    if problems:
        print(f"❌ {reduced.name} is not one-longest-protein-per-gene:")
        for line in problems:
            print(f"  - {line}")
        return 1
    print(f"✅ {reduced.name}: {len(kept)} record(s) = {len(by_gene)} gene(s), each the longest of "
          f"its gene, sequences unchanged from {source.name}")
    return 0


# A GFF3 with two genes: gene-G1 has two mRNAs (P1 400 aa, P2 900 aa), gene-G2 has one (P3 700 aa).
# CDS rows carry `gene=` so the two-route cross-check has something to compare.
_GFF3_OK = """\
chr1\tx\tgene\t1\t900\t.\t+\t.\tID=gene-G1
chr1\tx\tmRNA\t1\t900\t.\t+\t.\tID=rna-T1;Parent=gene-G1
chr1\tx\tCDS\t1\t300\t.\t+\t0\tID=cds-P1;Parent=rna-T1;gene=SYM_1;protein_id=P1
chr1\tx\tCDS\t400\t900\t.\t+\t0\tID=cds-P1;Parent=rna-T1;gene=SYM_1;protein_id=P1
chr1\tx\tmRNA\t1\t900\t.\t+\t.\tID=rna-T2;Parent=gene-G1
chr1\tx\tCDS\t1\t900\t.\t+\t0\tID=cds-P2;Parent=rna-T2;gene=SYM_1;protein_id=P2
chr1\tx\tgene\t2000\t2600\t.\t+\t.\tID=gene-G2
chr1\tx\tmRNA\t2000\t2600\t.\t+\t.\tID=rna-T3;Parent=gene-G2
chr1\tx\tCDS\t2000\t2600\t.\t+\t0\tID=cds-P3;Parent=rna-T3;gene=SYM_2;protein_id=P3
"""

_FASTA = {
    # the source: gene-G1 twice, gene-G2 once
    "src.faa": ">P1 short isoform\n" + "A" * 400 + "\n>P2 long isoform\n" + "A" * 900
               + "\n>P3 sole isoform\n" + "A" * 700 + "\n",
    "good.faa": ">P2 long isoform\n" + "A" * 900 + "\n>P3 sole isoform\n" + "A" * 700 + "\n",
    "shorter.faa": ">P1 short isoform\n" + "A" * 400 + "\n>P3 sole isoform\n" + "A" * 700 + "\n",
    "notreduced.faa": ">P1 short isoform\n" + "A" * 400 + "\n>P2 long isoform\n" + "A" * 900
                      + "\n>P3 sole isoform\n" + "A" * 700 + "\n",
    "lost.faa": ">P2 long isoform\n" + "A" * 900 + "\n",
    "stranger.faa": ">P2 long isoform\n" + "A" * 900 + "\n>P9 not in source\n" + "A" * 10 + "\n",
    "edited.faa": ">P2 long isoform\n" + "C" * 900 + "\n>P3 sole isoform\n" + "A" * 700 + "\n",
    "alien.faa": ">Z1 a\n" + "A" * 10 + "\n>Z2 b\n" + "A" * 20 + "\n",
    # header-encoded genes for the regex route: Q1/Q2 are one gene, Q3 has no recognisable gene
    "regex.faa": ">g1.m01.polypeptide a\n" + "A" * 100 + "\n>g1.m02.polypeptide b\n" + "A" * 300
                 + "\n>oddball c\n" + "A" * 50 + "\n",
}

# (label, argv template, expected exit code, substring the output must contain)
SELF_TEST = [
    ("reduce: longest isoform kept",
     ["{src}", "--out", "{tmp}/o.faa", "--gene-map-gff3", "{gff}"], 0, "1 record(s) dropped"),
    ("reduce: GFF3 has no protein_id",
     ["{src}", "--out", "{tmp}/o.faa", "--gene-map-gff3", "{tmp}/noprot.gff3"], 1,
     "no CDS feature carries"),
    ("reduce: Parent chain reaches no gene",
     ["{src}", "--out", "{tmp}/o.faa", "--gene-map-gff3", "{tmp}/orphan.gff3"], 1,
     "no gene feature anywhere up their"),
    ("reduce: protein absent from the GFF3",
     ["{tmp}/stranger.faa", "--out", "{tmp}/o.faa", "--gene-map-gff3", "{gff}"], 1,
     "not the same annotation release"),
    ("reduce: nothing to reduce is SAID so",
     ["{tmp}/good.faa", "--out", "{tmp}/o.faa", "--gene-map-gff3", "{gff}"], 0,
     "NOTHING WAS REDUCED"),
    ("cross-check: routes split a gene",
     ["{src}", "--out", "{tmp}/o.faa", "--gene-map-gff3", "{tmp}/split.gff3"], 1,
     "contradict each other"),
    ("cross-check: a shared symbol is allowed",
     ["{src}", "--out", "{tmp}/o.faa", "--gene-map-gff3", "{tmp}/merge.gff3"], 0,
     "share a symbol with another"),
    ("args: both grouping routes",
     ["{src}", "--out", "{tmp}/o.faa", "--gene-map-gff3", "{gff}",
      "--gene-regex", "(?P<gene>.)"], 2, "not allowed with"),
    ("args: no grouping route",
     ["{src}", "--out", "{tmp}/o.faa"], 2, "one of the arguments"),
    ("args: --out and --check together",
     ["{src}", "--out", "{tmp}/o.faa", "--check", "{tmp}/good.faa",
      "--gene-map-gff3", "{gff}"], 2, "not allowed with"),
    ("args: regex without a `gene` group",
     ["{src}", "--out", "{tmp}/o.faa", "--gene-regex", "(.*)"], 1, "no named group"),
    ("args: regex that will not compile",
     ["{src}", "--out", "{tmp}/o.faa", "--gene-regex", "(?P<gene>["], 1, "not a valid regex"),
    ("regex route: unmatched header is its own gene",
     ["{tmp}/regex.faa", "--out", "{tmp}/o.faa",
      "--gene-regex", r"^(?P<gene>[^ ]+?)\.m[0-9]+\.polypeptide"], 0,
     "1 of 3 header(s) did not match"),
    ("gzip decided by magic bytes, not suffix",
     ["{tmp}/gz_named_plain.faa", "--out", "{tmp}/o.faa", "--gene-map-gff3", "{gff}"], 0,
     "1 record(s) dropped"),
    ("check: a correct reduction passes",
     ["{src}", "--check", "{tmp}/good.faa", "--gene-map-gff3", "{gff}"], 0, "each the longest"),
    ("check: a shorter isoform is caught",
     ["{src}", "--check", "{tmp}/shorter.faa", "--gene-map-gff3", "{gff}"], 1,
     "kept a shorter isoform"),
    ("check: an unreduced file is caught",
     ["{src}", "--check", "{tmp}/notreduced.faa", "--gene-map-gff3", "{gff}"], 1,
     "do not contribute exactly one record"),
    ("check: a dropped gene is caught",
     ["{src}", "--check", "{tmp}/lost.faa", "--gene-map-gff3", "{gff}"], 1,
     "do not contribute exactly one record"),
    ("check: a foreign record is caught",
     ["{src}", "--check", "{tmp}/stranger.faa", "--gene-map-gff3", "{gff}"], 1,
     "are not in the source"),
    ("check: an edited sequence is caught",
     ["{src}", "--check", "{tmp}/edited.faa", "--gene-map-gff3", "{gff}"], 1,
     "must not edit sequences"),
    ("check: refuses the regex route as circular",
     ["{src}", "--check", "{tmp}/good.faa", "--gene-regex", "(?P<gene>P)"], 1,
     "--check needs --gene-map-gff3"),
    ("check: an empty gene map is NOT a pass",
     ["{tmp}/alien.faa", "--check", "{tmp}/alien.faa", "--gene-map-gff3", "{gff}"], 1,
     "nothing was actually checked"),
]


def self_test() -> int:
    """Run the real CLI against fixtures, because most guards live in argument handling.

    ⛔ EVERY CASE ASSERTS THE EXIT CODE **AND** A SUBSTRING OF THE MESSAGE. A guard that fires with
    the wrong explanation is the failure this project keeps hitting: the exit code is what CI reads,
    the message is what the next person acts on, and only checking the first lets the second rot.
    """
    import subprocess
    import tempfile

    bad = 0
    with tempfile.TemporaryDirectory() as td:
        tmp = pathlib.Path(td)
        gff = tmp / "ok.gff3"
        gff.write_text(_GFF3_OK)
        for name, text in _FASTA.items():
            (tmp / name).write_text(text)
        # a gzip stream deliberately NOT named .gz, to prove the magic-byte reader
        with gzip.open(tmp / "gz_named_plain.faa", "wt") as fh:
            fh.write(_FASTA["src.faa"])
        (tmp / "noprot.gff3").write_text(
            "".join(ln + "\n" for ln in _GFF3_OK.splitlines() if "protein_id=" not in ln))
        (tmp / "orphan.gff3").write_text(
            _GFF3_OK.replace(";Parent=gene-G1", "").replace(";Parent=gene-G2", ""))
        # one gene feature whose CDS rows claim two different symbols
        (tmp / "split.gff3").write_text(_GFF3_OK.replace("gene=SYM_1;protein_id=P2",
                                                         "gene=SYM_OTHER;protein_id=P2"))
        # two gene features sharing one symbol -- legitimate, must not refuse
        (tmp / "merge.gff3").write_text(_GFF3_OK.replace("gene=SYM_2", "gene=SYM_1"))

        for label, argv, want_rc, want_text in SELF_TEST:
            args = [x.format(tmp=tmp, gff=gff, src=tmp / "src.faa") for x in argv]
            p = subprocess.run([sys.executable, __file__, *args],
                               capture_output=True, text=True)
            got = p.stdout + p.stderr
            ok = p.returncode == want_rc and want_text in got
            bad += not ok
            note = ""
            if not ok:
                note = (f" -- wanted rc={want_rc} and {want_text!r}, got rc={p.returncode}: "
                        f"{got.strip().splitlines()[-1] if got.strip() else '<no output>'}")
            print(f"  {'pass' if ok else '⛔ FAIL'}  {label:46} rc={p.returncode}{note}")

        # the happy path must keep the LONG isoform, not merely one of them
        p = subprocess.run([sys.executable, __file__, str(tmp / "src.faa"), "--out",
                            str(tmp / "hp.faa"), "--gene-map-gff3", str(gff)],
                           capture_output=True, text=True)
        kept = sorted(h.split()[0] for h, _ in records(tmp / "hp.faa"))
        ok = p.returncode == 0 and kept == ["P2", "P3"]
        bad += not ok
        print(f"  {'pass' if ok else '⛔ FAIL'}  {'the kept record is the longer one':46} "
              f"kept={kept}")

    n = len(SELF_TEST) + 1
    print(f"\n  {n - bad}/{n} self-test case(s) behaved correctly")
    return 1 if bad else 0


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", type=pathlib.Path, nargs="?",
                    help="the FULL proteome; omit only with --self-test")
    ap.add_argument("--self-test", action="store_true",
                    help="prove each guard fires, with the right message")
    # ⚠ INTERCEPTED BEFORE parse_args, deliberately. Both option groups below are
    # `required=True`, so a bare --self-test would be refused by argparse before reaching any
    # branch -- and those groups being required is itself something the self-test asserts, so it
    # cannot be relaxed to make room for this flag.
    if "--self-test" in argv:
        return self_test()
    what = ap.add_mutually_exclusive_group(required=True)
    what.add_argument("--out", type=pathlib.Path,
                      help="written gzipped if the name ends .gz")
    what.add_argument("--check", type=pathlib.Path, metavar="REDUCED",
                      help="write nothing; instead assert that this already-reduced proteome IS "
                           "one-longest-protein-per-gene against the source FASTA. The audit for a "
                           "proteome someone else reduced, or one this wrote long ago.")
    how = ap.add_mutually_exclusive_group(required=True)
    how.add_argument("--gene-regex",
                     help="regex with a named group `gene`, matched against the FASTA header. "
                          "⚠ a pattern that matches nothing leaves the file unchanged and looks "
                          "like success, which is why one of these two is REQUIRED.")
    how.add_argument("--gene-map-gff3", type=pathlib.Path,
                     help="take the gene from this GFF3's CDS->mRNA->gene chain, keyed by the "
                          "FASTA header's first token. For proteomes whose headers carry no gene "
                          "at all, such as NCBI's `*_protein.faa.gz`.")
    a = ap.parse_args(argv)
    if a.fasta is None:
        ap.error("the input proteome is required")

    pat = gmap = None
    if a.gene_regex is not None:
        try:
            pat = re.compile(a.gene_regex)
        except re.error as e:
            sys.exit(f"--gene-regex is not a valid regex: {e}")
        if "gene" not in (pat.groupindex or {}):
            sys.exit(f"--gene-regex {a.gene_regex!r} has no named group `gene`. Write it as "
                     f"`(?P<gene>...)` so this knows which part identifies the gene.")
    else:
        gmap = gene_map_from_gff3(a.gene_map_gff3)

    if a.check is not None:
        if gmap is None:
            sys.exit("--check needs --gene-map-gff3: checking one-per-gene against a --gene-regex "
                     "would re-derive the grouping from the same rule that produced the file, "
                     "which cannot fail. Check against the annotation instead.")
        return check_reduced(a.fasta, a.check, gmap)

    best: dict[str, tuple[int, str, str]] = {}
    total = unmatched = 0
    missing: list[str] = []
    for header, seq in records(a.fasta):
        total += 1
        if gmap is not None:
            acc = header.split()[0]
            gene = gmap.get(acc)
            if gene is None:
                # ⛔ REFUSE, don't fall back. See the module docstring: a protein the annotation
                # does not mention means the FASTA and the GFF3 are different releases.
                missing.append(acc)
                continue
        else:
            m = pat.search(header)
            if not m:
                unmatched += 1
                # ⛔ AN UNMATCHED HEADER IS ITS OWN GENE, NOT A DISCARD. Dropping it would silently
                # shorten the proteome, and a proteome that is quietly short is the failure this
                # whole exercise is about. It is counted and reported instead.
                gene = header.split()[0]
            else:
                gene = m.group("gene")
        prev = best.get(gene)
        if prev is None or len(seq) > prev[0]:
            best[gene] = (len(seq), header, seq)

    if missing:
        sys.exit(f"{a.fasta}: {len(missing)} of {total} protein(s) are absent from "
                 f"{a.gene_map_gff3}, e.g. {', '.join(missing[:5])}. The FASTA and the GFF3 are "
                 f"not the same annotation release -- reducing the rest would produce a proteome "
                 f"that looks complete and is silently mixed.")
    if unmatched:
        print(f"  ⚠ {unmatched} of {total} header(s) did not match --gene-regex and were each "
              f"treated as their own gene. If that is most of the file, the pattern is wrong.")

    groups = total - len(best)
    print(f"  {total} record(s) -> {len(best)} gene(s); {groups} record(s) dropped as shorter "
          f"isoforms")
    if groups == 0:
        # Not an error: plenty of published proteomes are natively one-per-gene. But it is
        # indistinguishable from a pattern that matched nothing, so say which it was.
        detail = (f"the gene map grouped nothing ({len(set(gmap.values()))} gene(s) for "
                  f"{len(gmap)} protein(s))" if gmap is not None else
                  f"--gene-regex grouped nothing ({total - unmatched} header(s) did match)")
        print(f"  ⚠ NOTHING WAS REDUCED. Either this proteome is already one protein per gene, or "
              f"{detail}. Check before assuming the former.")

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
