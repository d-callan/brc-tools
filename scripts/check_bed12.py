#!/usr/bin/env python3
"""Check that a BED12 is structurally what its consumers assume.

    python3 scripts/check_bed12.py                 # every committed 12-column BED
    python3 scripts/check_bed12.py a.bed12 b.bed   # named files
    python3 scripts/check_bed12.py --self-test     # prove each rule bites

⛔ NOTHING IN THIS REPOSITORY CHECKED BED12 STRUCTURE, AND TWO PIPELINES KEY ON IT. WF-A publishes
`anchor_bed12s` from `gffread --bed` (cs10: 40,185 rows), TOGA2 takes it as its reference gene set,
and `phase_c4_merge` keys its reference genes off COLUMN 4. A malformed block list does not error --
it projects the wrong exon boundaries and reports success, which is the failure class this
repository keeps finding.

⚠ EACH RULE CARRIES ITS SOURCE, because "BED12" is a spec plus a convention and the two are not the
same thing. Two of these rules are NOT in the UCSC field descriptions; they are what UCSC's own
`validateFiles` enforces and what every consumer assumes, so they are checked and labelled as such
rather than presented as spec.

THE COORDINATE CONVENTION, from bedtools' general-usage page (which is unusually explicit):
`start` is "the zero-based starting position", `end` is "the one-based ending position", so
`start=9, end=20` spans bases 10 through 20 inclusive. `blockStarts` are offsets relative to the
feature's start. `thickStart`, `thickEnd`, `itemRgb`, `blockCount` and `blockSizes` are "allowed yet
ignored by bedtools" -- ignored by bedtools does not mean ignored by TOGA2, which is why they are
checked here.
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]

#: rule -> where it comes from. Printed with --rules so the distinction is not lost.
SOURCES = {
    "12 fields": "UCSC BED12 / bedtools: chrom start end name score strand thickStart thickEnd "
                 "itemRgb blockCount blockSizes blockStarts",
    "0 <= start < end": "bedtools: start is zero-based, end is the one-based end",
    "strand is + or -": "bedtools: \"either '+' or '-'\". A '.' is reported, not failed: UCSC "
                        "permits it and a gene model simply should not use it",
    "blockCount matches both lists": "UCSC: \"The number of items in this list should correspond "
                                     "to blockCount\"",
    "block sizes are positive": "a zero-length exon is not a feature",
    "blocks ascending and non-overlapping": "implied by blockStarts being offsets into one feature",
    "blockStarts[0] == 0": "⚠ CONVENTION, not the field description -- enforced by UCSC "
                           "validateFiles and assumed by every reader",
    "last block ends at end-start": "⚠ CONVENTION, as above: the blocks must span the feature",
    "start <= thickStart <= thickEnd <= end": "the thick range is a sub-range of the feature",
}


def check_text(text: str, label: str) -> tuple[int, list[str], dict[str, int]]:
    """Return (rows checked, failures, notes) for one BED12 body."""
    fails: list[str] = []
    notes: dict[str, int] = {}
    rows = 0
    for n, line in enumerate(text.splitlines(), 1):
        if not line.strip():
            continue
        if line.startswith(("#", "track", "browser")):
            # ⚠ SKIPPED BUT COUNTED. UCSC permits these; silently ignoring them would let a file
            # that is entirely track lines report "0 rows, all good".
            notes["track/browser/comment lines"] = notes.get("track/browser/comment lines", 0) + 1
            continue
        rows += 1
        f = line.split("\t")
        if len(f) != 12:
            fails.append(f"{label}:{n}: {len(f)} fields, need exactly 12")
            continue
        # name/score/itemRgb are unpacked for documentation: this is the 12-field order, and
        # naming them is how a reader checks the positional unpacking against the spec.
        chrom, s_, e_, _name, _score, strand, ts_, te_, _rgb, bc_, bs_, bt_ = f
        try:
            s, e, ts, te, bc = int(s_), int(e_), int(ts_), int(te_), int(bc_)
            sizes = [int(x) for x in bs_.rstrip(",").split(",") if x != ""]
            starts = [int(x) for x in bt_.rstrip(",").split(",") if x != ""]
        except ValueError:
            fails.append(f"{label}:{n}: a numeric field is not an integer: "
                         f"start={s_!r} end={e_!r} thick={ts_!r},{te_!r} count={bc_!r} "
                         f"sizes={bs_!r} starts={bt_!r}")
            continue
        if not chrom:
            fails.append(f"{label}:{n}: empty chrom")
        if not (0 <= s < e):
            fails.append(f"{label}:{n}: start={s} end={e}; need 0 <= start < end")
        if strand == ".":
            notes["strand '.'"] = notes.get("strand '.'", 0) + 1
        elif strand not in ("+", "-"):
            fails.append(f"{label}:{n}: strand {strand!r}, need '+' or '-'")
        if bc != len(sizes) or bc != len(starts):
            fails.append(f"{label}:{n}: blockCount={bc} but {len(sizes)} size(s) and "
                         f"{len(starts)} start(s)")
            continue
        if not sizes:
            fails.append(f"{label}:{n}: no blocks")
            continue
        if any(x <= 0 for x in sizes):
            fails.append(f"{label}:{n}: a block size is not positive: {sizes[:6]}")
        if starts[0] != 0:
            fails.append(f"{label}:{n}: blockStarts[0]={starts[0]}, must be 0 -- the offsets are "
                         f"relative to the feature start, so the first block starts there")
        if starts[-1] + sizes[-1] != e - s:
            fails.append(f"{label}:{n}: last block ends at {starts[-1] + sizes[-1]} but the "
                         f"feature is {e - s} long; the blocks must span it")
        for i in range(len(starts) - 1):
            if starts[i] + sizes[i] > starts[i + 1]:
                fails.append(f"{label}:{n}: block {i} ends at {starts[i] + sizes[i]} and block "
                             f"{i + 1} starts at {starts[i + 1]}; blocks overlap or are unsorted")
                break
        if not (s <= ts <= te <= e):
            fails.append(f"{label}:{n}: thick range {ts}-{te} is not inside {s}-{e}")
    return rows, fails, notes


def committed_bed12() -> list[pathlib.Path]:
    """Every git-tracked file whose first data line has 12 tab-separated fields.

    ⚠ FOUND BY SHAPE, NOT BY SUFFIX. Two of the four in this repo are named `.bed`, not `.bed12`
    (`phase_c2_triage/test-data/reference.bed` and `phase_c4_merge/test-data/
    toga2_query_annotation.bed`), so a suffix glob would check half of them and report a pass.
    """
    out = subprocess.run(["git", "-C", str(ROOT), "ls-files"],
                         capture_output=True, text=True, check=True).stdout.split()
    found = []
    for rel in out:
        if not rel.lower().endswith((".bed", ".bed12")):
            continue
        p = ROOT / rel
        try:
            text = p.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for line in text.splitlines():
            if line.strip() and not line.startswith(("#", "track", "browser")):
                if len(line.split("\t")) == 12:
                    found.append(p)
                break
    return found


#: (label, one row, does it fail?) -- the row is otherwise valid so each case isolates one rule.
GOOD = "chr1\t100\t200\tg1\t0\t+\t100\t200\t0\t2\t40,60,\t0,40,"
SELF_TEST = [
    ("a valid row", GOOD, False),
    ("11 fields", GOOD.rsplit("\t", 1)[0], True),
    ("13 fields", GOOD + "\textra", True),
    ("start == end", "chr1\t100\t100\tg1\t0\t+\t100\t100\t0\t1\t0,\t0,", True),
    ("start > end", "chr1\t200\t100\tg1\t0\t+\t200\t200\t0\t1\t1,\t0,", True),
    ("negative start", "chr1\t-1\t200\tg1\t0\t+\t0\t200\t0\t1\t201,\t0,", True),
    ("bad strand", GOOD.replace("\t+\t", "\t?\t"), True),
    ("strand '.' is a note, not a failure", GOOD.replace("\t+\t", "\t.\t"), False),
    ("non-integer start", GOOD.replace("\t100\t200\t", "\tx\t200\t", 1), True),
    ("blockCount disagrees with the lists", GOOD.replace("\t2\t", "\t3\t"), True),
    ("blockStarts[0] != 0", "chr1\t100\t200\tg1\t0\t+\t100\t200\t0\t2\t40,60,\t1,40,", True),
    ("last block does not reach the end", "chr1\t100\t200\tg1\t0\t+\t100\t200\t0\t2\t40,59,\t0,40,",
     True),
    ("blocks overlap", "chr1\t100\t200\tg1\t0\t+\t100\t200\t0\t2\t50,60,\t0,40,", True),
    ("zero-length block", "chr1\t100\t200\tg1\t0\t+\t100\t200\t0\t2\t0,100,\t0,100,", True),
    ("thickStart before start", "chr1\t100\t200\tg1\t0\t+\t99\t200\t0\t1\t100,\t0,", True),
    ("thickEnd after end", "chr1\t100\t200\tg1\t0\t+\t100\t201\t0\t1\t100,\t0,", True),
    ("thickEnd before thickStart", "chr1\t100\t200\tg1\t0\t+\t150\t120\t0\t1\t100,\t0,", True),
    ("empty chrom", "\t100\t200\tg1\t0\t+\t100\t200\t0\t1\t100,\t0,", True),
]


def self_test() -> int:
    """Each rule must reject the row that violates it and accept the one that does not."""
    bad = 0
    for label, row, should_fail in SELF_TEST:
        _, fails, _ = check_text(row, "case")
        ok = bool(fails) == should_fail
        bad += not ok
        want = "rejection" if should_fail else "acceptance"
        note = "" if ok else f" -- expected {want}"
        print(f"  {'pass' if ok else '⛔ FAIL'}  {label:44} "
              f"{'rejected' if fails else 'accepted'}{note}")
    # ⛔ A FILE OF NOTHING BUT TRACK LINES MUST NOT PASS AS "all good". 0 rows is not a clean file.
    rows, fails, notes = check_text("track name=x\n# a comment\n", "case")
    ok = rows == 0 and notes.get("track/browser/comment lines") == 2
    bad += not ok
    print(f"  {'pass' if ok else '⛔ FAIL'}  {'track/comment lines counted, 0 data rows':44} "
          f"rows={rows} notes={notes}")
    n = len(SELF_TEST) + 1
    print(f"\n  {n - bad}/{n} self-test case(s) behaved correctly")
    return 1 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="*", type=pathlib.Path)
    ap.add_argument("--self-test", action="store_true", help="prove each rule rejects its violation")
    ap.add_argument("--rules", action="store_true", help="print each rule and where it comes from")
    a = ap.parse_args()
    if a.rules:
        for rule, src in SOURCES.items():
            print(f"  {rule}\n      {src}")
        return 0
    if a.self_test:
        return self_test()

    paths = a.files or committed_bed12()
    if not paths:
        # ⛔ NOTHING TO CHECK IS NOT A PASS, for the same reason a skip is not a pass anywhere else
        # here: it is indistinguishable from a clean run and CI reads only the exit code.
        sys.exit("no BED12 files found. Pass paths explicitly, or check why the repository scan "
                 "returned nothing -- a silent zero here reads exactly like a clean check.")
    total_rows, total_fails = 0, 0
    for p in paths:
        rows, fails, notes = check_text(p.read_text(encoding="utf-8", errors="replace"),
                                        p.relative_to(ROOT) if p.is_relative_to(ROOT) else p)
        total_rows += rows
        total_fails += len(fails)
        extra = "".join(f", {v} {k}" for k, v in sorted(notes.items()))
        shown = p.relative_to(ROOT) if p.is_relative_to(ROOT) else p
        print(f"  {'⛔' if fails else 'ok'}  {shown!s:58} {rows:>7,} row(s){extra}")
        for f in fails[:5]:
            print(f"        {f}")
        if len(fails) > 5:
            print(f"        ... and {len(fails) - 5} more")
        if rows == 0 and not fails:
            total_fails += 1
            print("        ⛔ 0 data rows -- an empty BED12 is not a valid gene set")
    print(f"\n  {len(paths)} file(s), {total_rows:,} row(s), {total_fails} failure(s)")
    return 1 if total_fails else 0


if __name__ == "__main__":
    sys.exit(main())
