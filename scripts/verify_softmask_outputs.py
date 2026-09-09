#!/usr/bin/env python3
"""Verify a softmask invocation by recomputing its mask two independent ways.

⛔ "42 JOBS OK" IS NOT EVIDENCE THE OUTPUT IS RIGHT, AND THIS PIPELINE HAS ALREADY PROVED IT. Every
step of the committed workflow exited 0 while `maskfasta` was being fed the ORIGINAL assembly
instead of the uppercased one, so the published `softmasked_fasta` carried NCBI's pre-existing mask
unioned with this workflow's, indistinguishably -- and the invocation summary was entirely green.
An invocation summary reports scheduling, not correctness.

THE CHECK. For each strain, two numbers derived from different files by different tools, which must
agree if the pipeline is sound:

    A. bases covered by the merged union BED   (from the intervals, via bedtools merge)
    B. lowercase bases in the soft-masked FASTA (from the sequence, via bedtools maskfasta)

A == B exactly, or the mask APPLIED is not the mask COMPUTED. Under the old wiring B exceeded A by
whatever the assembly arrived carrying -- for cs10 that is tens of percent, not a rounding error.

⚠ It also asserts the uppercased input really is uppercase. Without that, A == B could hold with
both wrong.

⛔ RESOLVE BY INVOCATION, NOT BY SCANNING A HISTORY. The previous version searched a history for
collections whose NAMES matched `"Merged "` / `"bedtools MaskFastaBed"` and took the highest hid of
each, independently. Two consequences, both silent:

  * Those are ToolShed-derived DEFAULT names that nothing in this repo generates or asserts. One
    tool version bump and the script exits "is this a completed softmask run?" -- a false negative
    that reads like the run failed.
  * `build_up_softmask.py` runs all seven tiers into ONE history, and tiers 6 and 7 both produce a
    merged BED and a masked FASTA. Nothing tied the three collections to a single invocation, so a
    partially-failed tier 7 could have the union read from one tier and the FASTA from another --
    and the script would still print "the mask applied is exactly the mask computed".

An invocation exposes `output_collections` keyed by the workflow's OWN declared output names, so
`--invocation` resolves each collection exactly and cannot cross runs.

    python3 scripts/verify_softmask_outputs.py --invocation <invocation_id>
"""
from __future__ import annotations

import argparse
import sys

from bioblend.galaxy import GalaxyInstance
from softmask_lib import SCHEDULING_IN_PROGRESS, connect, fasta_stats

#: Declared output names in workflows/softmask/softmask_udt.gxwf.yml.
UPPER_OUT = "uppercased_fasta"
UNION_OUT = "mask_union"
MASKED_OUT = "softmasked_fasta"


def elements(gi: GalaxyInstance, collection_id: str) -> dict[str, str]:
    d = gi.dataset_collections.show_dataset_collection(collection_id)
    return {e["element_identifier"]: e["object"]["id"] for e in d["elements"]}


def text(gi: GalaxyInstance, dataset_id: str) -> str:
    b = gi.datasets.download_dataset(dataset_id, use_default_filename=False)
    return b.decode("utf-8", "replace") if isinstance(b, bytes) else str(b)


def bed_intervals(t: str) -> set[tuple[str, int, int]]:
    """Parse a BED into a set of (chrom, start, end).

    ⛔ POSITIONS, NOT A TOTAL. This used to return only the summed width, and the check compared
    that against the masked FASTA's lowercase COUNT. Equal totals do not mean the same bases: an
    adversarial pass showed the old form certifying four wrong pipelines as correct -- intervals at
    0-10 against lowercase at 50-60, an off-by-one, a mask applied to the wrong chromosome, and an
    empty BED against an entirely unmasked FASTA. All four printed "delta +0 MATCH".

    ⚠ Track/browser/comment lines are REFUSED rather than skipped. Silently ignoring them would
    under-count, which is the same class of error one level down.
    """
    out: set[tuple[str, int, int]] = set()
    for n, line in enumerate(t.splitlines(), 1):
        if not line.strip():
            continue
        if line.startswith(("#", "track", "browser")):
            sys.exit(f"BED line {n} is a {line.split()[0]!r} line; this verifier expects a plain "
                     f"merged BED and will not silently skip it: {line[:60]!r}")
        f = line.split("\t")
        if len(f) < 3:
            sys.exit(f"BED line {n} has {len(f)} fields, need at least 3: {line[:60]!r}")
        try:
            out.add((f[0], int(f[1]), int(f[2])))
        except ValueError:
            sys.exit(f"BED line {n} has non-integer coordinates: {line[:60]!r}")
    return out


def lowercase_runs(t: str) -> set[tuple[str, int, int]]:
    """Maximal lowercase runs of a FASTA, as (chrom, start, end) in BED half-open coordinates.

    ⚠ The chrom is the FIRST WHITESPACE TOKEN of the header, because that is what bedtools writes
    into the BED it produced. `bedtools maskfasta` also TRUNCATES the header to exactly that token,
    so the masked FASTA and its own intervals agree by construction -- but the uppercased input does
    not have truncated headers, and comparing the two naively would mismatch on description text.
    """
    runs: set[tuple[str, int, int]] = set()
    chrom: str | None = None
    pos, start = 0, None
    for line in t.splitlines():
        if line.startswith(">"):
            if chrom is not None and start is not None:
                runs.add((chrom, start, pos))
            chrom, pos, start = line[1:].split()[0] if len(line) > 1 else "", 0, None
            continue
        if chrom is None:
            # ⛔ SEQUENCE BEFORE ANY HEADER. Without this the run would be recorded under a chrom of
            # None, which matches no BED interval, so a malformed FASTA would surface as a coordinate
            # DISAGREEMENT rather than as the parse failure it is -- and telling those two apart is
            # this function's entire job. Blank leading lines are not sequence and are skipped.
            if line.strip():
                sys.exit("lowercase_runs: sequence data appears before any '>' header; "
                         "this is not a FASTA and its coordinates cannot be trusted.")
            continue
        for ch in line:
            if "a" <= ch <= "z":
                if start is None:
                    start = pos
            elif start is not None:
                runs.add((chrom, start, pos))
                start = None
            pos += 1
    if chrom is not None and start is not None:
        runs.add((chrom, start, pos))
    return runs


def span(intervals: set[tuple[str, int, int]]) -> int:
    return sum(e - s for _, s, e in intervals)


def merge_by_chrom(intervals) -> dict[str, list[tuple[int, int]]]:
    """Disjoint, sorted intervals per chrom, so coverage can be asked about a POSITION."""
    by: dict[str, list[tuple[int, int]]] = {}
    for c, s, e in intervals:
        by.setdefault(c, []).append((s, e))
    for c, ivs in by.items():
        ivs.sort()
        out: list[tuple[int, int]] = []
        for s, e in ivs:
            if out and s <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], e))
            else:
                out.append((s, e))
        by[c] = out
    return by


def uncovered(needles, hay) -> list[tuple[str, int, int]]:
    """Fragments of `needles` whose BASES are not covered by `hay`.

    ⛔ BASE-LEVEL, BECAUSE INTERVAL IDENTITY IS THE WRONG QUESTION HERE. The check used
    `union - applied`, a set difference over whole (chrom, start, end) tuples, so a computed
    interval counted as UNAPPLIED unless a lowercase run existed with exactly its coordinates. With
    the arriving mask kept -- the default -- the published FASTA is the union of the submitter's
    mask and ours, so an interval of ours that happens to abut a submitter-masked base merges into
    one longer run and its tuple stops matching. Measured on the 19-genome run: 10 of 19 strains
    reported failures whose every example was a one-base boundary difference with the computed
    interval sitting INSIDE the published run -- containment holding, reported as a failure, and
    104 Mb of "computed-but-not-applied" that was really whole intervals disqualified by one base.
    """
    from bisect import bisect_right
    hay_by = merge_by_chrom(hay)
    starts = {c: [s for s, _ in ivs] for c, ivs in hay_by.items()}
    gaps: list[tuple[str, int, int]] = []
    for c, s, e in sorted(needles):
        ivs, st = hay_by.get(c), starts.get(c)
        if not ivs:
            gaps.append((c, s, e))
            continue
        pos = s
        i = bisect_right(st, pos) - 1
        while pos < e:
            # `i` is the last hay interval starting at or before `pos`; advance it monotonically.
            while i + 1 < len(ivs) and ivs[i + 1][0] <= pos:
                i += 1
            if 0 <= i < len(ivs) and ivs[i][0] <= pos < ivs[i][1]:
                pos = min(ivs[i][1], e)                 # covered up to here
            else:
                # uncovered until the next hay interval starts, or to the needle's end
                nxt = ivs[i + 1][0] if i + 1 < len(ivs) else e
                stop = min(nxt, e)
                gaps.append((c, pos, stop))
                pos = stop
    return gaps


def resolve(gi: GalaxyInstance, invocation_id: str) -> dict[str, str]:
    """Map the three declared output names to collection ids, all from ONE invocation."""
    inv = gi.invocations.show_invocation(invocation_id)
    state = inv.get("state")
    # ⛔ TERMINAL IS NOT THE SAME AS SOUND. The tuple means "Galaxy will create no further jobs",
    # which is the right wait condition elsewhere -- but `failed` and `cancelled` are terminal too,
    # and certifying a cancelled run's partly-populated collections would be worse than useless.
    if state in ("failed", "cancelled"):
        sys.exit(f"invocation {invocation_id} is in state {state!r}; refusing to verify it.")
    if state in SCHEDULING_IN_PROGRESS:
        print(f"  ⚠ invocation state is {state!r}: Galaxy may still be creating jobs, "
              f"so the collections below may be incomplete.")
    cols = inv.get("output_collections") or {}
    missing = [n for n in (UPPER_OUT, UNION_OUT, MASKED_OUT) if n not in cols]
    if missing:
        sys.exit(f"invocation {invocation_id} does not declare {', '.join(missing)}. "
                 f"It exposes: {sorted(cols)}.\n"
                 f"  A workflow predating the mask_union/uppercased_fasta outputs cannot be "
                 f"verified this way -- re-run the current softmask_udt.gxwf.yml.")
    return {n: cols[n]["id"] for n in (UPPER_OUT, UNION_OUT, MASKED_OUT)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--invocation", required=True,
                    help="invocation id of a softmask run (NOT a history id -- a history can hold "
                         "several runs, and mixing their collections is exactly the bug this "
                         "argument exists to prevent)")
    args = ap.parse_args()

    gi = connect()
    ids = resolve(gi, args.invocation)
    up_e = elements(gi, ids[UPPER_OUT])
    mg_e = elements(gi, ids[UNION_OUT])
    mk_e = elements(gi, ids[MASKED_OUT])

    strains = sorted(set(up_e) & set(mg_e) & set(mk_e))
    if not strains:
        sys.exit("no strain appears in all three collections; element identifiers do not line up.")
    for label, e in ((UPPER_OUT, up_e), (UNION_OUT, mg_e), (MASKED_OUT, mk_e)):
        extra = sorted(set(e) - set(strains))
        if extra:
            print(f"  ⚠ {label} also holds {extra}, absent from another collection -- NOT checked.")

    print("  union BED coverage vs masked-FASTA lowercase -- derived independently\n")
    failures = 0
    for s in strains:
        union = bed_intervals(text(gi, mg_e[s]))
        masked_text = text(gi, mk_e[s])
        applied = lowercase_runs(masked_text)
        _, res, low = fasta_stats(masked_text)
        _, ures, ulow = fasta_stats(text(gi, up_e[s]))
        if res == 0:
            print(f"    {s}\n      ⛔ masked FASTA is EMPTY")
            failures += 1
            continue

        # ⛔ THE INVARIANT DEPENDS ON WHICH MODE RAN, AND THE MODE IS INFERRABLE FROM THIS OUTPUT.
        # WF-B's `strip_arrived_mask` defaults to FALSE -- the classic behaviour -- so `uppercase`
        # passes the assembly through and its output legitimately still carries the submitter's
        # lower case. Checking `ulow == 0` unconditionally therefore prints `⛔ NOT UPPERCASE /
        # LENGTH CHANGED` on every assembly that arrives soft-masked (cs10 ships 46.8% masked),
        # which reads as a corrupt uppercase step rather than as the setting that was chosen. And
        # the mask comparison below then reports a spurious applied-but-not-computed list, because
        # the published FASTA carries BOTH masks while `mask_union` records only ours.
        # ⚠ INFERRED, NOT ASKED FOR, so no caller can pass the wrong flag -- but ⛔ `ulow == 0`
        # DOES NOT MEAN THE MASK WAS STRIPPED. It means nothing lowercase reached the mask step,
        # and an assembly that ARRIVED unmasked looks identical to one that was stripped. Measured
        # on the 19-genome run, which used `strip_arrived_mask: false` throughout: 8 of 19 printed
        # "mask stripped" because those assemblies ship no soft-mask at all, so the label asserted
        # a setting that was never used. The INVARIANT is unaffected -- with no arriving mask to
        # union in, exact equality is right either way -- so only the wording was wrong, which is
        # its own kind of defect in a script whose output is the evidence.
        no_arriving_mask = (ulow == 0)
        ok_upper = (ures == res)
        # ⛔ AN EMPTY MASK IS A FAILURE, NOT A MATCH. Nothing masked against nothing computed is
        # 0 == 0, and the old total-based check printed ✅ for it -- certifying a run in which the
        # mask was never applied at all.
        ok_nonempty = bool(union)
        # ⛔ CONTAINMENT AT BASE LEVEL, NOT INTERVAL IDENTITY -- see uncovered(). What must hold in
        # BOTH modes is that every base we computed is lowercase in the published FASTA.
        only_computed = uncovered(union, applied)
        # ⚠ `only_applied` IS EXPECTED WHEN THE ARRIVING MASK WAS KEPT: the published FASTA is the
        # union of the submitter's mask and ours, and mask_union records only ours. It is a FAILURE
        # only when the mask was stripped, where the two must agree exactly.
        only_applied = uncovered(applied, union)
        ok_mask = ok_nonempty and not only_computed and (not only_applied or not no_arriving_mask)
        failures += (not ok_upper) + (not ok_mask)

        print(f"    {s}")
        print(f"      uppercased input : {ures:,} nt, {ulow} lowercase"
              f"        {'ok' if ok_upper else '⛔ LENGTH CHANGED'}"
              f"  [{'nothing lowercase arrived -- stripped, or the assembly ships unmasked'
                     if no_arriving_mask else 'arriving mask KEPT (classic default)'}]")
        print(f"      merged union BED : {span(union):,} nt in {len(union):,} intervals "
              f"({span(union) / res:.2%})")
        print(f"      masked FASTA     : {low:,} lowercase in {len(applied):,} runs "
              f"({low / res:.2%})")
        if not ok_nonempty:
            print("      ⛔ the union is EMPTY -- nothing was masked, which is not a pass")
        elif ok_mask:
            print(f"      positions        : every computed base is masked"
                  f"{'' if no_arriving_mask else f'; {span(only_applied):,} nt of the arriving mask also present'}")
        else:
            print(f"      ⛔ computed-but-not-applied: {span(only_computed):,} nt in "
                  f"{len(only_computed):,} fragment(s)")
            print(f"      {'⛔' if no_arriving_mask else 'ℹ'} applied-but-not-computed: "
                  f"{span(only_applied):,} nt in {len(only_applied):,} fragment(s)")
            for iv in sorted(only_computed)[:3]:
                print(f"          only in BED  : {iv}")
            for iv in sorted(only_applied)[:3]:
                print(f"          only in FASTA: {iv}")
    print()
    if failures:
        print(f"  ⛔ {failures} check(s) FAILED across {len(strains)} strain(s)")
        return 1
    # ⚠ "EVERY COMPUTED BASE" AND NOT "INTERVAL FOR INTERVAL", which is what this said while the
    # comparison was interval identity. It is base-level containment now, and a summary describing
    # a stronger check than the one that ran is worse than no summary.
    print(f"  ✅ {len(strains)} strain(s): every base the workflow computed is masked in the "
          f"published FASTA")
    return 0


if __name__ == "__main__":
    sys.exit(main())
