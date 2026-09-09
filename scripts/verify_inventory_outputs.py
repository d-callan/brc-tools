#!/usr/bin/env python3
"""Verify a WF-A (inventory) invocation against the PANEL it was run on.

⛔ AN INVOCATION SUMMARY REPORTS SCHEDULING, NOT CORRECTNESS -- and WF-A is the workflow where that
gap is widest, because six of its thirteen outputs are files whose CONTENT is a count. A relabel map
with the wrong number of rows, a sourmash matrix over fewer genomes than were staged, a self-pairs
list missing a strain: every one of those is a green run and a wrong answer, and every one is
consumed by a later workflow that will not notice either.

⚠ THE PANEL IS THE SECOND ROUTE, WHICH IS WHY IT IS REQUIRED. Checking the outputs against each
other only proves they agree; checking them against the panel definition proves they describe the
genomes someone meant to run. Both are done here: `n` is taken from the invocation's OWN assemblies
input and then required to equal the panel's asserted count, so a panel that quietly resolved short
fails even though every output is internally consistent.

⛔ AND THE THREE COLLECTIONS HAVE THREE DIFFERENT LENGTHS ON PURPOSE. `proteomes` is a subset (only
members with an annotation) and `anchor_gene_gff3s` a subset of that, so "all collections have n
elements" is the WRONG assertion -- it would have to be relaxed on the first real panel, and a
relaxed assertion is one nobody trusts. Each output is checked against the collection it is keyed
by, named in KEYED_BY below.

    python3 scripts/verify_inventory_outputs.py --invocation <id> --panel panels/x.panel.yml
    python3 scripts/verify_inventory_outputs.py --self-test    # prove each check bites

Controls from an earlier run are passed in, not baked in:

    --expect anchor_bed12s:cs10_NCBI_RefSeq_softmasked=40185

⚠ A CONTROL IS THE ONLY THING HERE THAT CATCHES AN UNCHANGED BRANCH GOING WRONG. Growing the panel
must not touch the anchor branch at all, so an anchor row count that moved is a regression none of
the count checks can see -- they would happily accept a different number of rows.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
from typing import TYPE_CHECKING

import yaml

# ⚠ bioblend AND softmask_lib ARE IMPORTED LAZILY, INSIDE THE FUNCTIONS THAT TALK TO A SERVER.
# CI installs pyyaml and galaxy-tool-util, not bioblend, and the whole point of --self-test is that
# it needs no server -- a module-level import would make the one part CI *can* run the one part it
# cannot.
if TYPE_CHECKING:                                                     # pragma: no cover
    from bioblend.galaxy import GalaxyInstance

#: Declared output names in workflows/inventory/inventory_udt.gxwf.yml, each with the collection
#: whose identifiers it must be keyed by -- `None` for a single dataset.
KEYED_BY = {
    "signatures": "assemblies",
    "sizes": "assemblies",
    "fasta_index": "assemblies",
    "busco_summaries": "proteomes",
    "anchor_bed12s": "anchors",
    "anchor_isoforms": "anchors",
    "similarity_matrix": None,
    "sourmash_heatmap": None,
    "sourmash_dendrogram": None,
    "qc_report": None,
    "panel_identifiers": None,
    "self_pairs": None,
    "relabel_map": None,
}


def rows(s: str) -> list[str]:
    """Non-blank lines.

    ⚠ COMMENTS ARE NOT STRIPPED. None of these formats has them, so a `#` line in a relabel map is
    a defect; filtering it out would hide it and still get the count right.
    """
    return [ln for ln in s.splitlines() if ln.strip()]


# --------------------------------------------------------------------------------------------
# The content checks, as PURE FUNCTIONS of (text, identifier set).
#
# ⛔ SEPARATED FROM THE API CALLS SO THEY CAN BE PROVEN TO BITE. A checker whose logic is welded to
# a live server is a checker nobody has watched go red, and this project has shipped two guards that
# were inert. Everything below is exercised by --self-test against fixtures.
# --------------------------------------------------------------------------------------------


def check_identifier_list(txt: str, ids: set[str]) -> list[str]:
    got = set(rows(txt))
    if got == ids:
        return []
    return [f"panel_identifiers lists {len(got)} name(s), the assemblies collection has "
            f"{len(ids)}: extra {sorted(got - ids)[:4]}, missing {sorted(ids - got)[:4]}"]


def check_self_pairs(txt: str, ids: set[str]) -> list[str]:
    """One `X_X` per strain.

    ⚠ CONSTRUCTED AND COMPARED, NOT PARSED. Splitting `A_B` back into names is ambiguous once a
    strain name contains `_` (Salk_ANCa) or `.` (T2T_GCA_054642775.1_IMPD), and a parser that
    guesses the split reports the wrong strain as missing. Building what the producer should have
    written and diffing the two sets has no such failure mode.
    """
    r = rows(txt)
    want, got = {f"{i}_{i}" for i in ids}, set(r)
    out = []
    if len(r) != len(got):
        out.append(f"self_pairs has {len(r) - len(got)} duplicate row(s)")
    if got != want:
        out.append(f"self_pairs is not exactly one X_X per strain: "
                   f"extra {sorted(got - want)[:3]}, missing {sorted(want - got)[:3]}")
    return out


def check_relabel_map(txt: str, ids: set[str]) -> list[str]:
    """`A_B<TAB>A.B` for every ordered pair of DISTINCT identifiers. Same construction argument."""
    r = rows(txt)
    want = {f"{a}_{b}": f"{a}.{b}" for a in ids for b in ids if a != b}
    diag = {f"{i}_{i}" for i in ids}
    got, dup, malformed = {}, [], []
    for ln in r:
        parts = ln.split("\t")
        if len(parts) != 2:
            malformed.append(ln)
            continue
        if parts[0] in got:
            dup.append(parts[0])
        got[parts[0]] = parts[1]

    out = []
    if malformed:
        out.append(f"relabel_map has {len(malformed)} row(s) that are not two tab-separated "
                   f"fields, e.g. {malformed[:2]}")
    if dup:
        out.append(f"relabel_map has {len(dup)} duplicate key(s), e.g. {dup[:3]}")
    # ⛔ THE DIAGONAL IS CALLED OUT SEPARATELY, because it is the one failure mode with a history
    # here: A_A rows make the map longer than the collection, and WF-C's strict relabel then
    # refuses with a message about the COUNT, which sends the reader to the collection.
    present_diag = sorted(k for k in got if k in diag)
    if present_diag:
        out.append(f"relabel_map contains {len(present_diag)} A_A diagonal row(s), e.g. "
                   f"{present_diag[:2]}; WF-C strict mode refuses a map longer than the "
                   f"collection it relabels")
    wrong = sorted(k for k, v in got.items() if k in want and v != want[k])
    if wrong:
        out.append(f"{len(wrong)} relabel_map row(s) map to the wrong label, e.g. "
                   f"{[(k, got[k], want[k]) for k in wrong[:2]]}")
    absent = sorted(set(want) - set(got))
    strays = sorted(set(got) - set(want) - set(present_diag))
    if absent or strays:
        out.append(f"relabel_map has {len(got)} distinct key(s), expected n(n-1) = {len(want)} "
                   f"for n = {len(ids)}: missing {absent[:3]}, unexpected {strays[:3]}")
    return out


def check_similarity_matrix(txt: str, n: int) -> list[str]:
    """`sourmash compare --csv`: one header row of names, then one row per genome."""
    r = rows(txt)
    cols = len(r[0].split(",")) if r else 0
    if cols == n and len(r) - 1 == n:
        return []
    return [f"similarity_matrix is {max(len(r) - 1, 0)} x {cols}, expected {n} x {n}"]


def check_sizes_fai(sizes_txt: str, fai_txt: str, ident: str) -> list[str]:
    a, b = rows(sizes_txt), rows(fai_txt)
    if not a:
        return [f"{ident}: sizes is empty"]
    if len(a) != len(b):
        return [f"{ident}: {len(a)} sizes row(s) vs {len(b)} .fai row(s)"]
    for x in a:
        f = x.split("\t")
        if len(f) < 2 or not f[1].isdigit() or int(f[1]) <= 0:
            return [f"{ident}: a sizes row has no positive length ({x!r})"]
    return []


def check_signature(txt: str, ident: str) -> list[str]:
    """⛔ ONE JOB SKETCHES THE WHOLE PANEL, so a member with no hashes is not a missing element --
    it is a signature file that exists, is valid JSON, and describes nothing. 180 kb of N sketches
    to zero mins and the matrix that comes out is not an error."""
    try:
        sig = json.loads(txt)
    except json.JSONDecodeError as exc:
        return [f"{ident}: not JSON ({exc.msg})"]
    recs = sig if isinstance(sig, list) else [sig]
    mins = sum(len(s.get("mins") or [])
               for rec in recs for s in (rec.get("signatures") or []))
    return [] if mins else [f"{ident}: 0 hashes"]


# --------------------------------------------------------------------------------------------
# Server side
# --------------------------------------------------------------------------------------------


def elements(gi: GalaxyInstance, cid: str) -> dict[str, str]:
    d = gi.dataset_collections.show_dataset_collection(cid)
    return {e["element_identifier"]: e["object"]["id"] for e in d["elements"]}


def text(gi: GalaxyInstance, did: str) -> str:
    b = gi.datasets.download_dataset(did, use_default_filename=False)
    return b.decode("utf-8", "replace") if isinstance(b, bytes) else str(b)


def load_panel(path: pathlib.Path) -> dict:
    d = yaml.safe_load(path.read_text(encoding="utf-8"))
    for k in ("expected_assemblies", "proteomes", "anchors"):
        if k not in d:
            sys.exit(f"{path}: panel has no `{k}`")
    return d


def check(gi: GalaxyInstance, inv_id: str, panel: dict, expect: dict[str, int]) -> int:
    from softmask_lib import JOBS_UNFINISHED, SCHEDULING_IN_PROGRESS

    inv = gi.invocations.show_invocation(inv_id)
    if inv["state"] in SCHEDULING_IN_PROGRESS:
        sys.exit(f"invocation {inv_id} is still scheduling ({inv['state']}); nothing to verify yet")

    # ⛔ `state == "scheduled"` MEANS THE STEPS WERE CREATED, NOT THAT THEY RAN. Every output
    # dataset exists from the moment its job is queued, and downloading a queued dataset returns
    # nothing -- which lands in the count checks below as "0 rows", i.e. as a defect report about a
    # run that has not finished, or worse passes a check comparing two empty things.
    states = gi.invocations.get_invocation_summary(inv_id).get("states", {})
    pending = {k: v for k, v in states.items() if k in JOBS_UNFINISHED}
    if pending:
        sys.exit(f"invocation {inv_id} still has unfinished job(s) {pending} out of {states}. "
                 f"Verifying now would read empty outputs; wait for the run to finish.")
    # ⚠ A FAILED JOB IS NOT A REASON TO STOP. It is reported and the output checks still run:
    # which outputs the failure actually damaged is the question worth answering, and "some job
    # failed" does not answer it.
    failed = {k: v for k, v in states.items() if k not in ("ok", *JOBS_UNFINISHED)}
    head = f"⚠ job state(s) not `ok`: {failed} -- the checks below say what that cost" \
        if failed else None

    by_label = {v["label"]: v for v in (inv.get("inputs") or {}).values()
                if isinstance(v, dict) and v.get("label")}
    if "assemblies" not in by_label:
        sys.exit(f"invocation {inv_id} has no `assemblies` input; is this a WF-A invocation?")

    ids = set(elements(gi, by_label["assemblies"]["id"]))
    n = len(ids)
    problems: list[str] = []
    notes: list[str] = []

    # ⛔ 0 == 0 IS NOT A PASS. Every count below would agree with itself on an empty panel.
    if n == 0:
        sys.exit("the assemblies input has NO elements, so every count check below would pass "
                 "trivially. Refusing.")
    if n != panel["expected_assemblies"]:
        problems.append(f"the invocation ran on {n} assemblies, the panel asserts "
                        f"{panel['expected_assemblies']}")

    want = {"assemblies": ids,
            "proteomes": set(panel["proteomes"]),
            "anchors": set(panel["anchors"])}

    outs = dict(inv.get("output_collections") or {})
    singles = dict(inv.get("outputs") or {})
    absent = [k for k in KEYED_BY if k not in outs and k not in singles]
    if absent:
        problems.append(f"{len(absent)} declared output(s) absent from the invocation: "
                        f"{', '.join(sorted(absent))}")

    for name, keyed in KEYED_BY.items():
        if keyed is None or name not in outs:
            continue
        got = set(elements(gi, outs[name]["id"]))
        if got != want[keyed]:
            problems.append(f"{name} is keyed by `{keyed}` ({len(want[keyed])}) but holds "
                            f"{len(got)}: extra {sorted(got - want[keyed])[:4]}, "
                            f"missing {sorted(want[keyed] - got)[:4]}")
        else:
            notes.append(f"{name}: {len(got)} element(s), identifiers match `{keyed}`")

    for name, fn, ok_msg in (
            ("panel_identifiers", lambda t: check_identifier_list(t, ids),
             f"{n} name(s), matches the collection"),
            ("self_pairs", lambda t: check_self_pairs(t, ids),
             f"{n} row(s), exactly one X_X per strain"),
            ("relabel_map", lambda t: check_relabel_map(t, ids),
             f"{n * (n - 1)} row(s) = n(n-1), every ordered distinct pair once, no diagonal"),
            ("similarity_matrix", lambda t: check_similarity_matrix(t, n), f"{n} x {n}")):
        if name not in singles:
            continue
        found = fn(text(gi, singles[name]["id"]))
        problems.extend(found)
        if not found:
            notes.append(f"{name}: {ok_msg}")

    if "sizes" in outs and "fasta_index" in outs:
        sz, fai = elements(gi, outs["sizes"]["id"]), elements(gi, outs["fasta_index"]["id"])
        short: list[str] = []
        for ident in sorted(set(sz) & set(fai)):
            short.extend(check_sizes_fai(text(gi, sz[ident]), text(gi, fai[ident]), ident))
        if short:
            problems.append(f"{len(short)} strain(s) with a bad sizes/.fai pair: {short[:3]}")
        else:
            notes.append(f"sizes and .fai agree row for row on all {len(sz)} strain(s)")

    if "signatures" in outs:
        bad: list[str] = []
        for ident, did in sorted(elements(gi, outs["signatures"]["id"]).items()):
            bad.extend(check_signature(text(gi, did), ident))
        if bad:
            problems.append(f"{len(bad)} signature(s) carry no hashes: {bad[:3]}")
        else:
            notes.append("every signature carries at least one hash")

    for spec, wanted in expect.items():
        out, _, ident = spec.partition(":")
        if out not in outs:
            problems.append(f"--expect names output `{out}`, which this invocation does not have")
            continue
        els = elements(gi, outs[out]["id"])
        if ident not in els:
            problems.append(f"--expect names `{ident}` in {out}, which holds {sorted(els)[:4]}")
            continue
        got = len(rows(text(gi, els[ident])))
        if got != wanted:
            problems.append(f"CONTROL FAILED: {out}/{ident} has {got} row(s), expected {wanted}. "
                            f"This branch does not depend on panel size, so a change here is a "
                            f"regression, not growth")
        else:
            notes.append(f"control {out}/{ident}: {got} rows, unchanged")

    if head:
        print(f"  {head}")
    for line in notes:
        print(f"  {line}")
    if problems:
        print(f"\n⛔ {len(problems)} problem(s) with invocation {inv_id}:")
        for p in problems:
            print(f"  - {p}")
        return 1
    print(f"\n✅ WF-A over {n} assemblies, {len(panel['proteomes'])} proteome(s) and "
          f"{len(panel['anchors'])} anchor(s): every output is keyed by the collection it should "
          f"be, and every count that describes the panel equals the panel")
    return 0


# --------------------------------------------------------------------------------------------

#: Three strains, chosen so the ambiguous-split trap is live: `a_b` is a name, and so is `a` --
#: a parser splitting `a_b_a` on the first `_` would read the pair as (a, b_a).
_IDS = {"a", "a_b", "c.1_d"}


def _pairs(ids: set[str], *, diagonal: bool = False) -> str:
    return "".join(f"{x}_{y}\t{x}.{y}\n" for x in sorted(ids) for y in sorted(ids)
                   if diagonal or x != y)


def _sig(mins: int) -> str:
    return json.dumps([{"signatures": [{"mins": list(range(mins))}]}])


_SELF_TEST = [
    ("identifier list: correct", check_identifier_list, "a\na_b\nc.1_d\n", False),
    ("identifier list: one missing", check_identifier_list, "a\na_b\n", True),
    ("identifier list: a stranger", check_identifier_list, "a\na_b\nc.1_d\nz\n", True),
    ("self_pairs: correct", check_self_pairs, "a_a\na_b_a_b\nc.1_d_c.1_d\n", False),
    ("self_pairs: one missing", check_self_pairs, "a_a\na_b_a_b\n", True),
    ("self_pairs: duplicated", check_self_pairs,
     "a_a\na_a\na_b_a_b\nc.1_d_c.1_d\n", True),
    ("self_pairs: a cross pair", check_self_pairs,
     "a_a\na_b_a_b\nc.1_d_c.1_d\na_a_b\n", True),
    ("relabel_map: correct", check_relabel_map, _pairs(_IDS), False),
    ("relabel_map: with the A_A diagonal", check_relabel_map, _pairs(_IDS, diagonal=True), True),
    ("relabel_map: one row short", check_relabel_map,
     "".join(_pairs(_IDS).splitlines(keepends=True)[1:]), True),
    ("relabel_map: a one-column row", check_relabel_map,
     _pairs(_IDS) + "a_a_b\n", True),
    ("relabel_map: wrong label", check_relabel_map,
     _pairs(_IDS).replace("a.a_b", "a.WRONG", 1), True),
    ("relabel_map: duplicate key", check_relabel_map, _pairs(_IDS) + "a_a_b\ta.a_b\n", True),
    ("relabel_map: empty file", check_relabel_map, "", True),
]

_SELF_TEST_N = [
    ("similarity_matrix: 3 x 3", check_similarity_matrix, "a,b,c\n1,0,0\n0,1,0\n0,0,1\n", False),
    ("similarity_matrix: a row short", check_similarity_matrix, "a,b,c\n1,0,0\n0,1,0\n", True),
    ("similarity_matrix: a column short", check_similarity_matrix,
     "a,b\n1,0\n0,1\n0,0\n", True),
    ("similarity_matrix: empty", check_similarity_matrix, "", True),
]


def self_test() -> int:
    """Break each property separately and confirm it is reported."""
    bad = 0
    for label, fn, txt, should_fail in _SELF_TEST:
        found = fn(txt, set(_IDS))
        ok = bool(found) == should_fail
        bad += not ok
        print(f"  {'pass' if ok else '⛔ FAIL'}  {label:44} "
              f"{found[0][:70] if found else 'clean'}")
    for label, fn, txt, should_fail in _SELF_TEST_N:
        found = fn(txt, 3)
        ok = bool(found) == should_fail
        bad += not ok
        print(f"  {'pass' if ok else '⛔ FAIL'}  {label:44} "
              f"{found[0][:70] if found else 'clean'}")

    extra = [
        ("sizes/fai: agree", check_sizes_fai("c1\t100\nc2\t50\n", "c1\t100\t1\t60\t61\n"
                                             "c2\t50\t1\t60\t61\n", "s"), False),
        ("sizes/fai: fai one row short", check_sizes_fai("c1\t100\nc2\t50\n",
                                                         "c1\t100\t1\t60\t61\n", "s"), True),
        ("sizes: empty is NOT clean", check_sizes_fai("", "", "s"), True),
        ("sizes: zero length", check_sizes_fai("c1\t0\n", "c1\t0\t1\t60\t61\n", "s"), True),
        ("sizes: length not a number", check_sizes_fai("c1\tNA\n", "c1\tNA\t1\t60\t61\n", "s"),
         True),
        ("signature: has hashes", check_signature(_sig(5), "s"), False),
        ("signature: zero hashes", check_signature(_sig(0), "s"), True),
        ("signature: not JSON", check_signature("not json at all", "s"), True),
        ("signature: JSON but no signatures", check_signature("[{}]", "s"), True),
    ]
    for label, found, should_fail in extra:
        ok = bool(found) == should_fail
        bad += not ok
        print(f"  {'pass' if ok else '⛔ FAIL'}  {label:44} "
              f"{found[0][:70] if found else 'clean'}")

    n = len(_SELF_TEST) + len(_SELF_TEST_N) + len(extra)
    print(f"\n  {n - bad}/{n} self-test case(s) behaved correctly")
    return 1 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--invocation")
    ap.add_argument("--panel", type=pathlib.Path,
                    help="the panel the run was staged from -- the independent second route")
    ap.add_argument("--expect", action="append", default=[], metavar="OUT:IDENT=ROWS",
                    help="a control row count from an earlier run, e.g. "
                         "anchor_bed12s:cs10=40185. Repeatable.")
    ap.add_argument("--self-test", action="store_true", help="prove each check bites")
    a = ap.parse_args()
    if a.self_test:
        return self_test()
    if not (a.invocation and a.panel):
        ap.error("--invocation and --panel are both required (or use --self-test)")
    expect = {}
    for spec in a.expect:
        key, _, val = spec.rpartition("=")
        if not key or not val.isdigit():
            sys.exit(f"--expect {spec!r} is not OUT:IDENT=ROWS")
        expect[key] = int(val)
    from softmask_lib import connect

    return check(connect(), a.invocation, load_panel(a.panel), expect)


if __name__ == "__main__":
    sys.exit(main())
