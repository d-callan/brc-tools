#!/usr/bin/env python3
"""Hold together the facts that let WF-C detect a wrong relabel map.

WF-C renames its cross-product cells from `A_B` to the `A.B` ids Phase E parses, using Galaxy's
`__RELABEL_FROM_FILE__` over a 2-column map. That step is the ONLY place a mismatched map can be
caught, and only in strict mode:

    lib/galaxy/tools/__init__.py:5013   default = None if strict else element_identifier

With strict off, an element absent from the map keeps its own `A_B` identifier and the step still
reports success -- so a map built from the wrong collection yields a green run whose outputs Phase E
silently cannot key. Both editions ran that way until this was written.

⛔ STRICT MODE ALSO CHECKS THE ROW COUNT (line 4983), AND THAT IS WHY SEVERAL THINGS ARE ONE FACT.
The relabel steps sit downstream of the four self-pair filters, so the collection holds n**2-n
elements. Every producer of the map therefore has to OMIT the `A_A` diagonal, or the count check
refuses the run before it looks at a single identifier. Turning strict on and dropping the diagonal
are not separate improvements; either alone breaks the workflow.

⛔ AND THERE IS MORE THAN ONE PRODUCER, WHICH IS HOW THIS WENT WRONG. The first version of this
script checked only `tools/collection_relabel_map/relabel_map.py` and printed "0 problems" while
`execution/cluster/gen_wfc_config.py` -- the generator the documented one-click recipe calls -- was
still emitting n**2 rows, so every cluster run would have been refused. A guard over one of two
producers of the same file is worse than none, because it reads as coverage. Both are checked here,
and PRODUCERS is the list to extend when a third appears.

⚠ THE REPLAY IS TRANSCRIBED, NOT IMPORTED, SO IT MUST BE AS STRICT AS GALAXY AND NOT MERELY
SIMILAR. Checking this must not require a Galaxy install (see check_udt_provenance.py's note on the
same choice). But the first transcription reproduced only the count, column and lookup checks and
dropped `add_copied_value_to_new_elements` -- Galaxy's `.strip()`, its identifier regex and its
duplicate check. It therefore ACCEPTED a map Galaxy rejects: with `GCA_000001405.29+alt` in the
panel it reported "6 relabelled" while Galaxy refuses `GCA_000001405.29+alt.cs10` on
`^[\\w\\- \\.,]+$`. A `+`, `:`, `#` or `/` in an assembly name is a legal Galaxy element identifier
and common in accessions, so that was a confident pass for a step that fails at run time. The
transcription now carries those three checks and the line numbers to re-read when Galaxy changes.

    python3 scripts/check_relabel_strict.py             # non-zero if anything has drifted
    python3 scripts/check_relabel_strict.py --self-test # prove the checks are not inert
"""

from __future__ import annotations

import argparse
import itertools
import pathlib
import re
import subprocess
import sys
import tempfile

import yaml

ROOT = pathlib.Path(__file__).resolve().parents[1]
WF_DIR = ROOT / "workflows"

#: Every producer of WF-C's `relabel_map`, and how to run it. Each must emit one row per ORDERED
#: DISTINCT pair -- n**2-n -- because strict mode compares that count against the collection.
PRODUCERS = ("collection_relabel_map", "gen_wfc_config")

#: WF-C2's producer, checked on a DIFFERENT invariant. `collection_anchor_grid` does not emit a
#: cross product, so "n**2-n rows" says nothing about it. What strict mode needs there is that the
#: `keep` list it writes and the `relabel` map it writes describe the same cells in the same order,
#: because `keep` is what selects the elements the map is then applied to.
GRID_PRODUCER = "collection_anchor_grid"

#: The workflows whose relabel steps this script verifies, and how many each must have. Asserted,
#: because a guard that finds NOTHING must not pass: stubbing the scan to return [] once printed
#: "0 relabel step(s) ... 0 problems" and exited 0, which is what a renamed workflow, a `steps:`
#: list instead of a mapping, or a deleted step would have produced.
#:
#: ⚠ WF-C2 IS HERE ON A MEASUREMENT, NOT BY ANALOGY. `project_annotations.gxwf.yml`'s `p_chain` is
#: fed by `collection_anchor_grid` rather than by a cross product, so it was reported as unverified
#: until the count was actually checked: the tool writes `keep` and `relabel` from ONE pairs list,
#: `p_chain_keep` filters with `remove_if_absent` against that same `keep`, and the survivors
#: therefore equal the rows -- 21 and 21 for the documented 3-anchor, 8-strain panel. Strict
#: refuses only when WF-C left a grid cell unproduced (20 against 21), which is the loud failure
#: one wants: without it the projection silently runs on an incomplete grid. `project_annotations_
#: udt.gxwf.yml` is absent from this list because it has no relabel step at all -- a WF-C2 parity
#: question, not a strict-mode one.
VERIFIED_WORKFLOWS = {
    "align_chain.gxwf.yml": 2,
    "align_chain_udt.gxwf.yml": 2,
    "project_annotations.gxwf.yml": 1,
}
#: The only mode the replay below describes. `txt` relabels BY LINE ORDER with no key lookup
#: (__init__.py:5019-5026), which would be catastrophic against a map keyed on `A_B`, and
#: `tabular_extended` honours `from`/`to` columns this hardcodes to 1 and 2. `strict` would still
#: read true in either case, so the mode has to be asserted or the whole replay stops applying.
REQUIRED_HOW_SELECT = "tabular"

#: Six real panel members, because the identifiers themselves are part of what is being checked:
#: two carry a `.` (so "count the dots" is not a usable test for a renamed id) and most carry `_`
#: (so neither is "split on the underscore").
IDS = ["cs10_NCBI_RefSeq_softmasked", "ASM2916894v1", "T2T_GCA_054642775.1_IMPD",
       "T2T_GCA_054642815.1_IMPD", "JL_Father", "JL_Mother"]

#: Galaxy's own identifier rule, __init__.py:4968. Kept as the literal pattern so a diff against
#: Galaxy is a text comparison.
IDENTIFIER_RE = r"^[\w\- \.,]+$"


def galaxy_relabel(elements: list[str], rows: list[str], strict: bool):
    """RelabelFromFileTool.produce_outputs, tabular branch, transcribed from 26.1.

    Returns the new identifiers, or the text of the MessageException Galaxy would raise. The line
    numbers are the point: this is a copy, and a copy has to be re-checked when the original moves.
    """
    if strict and len(elements) != len(rows):                       # :4983
        return "ERROR: Relabel mapping file contains incorrect number of identifiers"
    mapping = {}
    for i, line in enumerate(rows, 1):
        cols = line.strip().split("\t")                             # :4999 (line.strip() first)
        if len(cols) != 2:                                          # :5001
            return f"ERROR: Relabel mapping file contains {len(cols)} columns on line {i}"
        mapping[cols[0]] = cols[1]
    out: list[str] = []
    seen: set[str] = set()
    for el in elements:
        new = mapping.get(el, None if strict else el)               # :5013
        if not new:
            return f"ERROR: Failed to find original identifier [{el}]"
        # ⛔ add_copied_value_to_new_elements, :4966-4973. Dropping these three made the replay
        # more permissive than Galaxy, which is the one direction a guard must never fail in.
        new = new.strip()
        if not re.match(IDENTIFIER_RE, new):
            return f"ERROR: Invalid new collection identifier [{new}]"
        if new in seen:
            return (f"ERROR: New identifier [{new}] appears twice in resulting collection, these "
                    f"values must be unique.")
        seen.add(new)
        out.append(new)
    return out


def emitted_map(producer: str, ids: list[str]) -> list[str]:
    """Run a real producer, so this checks shipped code and not a paraphrase of it.

    ⚠ BLANK LINES ARE KEPT. Galaxy counts rows with `readlines()` (`_read_text_file_lines`,
    :4141), so a trailing blank line counts toward the strict row total. Filtering them here would
    hide a producer that emits one.
    """
    with tempfile.TemporaryDirectory() as td:
        wd = pathlib.Path(td)
        if producer == "collection_relabel_map":
            (wd / "ids.txt").write_text("".join(f"{i}\n" for i in ids))
            subprocess.run([sys.executable,
                            str(ROOT / "tools/collection_relabel_map/relabel_map.py"),
                            "--ids-file", str(wd / "ids.txt"),
                            "--out", str(wd / "map.tabular")], check=True)
            text = (wd / "map.tabular").read_text()
        elif producer == "gen_wfc_config":
            # It takes anchors too; any non-empty subset will do, the map is over strains.
            subprocess.run([sys.executable, str(ROOT / "execution/cluster/gen_wfc_config.py"),
                            str(wd / "cfg"), ",".join(ids), ids[0]], check=True,
                           capture_output=True)
            text = (wd / "cfg/relabel_map.tsv").read_text()
        else:
            raise SystemExit(f"unknown producer {producer!r}")
    return text.splitlines()


def check_grid_producer() -> list[str]:
    """Run collection_anchor_grid and hold its two outputs to what strict mode needs of them.

    ⚠ THE INVARIANT IS NOT A ROW COUNT AGAINST n**2-n. It is that `keep` and `relabel` agree: the
    filter selects on `keep`, the map is applied to whatever survived, and Galaxy compares those
    two counts. Checking a cross-product formula here would have been a check of the wrong thing
    that happened to pass.
    """
    bad = []
    anchors = ["PvW1", "PAM", "PvSY56"]
    strains = ["PvP01", "PvW1", "PAM", "PvSY56", "Sal-I", "PvT01", "PvC01", "MHC087"]
    with tempfile.TemporaryDirectory() as td:
        wd = pathlib.Path(td)
        subprocess.run([sys.executable, str(ROOT / "tools/collection_anchor_grid/anchor_grid.py"),
                        "--anchors", " ".join(anchors), "--strains", " ".join(strains),
                        "--keep", str(wd / "keep.txt"), "--relabel", str(wd / "relabel.tsv"),
                        "--order", str(wd / "order.txt")], check=True, capture_output=True)
        keep = (wd / "keep.txt").read_text().splitlines()
        rows = (wd / "relabel.tsv").read_text().splitlines()
    want = len(anchors) * (len(strains) - 1)
    print(f"  {GRID_PRODUCER:24} {len(rows):>4} rows against {len(keep)} kept cells "
          f"(|anchors| x (n-1) = {want})")
    if len(rows) != len(keep):
        bad.append(f"{GRID_PRODUCER} writes {len(rows)} relabel rows for {len(keep)} kept cells. "
                   f"`remove_if_absent` selects on `keep`, so strict mode compares those two and "
                   f"refuses the step.")
    if len(keep) != want:
        bad.append(f"{GRID_PRODUCER} kept {len(keep)} cells, expected {want} -- the grid is "
                   f"|anchors| x (n-1) once the anchor self-cells are dropped.")
    first_col = [r.split("\t")[0] for r in rows]
    if first_col != keep:
        bad.append(f"{GRID_PRODUCER}'s relabel keys do not match its `keep` list. Every kept "
                   f"element must have a row, or strict mode fails the lookup: "
                   f"{sorted(set(keep) ^ set(first_col))[:3]}")
    # And the replay must accept the real pair, and refuse an incomplete upstream.
    got = galaxy_relabel(keep, rows, strict=True)
    if not isinstance(got, list) or len(got) != len(keep):
        bad.append(f"{GRID_PRODUCER}: strict mode would refuse its own complete output: {got}")
    short = galaxy_relabel(keep[:-1], rows, strict=True)
    if not isinstance(short, str):
        bad.append(f"{GRID_PRODUCER}: strict mode accepted a map with one element MISSING "
                   f"upstream, which is the incomplete grid it exists to refuse.")
    print(f"  {'pass' if not bad else '⛔ FAIL'}  {GRID_PRODUCER}: complete grid accepted "
          f"({len(keep)} cells), a missing upstream chain refused")
    return bad


def relabel_steps() -> list[tuple[str, str, dict]]:
    """(workflow file name, step name, the step's `how` state) for EVERY relabel step in the repo.

    ⚠ GLOBBED, NOT LISTED. A hardcoded pair of paths made the scan blind to a renamed or moved
    workflow -- and blind to `project_annotations.gxwf.yml`'s relabel step, which sat at
    `strict: false` unseen while this script reported everything fine.
    """
    found = []
    for wf in sorted(WF_DIR.glob("*/*.gxwf.yml")):
        doc = yaml.safe_load(wf.read_text())
        steps = doc.get("steps")
        if not isinstance(steps, dict):
            # gxformat2 also permits a LIST of steps; nothing here uses it, but silently reading
            # zero steps out of a list is how this guard would go quietly green.
            raise SystemExit(f"{wf.name} declares `steps` as {type(steps).__name__}, which this "
                             f"scan does not understand -- teach it the list form rather than "
                             f"letting it find no relabel steps and pass.")
        for step, body in steps.items():
            if (body or {}).get("tool_id") == "__RELABEL_FROM_FILE__":
                found.append((wf.name, step, ((body.get("state") or {}).get("how") or {})))
    return found


def check() -> list[str]:
    bad: list[str] = []
    cells = [f"{a}_{b}" for a, b in itertools.product(IDS, IDS) if a != b]

    maps = {}
    for producer in PRODUCERS:
        rows = emitted_map(producer, IDS)
        maps[producer] = rows
        print(f"  {producer:24} {len(rows):>4} rows for {len(IDS)} ids "
              f"(collection has {len(cells)})")
        if len(rows) != len(cells):
            bad.append(f"{producer} emits {len(rows)} rows for {len(cells)} elements. Strict mode "
                       f"checks the count first, so this refuses every run. Does it still skip "
                       f"the `A_A` diagonal?")
    # ⚠ The producers must agree with EACH OTHER too, not just with the count -- two files that
    # are both the right length but disagree on the id format would each pass alone.
    if len({tuple(sorted(v)) for v in maps.values()}) > 1:
        bad.append(f"the producers disagree on the map's content: {list(maps)}. They feed the same "
                   f"workflow input, so a run's behaviour would depend on which one made the file.")

    bad += check_grid_producer()

    steps = relabel_steps()
    verified = [s for s in steps if s[0] in VERIFIED_WORKFLOWS]
    other = [s for s in steps if s[0] not in VERIFIED_WORKFLOWS]
    print(f"  relabel steps: {len(verified)} in the verified workflow(s), {len(other)} elsewhere")
    for wf, want in VERIFIED_WORKFLOWS.items():
        got = sum(1 for w, _, _ in verified if w == wf)
        if got != want:
            bad.append(f"found {got} relabel step(s) in {wf}, expected {want}. Finding none is "
                       f"not a pass -- a renamed workflow or a deleted step would leave this "
                       f"guard green while the property it checks is gone.")
    for wf, step, how in verified:
        if how.get("strict") is not True:
            bad.append(f"{wf}: step `{step}` has strict={how.get('strict', '(unset)')!r}. With "
                       f"strict off, an element missing from the map keeps its `A_B` identifier "
                       f"and the step still succeeds. Galaxy's own default is False, so this has "
                       f"to stay written down.")
        if how.get("how_select") != REQUIRED_HOW_SELECT:
            bad.append(f"{wf}: step `{step}` has how_select="
                       f"{how.get('how_select', '(unset)')!r}, not {REQUIRED_HOW_SELECT!r}. `txt` "
                       f"relabels by LINE ORDER with no key lookup, and `tabular_extended` reads "
                       f"columns this check assumes are 1 and 2 -- strict would still be true "
                       f"while the map was being applied to the wrong elements.")
    for wf, step, how in other:
        print(f"  ⚠ NOT VERIFIED HERE  {wf}: step `{step}` has strict="
              f"{how.get('strict', '(unset)')!r}. Its map does not come from a producer this "
              f"script knows, so whether strict is safe is a separate question -- reported so "
              f"that it is visible, not passed.")

    rows = maps[PRODUCERS[0]]
    # ⛔ THE CASES, AND WHY THE LAST THREE ARE HERE. A checker that only confirmed the happy path
    # would pass just as well against a map with the diagonal, and would say nothing about whether
    # strict detects anything at all.
    wrong = [f"{a}_{b}\t{a}.{b}" for a, b in
             itertools.product([*IDS[:-1], "OTHER"], repeat=2) if a != b]
    cases = [
        ("correct map, strict", cells, rows, True, "list"),
        ("map WITH the diagonal, strict",
         cells, [f"{a}_{b}\t{a}.{b}" for a, b in itertools.product(IDS, IDS)], True, "error"),
        ("map from the WRONG collection (right row count), strict", cells, wrong, True, "error"),
        ("map from the WRONG collection, strict OFF", cells, wrong, False, "silent"),
        # The three add_copied_value_to_new_elements checks, each on its own.
        ("a new identifier Galaxy's regex rejects, strict",
         ["a_b"], ["a_b\tGCA_1+alt.cs10"], True, "error"),
        ("two rows mapping to the SAME new identifier, strict",
         ["a_b", "c_d"], ["a_b\tsame.id", "c_d\tsame.id"], True, "error"),
        ("a three-column row, strict", ["a_b"], ["a_b\tx.y\textra"], True, "error"),
    ]
    for label, els, mp, strict, expect in cases:
        got = galaxy_relabel(els, mp, strict)
        if expect == "list":
            ok = isinstance(got, list) and len(got) == len(els) and all("." in x for x in got)
            note = f"{len(got)} relabelled" if isinstance(got, list) else got
        elif expect == "error":
            ok = isinstance(got, str)
            note = got if isinstance(got, str) else f"NO ERROR -- {len(got)} relabelled"
        else:
            kept = [x for x in got if x in els] if isinstance(got, list) else []
            ok = isinstance(got, list) and bool(kept)
            note = f"{len(kept)} element(s) kept their `A_B` name with no error"
        print(f"  {'pass' if ok else '⛔ FAIL'}  {label}: {note}")
        if not ok:
            bad.append(f"{label}: expected {expect}, got {note}")

    # ⚠ A SINGLE IDENTIFIER MUST STAY AN HONEST ZERO, not a refusal. Its only cross-product cell is
    # its own diagonal, which WF-C filters out, so an empty map against an empty collection is the
    # right answer -- and this is the one place a 0-row output here is not a silent failure.
    one = [ln for ln in emitted_map(PRODUCERS[0], ["cs10"]) if ln]
    r = galaxy_relabel([], one, strict=True)
    ok = one == [] and isinstance(r, list)
    print(f"  {'pass' if ok else '⛔ FAIL'}  n=1: {len(one)} rows against 0 elements -> "
          f"{r if isinstance(r, str) else 'accepted'}")
    if not ok:
        bad.append(f"n=1 should yield an empty map that strict mode accepts; got {len(one)} "
                   f"row(s) and {r}")
    return bad


def self_test() -> int:
    """Prove the replay actually rejects each shape Galaxy rejects, and accepts what it accepts.

    ⚠ THIS IS THE CHECK THE FIRST VERSION LACKED. Its cases all exercised the happy path plus the
    count and lookup errors, so removing Galaxy's identifier regex from the transcription changed
    nothing that was tested -- and the regex was in fact absent.
    """
    cases = [
        ("plain map accepted", ["a_b"], ["a_b\tx.y"], True, True),
        ("`+` in the new id rejected (Galaxy's regex)", ["a_b"], ["a_b\tGCA_1+alt.cs10"], True,
         False),
        ("`:` in the new id rejected", ["a_b"], ["a_b\tx:y"], True, False),
        ("`/` in the new id rejected", ["a_b"], ["a_b\tx/y"], True, False),
        ("dot, dash, comma and space accepted", ["a_b"], ["a_b\tx-y .z,w"], True, True),
        ("duplicate new id rejected", ["a_b", "c_d"], ["a_b\tsame", "c_d\tsame"], True, False),
        ("three columns rejected", ["a_b"], ["a_b\tx\ty"], True, False),
        ("row count mismatch rejected", ["a_b"], ["a_b\tx.y", "c_d\tz.w"], True, False),
        ("missing element rejected under strict", ["a_b"], ["c_d\tz.w"], True, False),
        ("missing element KEPT when strict is off", ["a_b"], ["c_d\tz.w"], False, True),
        ("blank row counts as a row (Galaxy readlines)", ["a_b"], ["a_b\tx.y", ""], True, False),
    ]
    failed = 0
    for label, els, rows, strict, expect_ok in cases:
        got = galaxy_relabel(els, rows, strict)
        ok = isinstance(got, list) == expect_ok
        failed += not ok
        print(f"  {'pass' if ok else '⛔ FAIL'}  {label}: "
              f"{'accepted' if isinstance(got, list) else got}")
    print(f"\n  {len(cases) - failed}/{len(cases)} self-test case(s) behaved correctly")
    return 1 if failed else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--self-test", action="store_true",
                    help="check the transcription against the shapes Galaxy accepts and rejects")
    args = ap.parse_args()
    if args.self_test:
        return self_test()
    bad = check()
    for line in bad:
        print(f"⛔ {line}")
    print(f"\n{len(bad)} problem(s)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
