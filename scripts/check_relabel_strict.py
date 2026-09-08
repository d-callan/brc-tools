#!/usr/bin/env python3
"""Hold together the three facts that let WF-C detect a wrong relabel map.

WF-C renames its cross-product cells from `A_B` to the `A.B` ids Phase E parses, using Galaxy's
`__RELABEL_FROM_FILE__` over a 2-column map from `collection_relabel_map`. That step is the ONLY
place a mismatched map can be caught, and only in strict mode:

    lib/galaxy/tools/__init__.py:5013   default = None if strict else element_identifier

With strict off, an element absent from the map keeps its own `A_B` identifier and the step still
reports success -- so a map built from the wrong collection yields a green run whose outputs Phase E
silently cannot key. Both editions ran that way until this was written.

⛔ STRICT MODE ALSO CHECKS THE ROW COUNT (line 4983), AND THAT IS WHY THREE THINGS ARE ONE FACT.
The relabel step sits downstream of the four self-pair filters, so the collection holds n**2-n
elements. `relabel_map` therefore has to OMIT the `A_A` diagonal, or the count check refuses every
run before it looks at a single identifier. Turning strict on and dropping the diagonal are not two
improvements; either alone breaks the workflow. This script fails if they ever drift apart:

  * `tools/collection_relabel_map/relabel_map.py` emits no diagonal row,
  * both editions set `strict: true` on every `__RELABEL_FROM_FILE__` step,
  * and Galaxy's own logic, replayed here, accepts the real map and rejects the three ways it can
    be wrong.

⚠ THE REPLAY IS TRANSCRIBED, NOT IMPORTED. Checking this must not require a Galaxy install (see
check_udt_provenance.py's note on the same choice), so `galaxy_relabel` below is a transcription of
produce_outputs from 26.1 and carries the line numbers to re-read when Galaxy changes.

    python3 scripts/check_relabel_strict.py          # non-zero if any of the three has drifted
"""

from __future__ import annotations

import itertools
import pathlib
import subprocess
import sys
import tempfile

import yaml

ROOT = pathlib.Path(__file__).resolve().parents[1]
HELPER = ROOT / "tools/collection_relabel_map/relabel_map.py"
EDITIONS = ("align_chain.gxwf.yml", "align_chain_udt.gxwf.yml")
#: Six real panel members, because the identifiers themselves are part of what is being checked:
#: two carry a `.` (so "count the dots" is not a usable test for a renamed id) and most carry `_`
#: (so neither is "split on the underscore").
IDS = ["cs10_NCBI_RefSeq_softmasked", "ASM2916894v1", "T2T_GCA_054642775.1_IMPD",
       "T2T_GCA_054642815.1_IMPD", "JL_Father", "JL_Mother"]


def galaxy_relabel(elements: list[str], rows: list[str], strict: bool):
    """RelabelFromFileTool.produce_outputs, tabular branch, transcribed from 26.1.

    Returns the new identifiers, or the text of the MessageException Galaxy would raise.
    """
    if strict and len(elements) != len(rows):                       # __init__.py:4983
        return "ERROR: Relabel mapping file contains incorrect number of identifiers"
    mapping = {}
    for i, line in enumerate(rows, 1):
        cols = line.split("\t")
        if len(cols) != 2:                                          # __init__.py:4999
            return f"ERROR: Relabel mapping file contains {len(cols)} columns on line {i}"
        mapping[cols[0]] = cols[1]
    out = []
    for el in elements:
        new = mapping.get(el, None if strict else el)               # __init__.py:5013
        if not new:
            return f"ERROR: Failed to find original identifier [{el}]"
        out.append(new)
    return out


def emitted_map(ids: list[str]) -> list[str]:
    """Run the real helper, so this checks the shipped tool and not a paraphrase of it."""
    with tempfile.TemporaryDirectory() as td:
        wd = pathlib.Path(td)
        (wd / "ids.txt").write_text("".join(f"{i}\n" for i in ids))
        subprocess.run([sys.executable, str(HELPER), "--ids-file", str(wd / "ids.txt"),
                        "--out", str(wd / "map.tabular")], check=True)
        return [ln for ln in (wd / "map.tabular").read_text().splitlines() if ln]


def strict_flags() -> list[tuple[str, str, object]]:
    """(edition, step name, strict value) for every relabel step in both workflows."""
    found = []
    for name in EDITIONS:
        doc = yaml.safe_load((ROOT / "workflows/align_chain_project" / name).read_text())
        for step, body in doc["steps"].items():
            if body.get("tool_id") == "__RELABEL_FROM_FILE__":
                how = (body.get("state") or {}).get("how") or {}
                found.append((name, step, how.get("strict", "(unset)")))
    return found


def main() -> int:
    bad = []

    rows = emitted_map(IDS)
    # The collection the step actually receives: the flat cross product minus the diagonal that
    # the self_pairs list removes.
    cells = [f"{a}_{b}" for a, b in itertools.product(IDS, IDS) if a != b]
    print(f"n={len(IDS)}: collection {len(cells)} elements, map {len(rows)} rows")
    if len(rows) != len(cells):
        bad.append(f"relabel_map emits {len(rows)} rows for {len(cells)} elements. Strict mode "
                   f"checks the count first, so this refuses every run. Does it still skip the "
                   f"`A_A` diagonal?")

    steps = strict_flags()
    for edition, step, value in steps:
        if value is not True:
            bad.append(f"{edition}: step `{step}` has strict={value!r}. With strict off, an "
                       f"element missing from the map keeps its `A_B` identifier and the step "
                       f"still succeeds -- the mislabel this pairing exists to catch.")
    print(f"strict flags: {len(steps)} relabel step(s) across {len(EDITIONS)} edition(s), "
          f"{sum(1 for _, _, v in steps if v is True)} strict")

    # ⛔ THE FOUR CASES, AND WHY THE LAST TWO ARE HERE. A checker that only confirmed the happy
    # path would pass just as well against a map with the diagonal (the count check would refuse
    # at RUN time instead), and would say nothing about whether strict actually detects anything.
    cases = [
        ("correct map, strict", cells, rows, True, "list"),
        ("map WITH the diagonal, strict",
         cells, [f"{a}_{b}\t{a}.{b}" for a, b in itertools.product(IDS, IDS)], True, "error"),
        ("map from the WRONG collection (right row count), strict",
         cells, [f"{a}_{b}\t{a}.{b}" for a, b in
                 itertools.product([*IDS[:-1], "OTHER"], repeat=2) if a != b], True, "error"),
        ("map from the WRONG collection, strict OFF",
         cells, [f"{a}_{b}\t{a}.{b}" for a, b in
                 itertools.product([*IDS[:-1], "OTHER"], repeat=2) if a != b], False, "silent"),
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
    one = emitted_map(["cs10"])
    r = galaxy_relabel([], one, strict=True)
    print(f"  {'pass' if one == [] and r == [] else '⛔ FAIL'}  n=1: {len(one)} rows against 0 "
          f"elements -> {r if isinstance(r, str) else 'accepted'}")
    if one != [] or not isinstance(r, list):
        bad.append(f"n=1 should yield an empty map that strict mode accepts; got {len(one)} "
                   f"row(s) and {r}")

    for line in bad:
        print(f"⛔ {line}")
    print(f"\n{len(bad)} problem(s)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
