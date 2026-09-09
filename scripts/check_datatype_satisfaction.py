#!/usr/bin/env python3
"""Every workflow edge: does the producer's datatype satisfy the consumer's declared format?

⛔ DATATYPE SATISFACTION IS ONE-DIRECTIONAL, AND THAT IS THE WHOLE REASON THIS EXISTS. Galaxy
accepts a dataset for an input when the dataset's datatype CLASS is the input's class or a SUBCLASS
of it. So a subtype flows into a supertype input and never the reverse, and there is usually no
converter to rescue it:

    sourmash.sig IS-A json   -> True     a signature satisfies a `json` input
    json IS-A sourmash.sig   -> False    a `json` dataset does NOT satisfy a signature input
    converters between them  -> none

That asymmetry cost this repository a real defect. WF-A published its `signatures` collection --
the output documented as BRC-reusable, the one that "closes WF-A's one parity gap" -- as `json`,
while IUC's `sourmash_compare` requires `sourmash.sig`. The collection could not be piped into the
tool it was made for, and nothing said so: every workflow here was green, because no sourmash tool
is installed on usegalaxy.org at all.

⚠ AND THAT DEFECT WAS NOT ON AN EDGE THIS SCRIPT CAN SEE, which is worth stating up front so the
green result is not over-read. It was an edge to a tool OUTSIDE our workflows. What this catches is
the same mistake made INTERNALLY -- a producer whose declared type cannot satisfy a consumer we
ourselves wired it to. For the external case the question is different and needs a human: does a
more specific datatype exist for what this output actually IS, and does the intended reuser require
it? `--specificity` prints the material for that question; it cannot answer it.

⛔ A PARENT-TYPE DECLARATION IS NOT A SAFE DEFAULT. `data` and `txt` look maximally permissive and
are the opposite: they satisfy only consumers that ask for as little, and foreclose every consumer
that asks for the real thing. Declaring the subtype costs nothing -- it still satisfies the
permissive ones.

    python3 scripts/check_datatype_satisfaction.py             # needs a Galaxy for toolshed steps
    python3 scripts/check_datatype_satisfaction.py --offline   # local definitions only
    python3 scripts/check_datatype_satisfaction.py --specificity
    python3 scripts/check_datatype_satisfaction.py --self-test
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import sys
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET

import yaml

ROOT = pathlib.Path(__file__).resolve().parent.parent

#: Producers whose DECLARED output format is a placeholder rather than a promise, so comparing it
#: against a consumer is meaningless. Each entry names WHY, because the alternative is an
#: allowlist that grows whenever the check is inconvenient.
#:
#: ⚠ THESE ARE NOT EXEMPTIONS, THEY ARE "NOT JUDGEABLE FROM A DECLARATION". A tool with
#: `format_source="input"` emits whatever its input was; `format="auto"` is resolved at runtime.
#: The API reports both as a concrete-looking string (`data`, `auto`), which is exactly what makes
#: a naive comparison produce false alarms -- seven of them here, on edges the 19-genome softmask
#: run proved work by running 421/421 jobs through them.
DYNAMIC_FORMAT = {
    "cat1": "format_source=\"input1\" -- Galaxy's own Concatenate inherits its input's datatype",
    "bedtools_sortbed": "format_source=\"input\" in the IUC wrapper",
    "collapse_dataset": "format_source -- nml/collapse_collections inherits the element datatype",
    "batched_lastz": "format=\"auto\" -- resolved from the chosen output type at runtime",
}


def registry(url: str, key: str) -> tuple[dict, dict]:
    req = urllib.request.Request(f"{url.rstrip('/')}/api/datatypes/mapping",
                                 headers={"x-api-key": key})
    m = json.load(urllib.request.urlopen(req, timeout=90))
    return m["ext_to_class_name"], m["class_to_classes"]


def satisfies(prod_ext: str, want_exts: list[str], e2c: dict, c2c: dict) -> bool | None:
    """Galaxy's own rule. `None` means the server does not know one of the extensions."""
    pc = e2c.get(prod_ext)
    if pc is None:
        return None
    judged = False
    for w in want_exts:
        wc = e2c.get(w)
        if wc is None:
            continue
        judged = True
        if pc == wc or c2c.get(pc, {}).get(wc):
            return True
    return False if judged else None


def merge_exts(into: dict, name: str, exts: list[str]) -> None:
    """Accumulate the UNION of formats a param accepts across conditional cases.

    ⛔ NOT ASSIGNMENT, AND THE DIFFERENCE PRODUCED TEN FALSE ALARMS. A conditional's cases reuse
    one param NAME -- `bedtools_genomecoveragebed` has an `input_type` conditional whose BED case
    and BAM case both declare `input` -- so assigning let the last case seen clobber the first, and
    every `bed` producer feeding it was reported as unsatisfiable. The 19-genome softmask run had
    already put 421/421 jobs through exactly those edges.

    ⚠ THE UNION IS DELIBERATELY OPTIMISTIC. Which case a workflow selects decides which formats
    really apply, and that is in the step's state rather than the tool's interface; accepting a
    producer that satisfies ANY case can therefore miss a genuine mismatch under the case actually
    chosen. Over-reporting is worse here: a checker whose failures are usually spurious is a
    checker that gets ignored, and this one exists to be believed when it does fire.
    """
    if not exts:
        return
    have = into.setdefault(name, [])
    for e in exts:
        if e not in have:
            have.append(e)


def base_id(tool_id: str) -> str:
    """`.../repos/owner/repo/name/version` -> `name`; a bare id is returned unchanged."""
    return tool_id.split("/")[-2] if "/" in tool_id else tool_id


def local_defs() -> dict[str, dict]:
    """Declared formats from this repo's own UDTs and tool XMLs, keyed by tool id AND base id."""
    defs: dict[str, dict] = {}

    def put(tid: str, rec: dict) -> None:
        defs.setdefault(tid, rec)
        defs.setdefault(base_id(tid), rec)

    for p in sorted((ROOT / "udt").glob("*.gxtool.yml")):
        d = yaml.safe_load(p.read_text(encoding="utf-8"))
        rec: dict = {"in": {}, "out": {}, "src": p.name}
        for i in d.get("inputs") or []:
            if i.get("format"):
                rec["in"][i["name"]] = [x.strip() for x in str(i["format"]).split(",")]
        for o in d.get("outputs") or []:
            if o.get("format"):
                rec["out"][o["name"]] = [x.strip() for x in str(o["format"]).split(",")]
            for dd in o.get("discover_datasets") or []:
                if dd.get("format"):
                    rec["out"][o["name"]] = [str(dd["format"]).strip()]
        if d.get("id"):
            put(d["id"], rec)

    for p in sorted((ROOT / "tools").rglob("*.xml")):
        try:
            root = ET.parse(p).getroot()
        except ET.ParseError as exc:
            # ⛔ NOT SKIPPED. A tool XML that does not parse cannot be loaded by Galaxy either, and
            # a checker that quietly ignores it reports a clean pipeline built on a broken tool.
            sys.exit(f"{p}: does not parse as XML ({exc}). Galaxy cannot load it either.")
        if root.tag != "tool" or not root.get("id"):
            continue
        rec = {"in": {}, "out": {}, "src": str(p.relative_to(ROOT))}
        for el in root.iter():
            if el.tag == "param" and el.get("type") in ("data", "data_collection") \
                    and el.get("format"):
                merge_exts(rec["in"], el.get("name"),
                           [x.strip() for x in el.get("format").split(",")])
            elif el.tag in ("data", "collection"):
                # ⚠ format_source/metadata_source means the datatype is INHERITED, so there is no
                # declaration to compare. Recorded as dynamic rather than omitted, so the edge is
                # reported as unjudged instead of silently vanishing.
                if el.get("format_source") or el.get("metadata_source"):
                    rec["out"][el.get("name")] = None
                elif el.get("format"):
                    rec["out"][el.get("name")] = [el.get("format").strip()]
            elif el.tag == "discover_datasets" and el.get("format"):
                rec["out"].setdefault("__discovered__", [el.get("format").strip()])
        put(root.get("id"), rec)
    return defs


def server_defs(url: str, key: str, ids: set[str]) -> dict[str, dict]:
    from bioblend.galaxy import GalaxyInstance

    gi = GalaxyInstance(url.rstrip("/"), key=key)
    installed = {t["id"] for t in gi.tools.get_tools()}

    def walk(ps, pre=""):
        for p in ps:
            if p.get("type") == "conditional":
                for c in p.get("cases", []):
                    yield from walk(c.get("inputs", []), f"{pre}{p['name']}|")
            elif p.get("type") in ("section", "repeat"):
                yield from walk(p.get("inputs", []), f"{pre}{p['name']}|")
            elif p.get("type") in ("data", "data_collection"):
                yield f"{pre}{p['name']}", p.get("extensions") or []

    out: dict[str, dict] = {}
    for tid in sorted(ids):
        resolved = tid if tid in installed else None
        if resolved is None:
            cands = [i for i in installed if f"/{tid}/" in i or base_id(i) == tid]
            resolved = max(cands) if cands else None
        if resolved is None:
            continue
        try:
            d = gi.tools.show_tool(resolved, io_details=True)
        except Exception:                                                    # noqa: BLE001
            continue
        ins: dict[str, list[str]] = {}
        for n, e in walk(d.get("inputs", [])):
            merge_exts(ins, n, e)
        rec = {"in": ins,
               "out": {o["name"]: ([o["format"]] if o.get("format") else None)
                       for o in d.get("outputs", [])},
               "src": f"server:{resolved}"}
        out[tid] = rec
        out.setdefault(base_id(tid), rec)
    return out


def edges() -> list[tuple[str, str, str, str, str, str]]:
    """(workflow, producer_step, producer_tool, output, consumer_step.input, consumer_tool)."""
    found = []
    for wf in sorted((ROOT / "workflows").glob("*/*.gxwf.yml")):
        d = yaml.safe_load(wf.read_text(encoding="utf-8"))
        steps = d.get("steps") or {}
        tid_of = {label: st.get("tool_id") for label, st in steps.items()}
        for label, st in steps.items():
            for iname, src in (st.get("in") or {}).items():
                if not isinstance(src, str) or "/" not in src:
                    continue          # a workflow-level input, not a step output
                pstep, _, oname = src.rpartition("/")
                if pstep not in tid_of:
                    continue
                found.append((wf.name, pstep, tid_of[pstep] or "", oname,
                              f"{label}.{iname}", tid_of[label] or ""))
    return found


def run(e2c: dict, c2c: dict, defs: dict) -> tuple[list[str], list[str], int]:
    problems, notes, judged = [], [], 0
    dynamic_seen: dict[str, int] = {}
    unjudged = 0
    for wfn, pstep, ptid, oname, cons, ctid in edges():
        pdef, cdef = defs.get(ptid), defs.get(ctid)
        if not (pdef and cdef):
            unjudged += 1
            continue
        pexts = pdef["out"].get(oname, "missing")
        if pexts == "missing":
            pexts = pdef["out"].get("__discovered__")
        cexts = cdef["in"].get(cons.split(".", 1)[1])
        if pexts is None or base_id(ptid) in DYNAMIC_FORMAT:
            dynamic_seen[base_id(ptid)] = dynamic_seen.get(base_id(ptid), 0) + 1
            continue
        if not (pexts and cexts):
            unjudged += 1
            continue
        verdict = [satisfies(pe, cexts, e2c, c2c) for pe in pexts]
        if any(v is True for v in verdict):
            judged += 1
            continue
        if all(v is None for v in verdict):
            unjudged += 1
            continue
        judged += 1
        problems.append(f"{wfn}: {pstep}/{oname} emits {pexts}, but {cons} accepts {cexts} -- "
                        f"no declared producer type is that class or a subclass of it, and "
                        f"Galaxy will refuse the connection")
    for tid, n in sorted(dynamic_seen.items()):
        why = DYNAMIC_FORMAT.get(tid, "output declares format_source/metadata_source")
        notes.append(f"{n} edge(s) from `{tid}` not judged: {why}")
    if unjudged:
        notes.append(f"{unjudged} edge(s) had no declared format on one end, or a tool this run "
                     f"could not resolve")
    return problems, notes, judged


def specificity(e2c: dict, c2c: dict, defs: dict) -> None:
    """List our outputs whose declared type has subtypes -- material for a human question."""
    print("Outputs declaring a type that HAS subtypes. Only a human can say whether the data is\n"
          "actually one of them; a parent-type declaration is not wrong, only weaker.\n")
    seen = set()
    for _tid, rec in sorted(defs.items()):
        if not str(rec.get("src", "")).startswith(("udt/", "tools/")):
            continue
        for oname, exts in sorted(rec["out"].items()):
            for ext in exts or []:
                cls = e2c.get(ext)
                if cls is None or (ext, oname) in seen:
                    continue
                seen.add((ext, oname))
                subs = sorted(e for e, c in e2c.items() if e != ext and c2c.get(c, {}).get(cls))
                if subs:
                    print(f"  {rec['src']:44} {oname:20} {ext:12} "
                          f"{len(subs)} subtype(s)")


# --------------------------------------------------------------------------------------------

#: (label, producer exts, consumer exts, must_fail). Uses only the relations asserted below, so
#: the self-test needs no server.
_E2C = {"json": "J", "sourmash.sig": "S", "tabular": "T", "len": "L", "data": "D", "bed": "B"}
_C2C = {"S": {"J": True, "D": True}, "J": {"D": True}, "L": {"T": True, "D": True},
        "T": {"D": True}, "B": {"T": True, "D": True}, "D": {}}

_CASES = [
    ("subtype into supertype input", ["sourmash.sig"], ["json"], False),
    ("supertype into subtype input", ["json"], ["sourmash.sig"], True),
    ("exact match", ["json"], ["json"], False),
    ("len into a tabular input", ["len"], ["tabular"], False),
    ("tabular into a len input", ["tabular"], ["len"], True),
    ("root type into a specific input", ["data"], ["bed"], True),
    ("specific into a root input", ["bed"], ["data"], False),
    ("unrelated siblings", ["bed"], ["sourmash.sig"], True),
    ("one of several producer types matches", ["data", "bed"], ["tabular"], False),
    ("consumer lists several, one matches", ["len"], ["bed", "tabular"], False),
    ("consumer lists several, none match", ["json"], ["bed", "tabular"], True),
]


def self_test() -> int:
    bad = 0
    for label, pexts, cexts, must_fail in _CASES:
        verdict = [satisfies(pe, cexts, _E2C, _C2C) for pe in pexts]
        refused = not any(v is True for v in verdict)
        ok = refused == must_fail
        bad += not ok
        print(f"  {'pass' if ok else '⛔ FAIL'}  {label:40} "
              f"{'refused' if refused else 'accepted'}")

    # an extension the server has never heard of must be UNJUDGED, never a silent pass
    v = satisfies("no.such.ext", ["json"], _E2C, _C2C)
    ok = v is None
    bad += not ok
    print(f"  {'pass' if ok else '⛔ FAIL'}  {'unknown producer extension -> unjudged':40} {v!r}")
    v = satisfies("json", ["no.such.ext"], _E2C, _C2C)
    ok = v is None
    bad += not ok
    print(f"  {'pass' if ok else '⛔ FAIL'}  {'unknown consumer extension -> unjudged':40} {v!r}")

    # ⛔ THE CONDITIONAL-CASE MERGE, because getting this wrong invented ten failures
    got: dict = {}
    merge_exts(got, "input", ["bed", "bedgraph"])
    merge_exts(got, "input", ["bam"])
    merge_exts(got, "input", ["bed"])            # already present, must not duplicate
    ok = got["input"] == ["bed", "bedgraph", "bam"]
    bad += not ok
    print(f"  {'pass' if ok else '⛔ FAIL'}  {'conditional cases union, no duplicates':40} "
          f"{got['input']}")
    merge_exts(got, "empty", [])
    ok = "empty" not in got
    bad += not ok
    print(f"  {'pass' if ok else '⛔ FAIL'}  {'a param with no formats is not recorded':40} "
          f"{'absent' if ok else got.get('empty')}")

    # every DYNAMIC_FORMAT entry must carry a reason, or it is an unexplained allowlist
    unexplained = sorted(k for k, v in DYNAMIC_FORMAT.items() if not (v or "").strip())
    ok = not unexplained
    bad += not ok
    print(f"  {'pass' if ok else '⛔ FAIL'}  {'every dynamic-format entry says why':40} "
          f"{unexplained or 'all ' + str(len(DYNAMIC_FORMAT)) + ' explained'}")

    n = len(_CASES) + 5
    print(f"\n  {n - bad}/{n} self-test case(s) behaved correctly")
    return 1 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--offline", action="store_true",
                    help="judge only edges between this repo's own tools; toolshed consumers are "
                         "reported as unjudged rather than guessed")
    ap.add_argument("--specificity", action="store_true",
                    help="list outputs whose declared type has subtypes, for review")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()

    url, key = os.environ.get("GALAXY_URL", ""), os.environ.get("GALAXY_API_KEY", "")
    if not (url and key):
        # ⛔ THE REGISTRY IS NOT OPTIONAL AND IS NOT GUESSABLE. Subclass relations live on the
        # server; hardcoding a table here would be a second copy that silently disagrees with the
        # Galaxy actually being targeted, which is the whole class of bug this checks for.
        sys.exit("GALAXY_URL and GALAXY_API_KEY must be set: the subclass relations come from "
                 "/api/datatypes/mapping and cannot be inferred locally.")
    try:
        e2c, c2c = registry(url, key)
    except (urllib.error.URLError, OSError) as exc:
        sys.exit(f"could not read the datatype registry from {url}: {exc}")

    defs = local_defs()
    if not a.offline:
        wanted = {t for _, _, p, _, _, c in edges() for t in (p, c)
                  if t and t not in defs}
        try:
            defs.update(server_defs(url, key, wanted))
        except ImportError:
            print("  ⚠ bioblend is not installed; toolshed consumers stay unjudged")

    if a.specificity:
        specificity(e2c, c2c, defs)
        return 0

    problems, notes, judged = run(e2c, c2c, defs)
    for n in notes:
        print(f"  note: {n}")
    # ⛔ 0 == 0 IS NOT A PASS. With no judgeable edge there is nothing to be clean about, and this
    # would report success on a repo whose every format declaration had been deleted.
    if judged == 0:
        print("\n⛔ no edge had a declared, resolvable format on BOTH ends, so nothing was "
              "actually checked. Refusing to report a clean pipeline.")
        return 1
    if problems:
        print(f"\n⛔ {len(problems)} unsatisfiable edge(s) of {judged} judged:")
        for p in problems:
            print(f"  - {p}")
        return 1
    print(f"\n✅ {judged} edge(s) judged: every producer's declared datatype is the consumer's "
          f"class or a subclass of it")
    return 0


if __name__ == "__main__":
    sys.exit(main())
