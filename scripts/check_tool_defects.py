#!/usr/bin/env python3
"""Refuse the Galaxy tool-XML defects that neither planemo nor galaxy-tool-refactor detects.

    check_tool_defects.py [PATH ...] [--baseline FILE] [--write-baseline] [--self-test]

STANDARD LIBRARY ONLY, like check_udt_definitions.py and for the same reason: a defect landing in
tools/ must be noticed by CI even when nobody's planemo environment is present.

⚠ THIS RUNS AGAINST A BASELINE, because 26 findings were already here when it was written. CI
fails on a NEW finding and ALSO on a baseline line whose defect is gone -- so the file can only
shrink. Delete lines from it as you fix; never regenerate it to make a failure go away.

Every rule here was measured on a real wrapper in this project's orbit, and every one of them
produces a GREEN job that computes the wrong thing. Rules planemo or galaxy-tool-refactor already
cover are deliberately absent -- run those too, this is the residue.

⛔ BOOL-TRUTHY is the whole reason this file exists. A Galaxy boolean renders into Cheetah as the
*string* `truevalue` or `falsevalue`, never a Python bool. `falsevalue="no"` is a non-empty string,
so `#if $flag` is true whether the box is ticked or not. The flag becomes unconditional and the
job still exits 0. `falsevalue=""` is the only spelling for which a bare `#if` happens to work,
which is exactly why the bug survives review: most wrappers use it, so the pattern looks safe.

⚠ BOOL-TEST-ONE-SIDED catches it from the other side. Four live wrappers carried BOOL-TRUTHY and
every single <test> set the boolean to "true" -- the false branch, the only branch the bug breaks,
was never exercised. A boolean tested on one side only is an untested boolean.
"""

from __future__ import annotations

import argparse
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

# A bare Cheetah truth test: `#if $x`, `#if not $x`, `#if $a.x` -- with no comparison following.
_BARE_IF = re.compile(
    r"#if\s+(?:not\s+)?\$\{?((?:\w+\.)*)(\w+)\}?\s*(?:##.*)?$",
    re.MULTILINE,
)
_IS_TEST = re.compile(r"#if\s+\$\{?((?:\w+\.)*)(\w+)\}?\s+is\s+(?:True|False)\b")


def _text(node: ET.Element | None) -> str:
    return "".join(node.itertext()) if node is not None else ""


def _param_name(param: ET.Element) -> str | None:
    """The name Galaxy derives, from `name` or -- as Galaxy itself does -- from `argument`."""
    name = param.get("name")
    if name:
        return name
    argument = param.get("argument")
    if argument:
        return argument.lstrip("-").replace("-", "_")
    return None


def _bare_refs(command: str) -> set[str]:
    return {m.group(2) for m in _BARE_IF.finditer(command)} | {
        m.group(2) for m in _IS_TEST.finditer(command)
    }


def check_tool(path: Path) -> list[tuple[str, str, str]]:
    """Return (code, subject, message) per defect.

    ``subject`` is the param name, never the prose: the baseline is keyed on it so
    that rewording a message cannot silently invalidate every baselined line.
    """
    raw = path.read_text(errors="replace")
    out: list[tuple[str, str, str]] = []

    # XML-DASH first: a `--` inside a comment makes the document unparseable, so it has to be
    # found textually, before ET is allowed to fail on it.
    for match in re.finditer(r"<!--(.*?)(?:-->|$)", raw, re.DOTALL):
        if "--" in match.group(1):
            out.append(
                ("XML-DASH", "comment", "`--` inside an XML comment is fatal: "
                 f"{match.group(1)[:40].strip()!r}")
            )

    try:
        root = ET.fromstring(raw)
    except ET.ParseError as exc:
        out.append(("XML-PARSE", "document", str(exc)))
        return out
    if root.tag != "tool":
        return []

    command = _text(root.find("command"))
    bare = _bare_refs(command)
    filters = {
        f.text.strip()
        for out_el in (root.find("outputs") if root.find("outputs") is not None else [])
        for f in out_el.findall("filter")
        if f.text
    }

    tested: dict[str, set[str]] = {}
    for test in root.iter("test"):
        for param in test.findall("param"):
            if param.get("name") and param.get("value") is not None:
                tested.setdefault(param.get("name"), set()).add(param.get("value").lower())

    for param in root.iter("param"):
        name = _param_name(param)
        if not name:
            continue
        ptype = param.get("type", "")

        if ptype == "boolean":
            # ⛔ AN ABSENT falsevalue IS NOT AN EMPTY ONE. Galaxy defaults the pair to the
            # strings 'true'/'false' (galaxy.tool_util.parser.util.boolean_true_and_false_values,
            # checked against profiles 21.09 / 24.0 / 26.1), so a param declaring neither is the
            # WORST case: `#if $flag` is always true and there is no wrong-looking attribute to
            # notice. Reading `falsevalue` alone missed six live instances in another repo.
            declared = param.get("falsevalue")
            falsevalue = "false" if declared is None else declared
            if name in bare and falsevalue:
                why = (
                    f"falsevalue={declared!r} is a non-empty string, not a Python false"
                    if declared is not None
                    else "it declares no truevalue/falsevalue, so Galaxy defaults them to the "
                         "strings 'true'/'false'"
                )
                out.append(("BOOL-TRUTHY", name, f"`#if ${name}` is ALWAYS true: {why}"))
                if name in filters:
                    out.append(
                        ("BOOL-FILTER", name, f"`{name}` also gates an <output> <filter>; the command "
                                        f"emits unconditionally while the dataset may not exist")
                    )
            seen = tested.get(name, set())
            # ⛔ `if seen and ...` WAS THE BLIND SPOT. A boolean no <test> ever SETS yields an
            # empty `seen`, so the one-sided rule stayed silent and the ratchet read clean --
            # which is strictly worse than one-sided, because NEITHER branch is exercised.
            # Measured when this was added: 9 booleans across 4 hub wrappers were invisible.
            if not seen:
                out.append(
                    ("BOOL-TEST-NEVER", name,
                     f"boolean `{name}` is never set by any <test>, so neither branch is "
                     f"exercised -- add a case that sets it and assert the flag's presence "
                     f"with <assert_command>, and one that leaves it default")
                )
            elif not ({"true", "yes", "1"} & seen and {"false", "no", "0"} & seen):
                out.append(
                    ("BOOL-TEST-ONE-SIDED", name,
                     f"boolean `{name}` is only ever tested as {sorted(seen)} -- the other branch "
                     f"is never exercised")
                )

        elif ptype in ("integer", "float") and name in bare:
            out.append(
                ("NUM-ZERO", name, f"`#if ${name}` silently drops the value 0, which is a legal "
                             f"{ptype}; use `#if str(${name}):`")
            )

        # galaxy-tool-refactor's GTR034 reads param.get("name") and skips the param when it is
        # absent, so an argument-only param is never checked for being orphaned. Do it properly.
        # The param may be reached through its enclosing conditional/section --
        # `$ref_mode.path`, not just `$path`. Requiring the bare spelling produced two false
        # positives on a live wrapper.
        if (
            not param.get("name")
            and param.get("argument")
            and not re.search(rf"\${{?(?:\w+\.)*{re.escape(name)}\b", command)
        ):
            out.append(
                ("ARG-ORPHAN", name, f"`{name}` (from argument={param.get('argument')!r}) is never "
                               f"referenced in <command> -- GTR034 cannot see this one")
            )
    return out


def rel(path: Path) -> str:
    """Repo-relative POSIX path.

    ⚠ THE BASELINE IS KEYED ON THIS. Recording ``path.as_posix()`` embeds the checkout directory,
    so a baseline written locally matches nothing under /home/runner/work and CI reports every
    finding as new. Normalise before the key is built, never after.
    """
    try:
        return path.resolve().relative_to(Path.cwd().resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def iter_tools(paths: list[Path]):
    for path in paths:
        if path.is_file():
            yield path
        else:
            for found in sorted(path.rglob("*.xml")):
                if ".git" in found.parts or "test-data" in found.parts:
                    continue
                yield found


def _self_test() -> int:
    """Each fixture must trip exactly the rule it is named for, and the clean one must trip none."""
    cases = [
        ("BOOL-TRUTHY", """<tool id="t" name="T" version="1"><command><![CDATA[
            p #if $flag
            --go
            #end if]]></command><inputs>
            <param name="flag" type="boolean" truevalue="yes" falsevalue="no"/>
            </inputs><outputs/></tool>"""),
        ("BOOL-TRUTHY", """<tool id="t" name="T" version="1"><command><![CDATA[
            p #if $flag
            --go
            #end if]]></command><inputs>
            <param argument="--flag" type="boolean" checked="false"/>
            </inputs><outputs/></tool>"""),
        ("NUM-ZERO", """<tool id="t" name="T" version="1"><command><![CDATA[
            p #if $n
            --n $n
            #end if]]></command><inputs>
            <param name="n" type="integer" value="1"/></inputs><outputs/></tool>"""),
        ("ARG-ORPHAN", """<tool id="t" name="T" version="1"><command><![CDATA[p --x]]></command>
            <inputs><param argument="--min-score" type="float" value="1"/></inputs>
            <outputs/></tool>"""),
        ("XML-DASH", """<tool id="t" name="T" version="1"><!-- a -- b --><command>p</command>
            <inputs/><outputs/></tool>"""),
        # A param reached through its enclosing conditional (`$mode.path`) is NOT an orphan.
        # This spelling produced two false positives on a live wrapper before it was pinned.
        ("", """<tool id="t" name="T" version="1"><command><![CDATA[
            p #for $x in str($mode.path).split(',')
            --path '$x'
            #end for]]></command><inputs>
            <conditional name="mode"><param name="sel" type="select"><option value="a"/></param>
            <when value="a"><param argument="--path" type="text" value=""/></when>
            </conditional></inputs><outputs/></tool>"""),
        # ⛔ Found by the coverage check below the moment it became a real check: these two
        # rules had no fixture either, so the file had been claiming to cover three rules it
        # did not exercise at all.
        ("XML-PARSE", "<tool><command>unclosed"),
        # the boolean is bare in the command AND gates an <output> <filter>; BOOL-TRUTHY fires
        # alongside, which is fine -- the driver asserts the expected code is PRESENT
        ("BOOL-FILTER", """<tool id="t" name="T" version="1">
            <command><![CDATA[p #if $flag
            --go
            #end if]]></command><inputs>
            <param name="flag" type="boolean" truevalue="--go" falsevalue="no"/></inputs>
            <outputs><data name="o" format="txt"><filter>flag</filter></data></outputs></tool>"""),
        # ⛔ NEITHER branch exercised, which the one-sided rule could not see: it keyed on the
        # booleans a test SETS, so an untested boolean produced an empty set and stayed silent.
        ("BOOL-TEST-NEVER", """<tool id="t" name="T" version="1">
            <command><![CDATA[p $flag]]></command><inputs>
            <param name="flag" type="boolean" truevalue="--go" falsevalue=""/></inputs>
            <outputs/><tests><test><param name="other" value="1"/></test></tests></tool>"""),
        ("BOOL-TEST-ONE-SIDED", """<tool id="t" name="T" version="1">
            <command><![CDATA[p $flag]]></command><inputs>
            <param name="flag" type="boolean" truevalue="--go" falsevalue=""/></inputs>
            <outputs/><tests><test><param name="flag" value="true"/></test></tests></tool>"""),
        # falsevalue="" is the one spelling for which a bare #if is correct -- it must stay silent,
        # or the rule would fire on most of the ecosystem and be turned off.
        ("", """<tool id="t" name="T" version="1"><command><![CDATA[
            p #if $flag
            --go
            #end if]]></command><inputs>
            <param name="flag" type="boolean" truevalue="--go" falsevalue=""/></inputs>
            <outputs/><tests><test><param name="flag" value="true"/></test>
            <test><param name="flag" value="false"/></test></tests></tool>"""),
    ]
    failures = 0
    for index, (expected, xml) in enumerate(cases):
        tmp = Path(f"/tmp/_gtcheck_{index}_{expected or 'clean'}.xml")
        tmp.write_text(xml)
        codes = {code for code, _subject, _message in check_tool(tmp)}
        tmp.unlink()
        if expected and expected not in codes:
            print(f"  FAIL {expected}: not raised (got {sorted(codes) or 'nothing'})")
            failures += 1
        elif not expected and codes:
            print(f"  FAIL clean fixture raised {sorted(codes)}")
            failures += 1
        else:
            print(f"  ok   {expected or 'clean fixture stays silent'}")
    # ⛔ THIS LINE USED TO BE A PRINT, NOT A CHECK. It said "self-test covers every rule this
    # file defines" whenever `failures == 0`, which is a different statement entirely -- a rule
    # added with no fixture printed that sentence and was never exercised. BOOL-TEST-NEVER was
    # added in exactly that state. The codes are scraped from this file's own rule emissions, so
    # a new rule with no case here now FAILS instead of being congratulated.
    emitted = set(re.findall(r'out\.append\(\s*\(\s*"([A-Z0-9-]+)"',
                             Path(__file__).read_text()))
    exercised = {expected for expected, _xml in cases if expected}
    uncovered = sorted(emitted - exercised)
    if uncovered:
        print(f"  FAIL rules with no self-test fixture: {uncovered}")
        failures += len(uncovered)
    else:
        print(f"  ok   every one of the {len(emitted)} rule(s) this file defines has a fixture")
    print("  ok   self-test passed" if failures == 0 else f"  {failures} self-test failure(s)")
    return 1 if failures else 0


BASELINE_DEFAULT = Path("tools/.tool-defects-baseline.tsv")


def load_baseline(path: Path) -> set[tuple[str, str, str]]:
    if not path.exists():
        return set()
    entries = set()
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        file, code, subject = line.split("\t")
        entries.add((file, code, subject))
    return entries


def write_baseline(path: Path, found: list[tuple[str, str, str]]) -> None:
    header = (
        "# Galaxy tool defects present when check_tool_defects.py was introduced.\n"
        "# CI fails on anything NOT listed here, and on any line here whose defect is gone.\n"
        "# DELETE LINES AS YOU FIX. Never regenerate this to clear a failure -- that is exactly\n"
        "# the ratchet it exists to hold.\n"
    )
    body = "".join(f"{f}\t{c}\t{s}\n" for f, c, s in sorted(found))
    path.write_text(header + body)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("paths", nargs="*", type=Path, default=[Path("tools")])
    ap.add_argument("--baseline", type=Path, default=BASELINE_DEFAULT,
                    help="known findings to tolerate (default: %(default)s)")
    ap.add_argument("--write-baseline", action="store_true",
                    help="record the CURRENT findings as the baseline, then exit")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        return _self_test()

    found: list[tuple[str, str, str]] = []
    messages: dict[tuple[str, str, str], str] = {}
    for tool in iter_tools(args.paths):
        for code, subject, message in check_tool(tool):
            key = (rel(tool), code, subject)
            found.append(key)
            messages[key] = message

    if args.write_baseline:
        write_baseline(args.baseline, found)
        print(f"wrote {len(found)} finding(s) to {args.baseline}")
        return 0

    known = load_baseline(args.baseline)
    current = set(found)
    new = sorted(current - known)
    # A baselined defect that no longer fires is a line someone has to delete. Failing on it is
    # what makes this a ratchet rather than a permanent amnesty.
    fixed = sorted(known - current)

    for key in new:
        file, code, subject = key
        print(f"NEW      {file}\n         {code:20s} {messages[key]}")
    for file, code, subject in fixed:
        print(f"FIXED    {file}: {code} on {subject!r} is gone -- "
              f"delete this line from {args.baseline}")

    carried = len(current & known)
    print(f"\n{len(new)} new, {len(fixed)} fixed-but-still-baselined, {carried} carried")
    if not new and not fixed:
        print(f"OK -- at the baseline ({carried} known finding(s) left to fix)")
    return 1 if (new or fixed) else 0


if __name__ == "__main__":
    sys.exit(main())
