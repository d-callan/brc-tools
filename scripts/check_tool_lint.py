#!/usr/bin/env python3
"""Run galaxy-tool-refactor's strict ruleset over tools/, against a per-(file, code) baseline.

    check_tool_lint.py [PATH ...] [--baseline FILE] [--write-baseline]

⚠ NEEDS `pip install galaxy-tool-refactor-registry==0.3.10` -- unlike check_tool_defects.py,
which is stdlib only. It must be the REGISTRY distribution, not galaxy-tool-lint: the lint package
does not depend on the registry, so installing it alone leaves this import failing. The version is
PINNED, because a new rule in a later release would otherwise turn CI red on a day nobody touched
this repo, and a lint gate that goes red on its own gets switched off. Raising the pin is a
maintainer step that comes with re-baselining, since a release can move counts in either
direction.

NOTHING IS EXCLUDED ANY MORE. GTR025 (<requirements>) and GTR038 (<citations>) used to be, because
they read the unexpanded tree: every wrapper declaring those through `<expand macro=.../>` was
reported as declaring nothing -- 40 of 42 GTR025 findings here and 44 of 50 GTR038 were that false
positive, which would have buried the real ones. Fixed upstream in 0.3.8, so both are back and the
8 genuine findings between them are visible.

The baseline is keyed on (file, code) and a COUNT, not on a line number or a message: line numbers
move whenever anything above them is edited, and messages get reworded between releases. Either
would make the baseline rot into a file nobody trusts.
"""

from __future__ import annotations

import argparse
import collections
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

# Empty since 0.3.8 fixed the macro blindness. Kept as the seam: a rule that turns out to be
# systematically wrong on this repo gets parked here WITH ITS MEASUREMENT, never silently dropped.
EXCLUDED: frozenset[str] = frozenset()
BASELINE_DEFAULT = Path("tools/.tool-lint-baseline.tsv")


def strict_codes() -> frozenset[str]:
    from galaxy_tool_refactor_registry import facade

    strict = next(r for r in facade.list_rulesets() if r.name == "strict")
    return frozenset(strict.codes) - EXCLUDED


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
        candidates = [path] if path.is_file() else sorted(path.rglob("*.xml"))
        for found in candidates:
            if ".git" in found.parts or "test-data" in found.parts:
                continue
            try:
                # Only a document whose ROOT is <tool> is a wrapper. Matching the text "<tool "
                # instead sweeps in galaxy_config_*.xml, whose <tool> elements are toolbox
                # entries -- that mistake added 112 bogus findings the first time this was run.
                if ET.parse(found).getroot().tag == "tool":
                    yield found
            except ET.ParseError:
                continue


def collect(paths: list[Path]) -> tuple[collections.Counter, dict]:
    from galaxy_tool_refactor_registry import facade

    codes = strict_codes()
    counts: collections.Counter = collections.Counter()
    samples: dict[tuple[str, str], str] = {}
    for tool in iter_tools(paths):
        for violation in facade.detect(tool, codes=codes).violations:
            key = (rel(tool), violation.code)
            counts[key] += 1
            samples.setdefault(key, getattr(violation, "message", ""))
    return counts, samples


def load_baseline(path: Path) -> collections.Counter:
    counts: collections.Counter = collections.Counter()
    if not path.exists():
        return counts
    for line in path.read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            file, code, count = line.split("\t")
            counts[(file, code)] = int(count)
    return counts


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("paths", nargs="*", type=Path, default=[Path("tools")])
    ap.add_argument("--baseline", type=Path, default=BASELINE_DEFAULT)
    ap.add_argument("--write-baseline", action="store_true")
    args = ap.parse_args()

    try:
        counts, samples = collect(args.paths)
    except ImportError:
        print("check_tool_lint: needs `pip install galaxy-tool-refactor-registry==0.3.10`",
              file=sys.stderr)
        return 2

    if args.write_baseline:
        header = (
            "# galaxy-tool-refactor 0.3.10 strict findings per (file, code).\n"
            "# CI fails when a count RISES or a new pair appears.\n"
            "# Lower a number as you fix; never raise one.\n"
        )
        body = "".join(f"{f}\t{c}\t{n}\n" for (f, c), n in sorted(counts.items()))
        args.baseline.write_text(header + body)
        print(f"wrote {len(counts)} (file, code) pair(s), "
              f"{sum(counts.values())} finding(s) to {args.baseline}")
        return 0

    known = load_baseline(args.baseline)
    regressions = sorted(k for k in counts if counts[k] > known.get(k, 0))
    # A count that FELL is good news that needs the baseline lowered, or the ratchet slips back.
    improvements = sorted(k for k in set(known) | set(counts) if counts.get(k, 0) < known.get(k, 0))

    for file, code in regressions:
        was, now = known.get((file, code), 0), counts[(file, code)]
        print(f"WORSE    {file}\n         {code} {was} -> {now}: {samples[(file, code)][:84]}")
    for file, code in improvements:
        print(f"BETTER   {file}: {code} {known[(file, code)]} -> {counts.get((file, code), 0)}"
              f" -- lower it in {args.baseline}")

    total = sum(counts.values())
    print(f"\n{len(regressions)} worse, {len(improvements)} better, {total} finding(s) total")
    if not regressions and not improvements:
        print(f"OK -- at the baseline ({total} finding(s) left to fix)")
    return 1 if (regressions or improvements) else 0


if __name__ == "__main__":
    sys.exit(main())
