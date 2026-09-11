#!/usr/bin/env python3
"""Prefix the src names of a pairwise MAF with strain names.

multiz decides which species a row belongs to by reading the text before the
first '.' in the src field. axtToMaf writes whatever the AXT carried -- for our
alignments that is a bare contig accession like ``QMFC01000014.1`` -- so without
this step multiz sees one species per contig and the progressive fold produces
nonsense.

Each block of a pairwise MAF holds the target (the hinge) first and the query
second, so the row position determines which strain name to apply.

The rewrite is idempotent: a src that already carries the right strain prefix is
left alone, which keeps hand-written MAFs (and the tool's own test data) working
unchanged.
"""
import argparse
import re
import sys

# s <src> <start> <size> <strand> <srcSize> <text>; capture the src alone so the
# original column spacing on the rest of the line survives untouched
S_LINE = re.compile(r"^(s\s+)(\S+)(\s.*)$")
NAME_OK = re.compile(r"^[A-Za-z0-9_.-]+$")


def prefixed(src, strain):
    if src == strain or src.startswith(strain + "."):
        return src
    return f"{strain}.{src}"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hinge", required=True, help="strain name for the target (first) row")
    p.add_argument("--query", required=True, help="strain name for the query (second) row")
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    a = p.parse_args()

    for name in (a.hinge, a.query):
        if not NAME_OK.match(name):
            sys.exit(f"strain name must be letters, digits, underscore, hyphen or dot: {name!r}")

    row = 0          # position of the current s-line within its block
    blocks = 0
    rewritten = 0
    with open(a.input) as fh, open(a.output, "w") as out:
        for line in fh:
            if line.startswith("a"):
                row = 0
                blocks += 1
            m = S_LINE.match(line) if line.startswith("s") else None
            if m:
                strain = a.hinge if row == 0 else a.query
                new = prefixed(m.group(2), strain)
                if new != m.group(2):
                    line = f"{m.group(1)}{new}{m.group(3)}\n"
                    rewritten += 1
                row += 1
            out.write(line)

    print(f"{a.input}: {blocks} blocks, {rewritten} src names prefixed "
          f"(target={a.hinge}, query={a.query})", file=sys.stderr)


if __name__ == "__main__":
    main()
