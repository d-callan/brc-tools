#!/usr/bin/env python3
"""Emit `{a}_{b}\t{a}.{b}` for every ordered pair of DISTINCT collection element identifiers.
The cross-product cells are named `A_B` (underscore join); downstream Phase E expects
`A.B`. This 2-col TSV drives WF-C's __RELABEL_FROM_FILE__ step."""
import argparse
import itertools
import pathlib

ap = argparse.ArgumentParser()
# ⛔ `--ids-file` EXISTS BECAUSE A UDT CANNOT BUILD `--ids`. The UDT template used a shell
# command substitution wrapped around a Galaxy expression to turn the identifier FILE into a
# space-separated string. Galaxy owns that syntax in a shell_command and interpolates it ANYWHERE,
# heredocs included, so the outer substitution is consumed as a Galaxy expression, the command line
# cannot be built at all, and the job fails with exit=None and an empty command before anything
# runs. Measured on BOTH usegalaxy.org (26.1) and a local 25.0: "Error occurred while building
# command line". Those two tools had therefore never run on any Galaxy.
#
# ⚠ scripts/build_softmask_udts.py already refuses that construct -- its FORBIDDEN constant names
# the two-character opening sequence -- and build_inventory_udts.py did not. Passing the PATH and
# reading the file here needs no substitution at all, so nothing can be interpolated.
#
# ⚠ AND THIS COMMENT MUST NOT CONTAIN THE SEQUENCE IT DESCRIBES. This script is INLINED into the
# UDT's heredoc, so a literal example here would be interpolated exactly like real code -- and the
# generator's own guard rejected an earlier draft of this comment for precisely that reason.
ap.add_argument("--ids", help="space-separated element identifiers")
ap.add_argument("--ids-file", help="file of element identifiers, one per line")
ap.add_argument("--out", required=True)
a = ap.parse_args()
# ⚠ `is None`, NOT falsiness. The classic XML wrappers always pass --ids, and an EMPTY collection
# renders `--ids ''` -- which is a valid request for zero identifiers, and used to produce a
# zero-row file. Testing truthiness turned that into "one of --ids or --ids-file is required" and
# exit 2, a regression against the merged tool for the one input most likely to be automated.
if a.ids is None and a.ids_file is None:
    ap.error("one of --ids or --ids-file is required")
# ⛔ BOTH IS AN ERROR, NOT A PRECEDENCE. Silently preferring one meant `--ids "" --ids-file real.txt`
# produced a ZERO-ROW file from a fully populated list -- and a zero-row self-pair list removes
# nothing while a zero-row relabel map renames nothing, both without a word. Refuse the ambiguity.
if a.ids is not None and a.ids_file is not None:
    ap.error("give --ids or --ids-file, not both -- which one wins is not something to guess at")
# ⚠ `is not None` HERE TOO. `--ids ''` is falsy but PRESENT, and testing truthiness sent it to
# the file branch with ids_file=None -- a TypeError instead of the zero-row file the empty-string
# case is supposed to produce. The guard above and this selection must agree on what "given" means.
# ⚠ utf-8-SIG, NOT utf-8. A BOM on a hand-made identifier file survives into the first name --
# `\ufeffcs10` -- which then matches no collection element, silently no-opping exactly one row.
ids = ([x for x in a.ids.split() if x] if a.ids is not None
       else [x.strip() for x in
             pathlib.Path(a.ids_file).read_text(encoding="utf-8-sig").splitlines()
             if x.strip()])
# ⛔ AN EMPTY --ids-file IS THE LAST SILENT ZERO, AND IT IS THE ONLY PATH THE UDTs USE.
# `--ids ''` is a deliberate request for zero and stays one (see above); an empty FILE is not a
# request, it is an empty collection arriving from `collection_element_identifiers`. It produced a
# 0-byte output with exit 0, which the workflow then renamed and tagged `WF-C input: self_pairs` /
# `wfc_relabel_map` -- so WF-C could be handed a no-op diagonal filter and a no-op relabel map out
# of a green WF-A run. The sourmash step refuses an empty panel, but these are parallel branches
# and nothing else in WF-A looks at these two files.
if a.ids_file is not None and not ids:
    raise SystemExit(f"{a.ids_file} holds no identifiers. That is an empty collection, not a "
                     f"request for zero rows -- pass --ids '' if zero is genuinely intended. A "
                     f"0-byte output here reads downstream as 'nothing to exclude', not a failure.")
# ⛔ A TAB IN AN IDENTIFIER BREAKS THE COLUMN CONTRACT DOWNSTREAM. relabel_map emits
# `{a}_{b}<TAB>{a}.{b}`, so a tab inside a name yields a THREE-column row and
# `__RELABEL_FROM_FILE__` reads the wrong field. --ids could never carry one (it splits on
# whitespace); --ids-file can, so the check belongs here.
_bad = [i for i in ids if "\t" in i]
if _bad:
    raise SystemExit(f"identifier(s) contain a tab, which breaks the output's column contract: "
                     f"{_bad[:3]}")
# ⛔ THE DIAGONAL IS SKIPPED SO THE ROW COUNT MATCHES THE COLLECTION, WHICH IS WHAT LETS
# `__RELABEL_FROM_FILE__` RUN IN STRICT MODE. WF-C removes the `A_A` cells with the self_pairs
# list before relabelling, so a map carrying them has n**2 rows against n**2-n elements. Galaxy's
# strict mode checks BOTH the row count and every lookup, and the count check is the one that
# refuses first -- so while this emitted the diagonal, strict mode could not be turned on, and
# without it an element missing from the map keeps its own `A_B` identifier and no error is raised
# anywhere. Phase E then splits ids on `.` and silently cannot key those rows.
# scripts/check_relabel_strict.py fails if this and the workflows' `strict` ever drift apart.
#
# ⚠ ONE identifier therefore yields an EMPTY map, and that is correct rather than a silent zero:
# the cross product of a single element is just its diagonal, which WF-C filters out, so the
# collection reaching the relabel step is empty too and 0 == 0 passes strict mode honestly.
pathlib.Path(a.out).write_text(
    "".join(f"{x}_{y}\t{x}.{y}\n" for x, y in itertools.product(ids, ids) if x != y))
