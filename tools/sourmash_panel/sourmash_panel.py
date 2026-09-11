#!/usr/bin/env python3
"""Sketch every assembly in a collection, then compare -- labelling by element identifier.

The identifiers cannot be read from the collection inside a job (see the module docstring of
scripts/build_inventory_udts.py); they are supplied as a file, in collection order.
"""
import argparse
import concurrent.futures
import json
import pathlib
import re
import shutil
import subprocess
import sys

ap = argparse.ArgumentParser()
ap.add_argument("--rendered", required=True, help="file holding the rendered collection input")
ap.add_argument("--ids", required=True, help="element identifiers, one per line, in collection order")
ap.add_argument("--processes", default="1",
                help="how many assemblies to sketch at once, and --processes for each sourmash "
                     "compare. Pass $GALAXY_SLOTS; a UDT gets ONE core unless its definition "
                     "carries a `resource` requirement.")
ap.add_argument("--ksize", default="31")
ap.add_argument("--scaled", default="1000")
ap.add_argument("--containment", default="false",
                help="also emit an asymmetric containment matrix. OFF by default: it adds an "
                     "output, and WF-A's output set is consumed downstream.")
a = ap.parse_args()

# ⛔ THE KEYS ARE QUOTED. Galaxy renders a collection input as JSON -- {"path": "/..."} -- not as a
# JavaScript object literal, so a pattern written for a bare path key matches NOTHING and this
# script reported "0 assemblies" against a perfectly good render: the heredoc held six well-formed
# File records and every one was missed, which is why WF-A could never sketch anything. Parse it as
# JSON, and fall back to the pattern only if that fails, so a future change in render shape
# degrades loudly instead of silently finding zero.
_raw = pathlib.Path(a.rendered).read_text()
_json_ok = True
try:
    _recs = json.loads(_raw[_raw.index("["):_raw.rindex("]") + 1])
    paths = [r["path"] for r in _recs if isinstance(r, dict) and r.get("path")]
except Exception:                                   # noqa: BLE001 -- see below
    # ⚠ THE FALLBACK IS A DIAGNOSTIC, NOT A SECOND PARSER. Once it scraped even ONE path the
    # render fault below stopped firing and the count check took over -- so a render shaped
    # {"class":"File","path":...} instead of a list, which is exactly the "future change in render
    # shape" this fallback anticipates, was reported as "1 assemblies but 2 identifiers. The
    # identifier file must come from the SAME collection" and sent the operator to inspect the one
    # input that was correct. Remember that the JSON parse failed, and say so.
    _json_ok = False
    paths = re.findall(r'"?path"?\s*:\s*"?([^",}\s]+)', _raw)
# ⚠ utf-8-SIG, NOT utf-8, and the same identifier checks the two sibling helpers already make.
# A BOM survives `.strip()` and lands in the first name, so the matrix comes out labelled
# `﻿cs10` -- which then joins against nothing in the self-pair or relabel files built from the
# SAME file. Measured. A tab is rejected there and was accepted here; the CSV header took it.
ids = [x.strip() for x in
       pathlib.Path(a.ids).read_text(encoding="utf-8-sig").splitlines() if x.strip()]
_bad = [i for i in ids if "	" in i]
if _bad:
    sys.exit(f"identifier(s) contain a tab, which lands in the similarity matrix header and breaks "
             f"every column contract downstream: {_bad[:3]}")

# ⛔ THE ASSERTION IS THE POINT. Pairing by position is correct only while the two lists describe
# the same collection; if they do not, every row of the matrix is mislabelled and nothing
# downstream can tell. Refuse instead.
# ⛔ DIAGNOSE THE RENDER BEFORE BLAMING THE IDENTIFIERS. This comparison used to run first, so an
# unparseable or empty rendered block reported "0 assemblies but N identifiers -- the identifier file
# must come from the SAME collection", pointing the operator at the one input that was correct. A
# render this cannot read is a DIFFERENT fault and says so.
# ⛔ AN EMPTY COLLECTION IS NOT A RENDER FAULT, AND SAYING SO SENT THE OPERATOR THE WRONG WAY.
# These two states shared one message: an empty `assemblies` renders `[]`, json.loads SUCCEEDS,
# paths is empty, and the operator was told "could not parse the rendered collection input as
# JSON" -- asserting a parse failure that did not happen, about the one input that was fine. That
# is the same misdirection this block was rewritten to remove; it just moved.
if _json_ok and not paths:
    sys.exit("the collection input rendered cleanly and is EMPTY -- zero assemblies. Nothing to "
             "sketch, and a similarity matrix over nothing is not a result. Check the collection "
             "you passed as `assemblies`; the identifier file is not implicated.")
if not _json_ok:
    sys.exit(f"could not parse the rendered collection input as JSON; the fallback pattern "
             f"scraped {len(paths)} path(s). That is a "
             "RENDER problem, not an identifier problem -- the identifier file is not implicated. "
             "The block should be a JSON array of File records; check what the tool actually "
             "received before changing anything about the identifiers.")
# ⛔ DUPLICATE IDENTIFIERS SILENTLY DESTROY A ROW. Each signature is staged at `stage/{name}.sig`,
# so two elements sharing a name overwrite one file -- and the count assertion below still passes,
# because the COUNTS match. The result is a fully populated matrix in which one genome does not
# appear and another appears twice, which is exactly the "runs and mislabels its output" failure
# this script's assertions exist to prevent.
if len(set(ids)) != len(ids):
    _dupes = sorted({i for i in ids if ids.count(i) > 1})
    sys.exit(f"duplicate element identifier(s) {_dupes[:3]}: a repeat would put two genomes in one "
             f"row of the matrix. Counts alone cannot detect this.")
# ⚠ AND THESE NAMES ARE STILL REFUSED, THOUGH THE REASON HAS CHANGED. While signatures were staged
# at `stage/{name}.sig` this was a data-loss guard: `gA` and `./gA` are distinct strings naming the
# SAME file, so one sketch overwrote the other, the matrix came out `./gA,./gA` with an
# off-diagonal 1.0 that is a self-comparison, and a genome was absent -- exit 0. `../escaped` wrote
# the signature outside the job directory. Staging by INDEX (below) removed both hazards from the
# SKETCH, and an earlier version of this comment concluded that what was left was only a cosmetic
# label problem. That is no longer true in either half:
#   * the name still travels into `--name` and becomes a COLUMN LABEL in similarity.csv, where a
#     leading dash or an embedded path reads as a malformed strain and joins against nothing; and
#   * the name is a FILENAME again -- `signatures/{name}.sig`, published as the discovered
#     `signatures` collection at the end of this script. `gA` and `./gA` would collide there
#     exactly as they once did in `stage/`, and `../escaped` would write outside the work dir.
# So this guard is load-bearing for the data and not just for the labels; the count check beside
# the publish step is the second line of defence.
_bad = [i for i in ids if "/" in i or i in (".", "..") or i.startswith("-")]
if _bad:
    sys.exit(f"element identifier(s) {_bad[:3]} contain a path separator, are a directory alias, or "
             f"begin with a dash. Galaxy allows them; this tool will not stage them, because such a "
             f"name can silently overwrite another element's signature or write outside the job.")
if len(paths) != len(ids):
    sys.exit(f"refusing to guess: {len(paths)} assemblies but {len(ids)} identifiers. "
             f"The identifier file must come from the SAME collection, via the IUC "
             f"collection_element_identifiers tool.")
subprocess.run(["mkdir", "-p", "stage"], check=True)

# ⚠ `sourmash sketch` HAS NO PARALLELISM FLAG, so concurrency here means several invocations at
# once. THREADS, not processes: each sketch is already its own OS process, so the GIL is released
# for the whole of it and a thread pool adds neither pickling nor a second copy of anything.
# ⚠ SIZED FROM $GALAXY_SLOTS, WHICH IS 1 UNLESS THIS TOOL ASKS FOR MORE. Sketching is CPU-bound on
# k-mer hashing and cheap in memory -- at scaled=1000 a 900 Mb assembly is ~900k hashes, ~7 MB --
# so cores are the whole win here, and RAM matters for the compare passes rather than for this.
_workers = max(1, int(a.processes or 1))


def _sketch_one(_i, path, name):
    """Sketch ONE assembly and count its hashes; the ordered pass below does the judging."""
    sig = f"stage/{_i:04d}.sig"
    subprocess.run(["sourmash", "sketch", "dna", "-p", f"k={a.ksize},scaled={a.scaled}",
                    "--name", name, "-o", sig, path], check=True)
    _n = sum(len(s.get("mins") or []) for rec in json.loads(pathlib.Path(sig).read_text())
             for s in (rec.get("signatures") or []))
    return _i, path, name, sig, _n


# ⛔ SUBMITTED AT ONCE, COLLECTED IN PANEL ORDER, AND THAT ORDER IS A CORRECTNESS PROPERTY.
# `as_completed` would surface whichever sketch finished first, so on a panel with two bad members
# the refusal below would name a different genome on every run -- and that refusal is fatal for
# the WHOLE panel, so which genome it names is the operator's only lead. Iterating the futures in
# submission order blocks on 0, then 1, and re-raises the LOWEST-INDEXED failure, which is what
# the serial loop this replaced did.
# ⛔ EVERY FAILURE IS COLLECTED AND REPORTED TOGETHER, NOT JUST THE FIRST. This refusal is fatal
# for the WHOLE panel, and naming one bad genome sends the operator round a loop that costs a full
# re-run per offender: fix one, re-stage, wait hours, meet the next. Sketching has already been
# done for every member by the time anything is judged, so bailing early saves nothing -- the
# information is all in hand, and withholding it is the only cost.
#
# ⚠ FUTURES ARE STILL DRAINED IN SUBMISSION ORDER, so the report comes out in PANEL order rather
# than completion order. `as_completed` would list the same genomes in an order that changes from
# run to run, which makes two runs of one broken panel look like two different problems.
_sketched, _failed = [], []
with concurrent.futures.ThreadPoolExecutor(max_workers=_workers) as _ex:
    _futs = [(_i, name, _ex.submit(_sketch_one, _i, path, name))
             for _i, (path, name) in enumerate(zip(paths, ids, strict=True))]
    for _i, name, _f in _futs:
        try:
            _sketched.append(_f.result())
        except Exception as _exc:                                    # noqa: BLE001
            # ⚠ CAPTURED, NOT SWALLOWED. Each one is reprinted below and the run exits non-zero;
            # the broad catch is here so ONE crashed sketch cannot hide the other 218 results.
            _failed.append((_i, name, f"{type(_exc).__name__}: {_exc}"))

if _failed:
    _lines = "\n".join(f"  [{_i:04d}] {name}: {msg}" for _i, name, msg in _failed)
    sys.exit(f"{len(_failed)} of {len(ids)} assemblies failed to sketch, in panel order:\n"
             f"{_lines}\n"
             f"All of them are listed so the panel needs ONE repair pass, not one per genome.")

# ⚠ THE FLOOR IS A CONSTANT, DEFINED ONCE. It used to be assigned inside the per-genome loop,
# which read like a per-genome property; it is a property of `--scaled`.
_floor = 20
_empty, _thin = [], []
for _i, path, name, _sig, _n in _sketched:
    if _n == 0:
        _empty.append((_i, name, path))
    elif _n < _floor:
        _thin.append((_i, name, _n))
    print(f"sketched {name} <- {path}  ({_n} hashes)", file=sys.stderr)

for _i, name, _n in _thin:
    print(f"⚠ {name} sketched only {_n} hash(es) at scaled={a.scaled} (floor {_floor}). Its row "
          f"will look like an unrelated genome's whatever it actually is -- a 5 kb genome at 7 "
          f"hashes is byte-identical in the matrix to something that shares nothing. Treat its "
          f"similarities as unmeasured, or lower --scaled for the whole panel.", file=sys.stderr)

if _empty:
    _lines = "\n".join(f"  [{_i:04d}] {name} <- {path}" for _i, name, path in _empty)
    sys.exit(f"{len(_empty)} of {len(ids)} assemblies sketched to ZERO hashes at "
             f"scaled={a.scaled}, so each can only appear in the matrix as a genome sharing "
             f"nothing with anyone -- which is not a failure the run would otherwise report:\n"
             f"{_lines}\n"
             f"Lower --scaled, or drop those elements. ⚠ It is N-masking that empties a sketch, "
             f"NOT soft-masking: a 200 kb genome in all-lowercase gives 195 hashes, measured. "
             f"All offenders are listed so this takes ONE repair pass, not one per genome.")

# ⚠ THE PATHS THE WORKERS ACTUALLY WROTE, not a second reconstruction of them. This used to
# rebuild `stage/{i:04d}.sig` from the index, which is the same filename computed in two
# places -- change the naming in `_sketch_one` and the compare below would silently read a
# set of files that no longer exists. `_sketched` is in panel order, which is the order the
# matrix columns must be in.
sigs = [sig for _i, path, name, sig, _n in _sketched]
subprocess.run(["sourmash", "compare", "--processes", a.processes, "--ksize", a.ksize, "-o", "cmp", "--csv", "similarity.csv",
                *sigs], check=True)
subprocess.run(["sourmash", "plot", "--labels", "cmp"], check=True)

# ---- size-robust comparison, because Jaccard is confounded by assembly SIZE ------------------
# ⛔ JACCARD DIVIDES BY THE UNION, SO IT SCORES A SIZE DIFFERENCE AS DISTANCE. For |A| = 286 Mb
# against |B| = 800 Mb the ceiling is 286/800 = 0.36 EVEN IF A IS A PERFECT SUBSET OF B, and the
# same ceiling applies at the other end -- an unpurged 2,297 Mb assembly caps at ~0.35 against the
# same 800 Mb genome. Measured on the 23-genome panel already run: 641.8 Mb to 1,333.4 Mb, a 2.08x
# spread, so JL_DASH cannot exceed 0.48 against Salk_SRIb however related they are, while two
# similar-sized members sharing 90% score ~0.82. `similarity` is what WF-I's fold order consumes,
# so the failure is a wrong ordering with every job green.
#
# TWO MATRICES ARE WRITTEN, AND THE SECOND ONE IS THE REASON THE FIRST IS NOT ENOUGH:
#
#   containment.csv       |A n B| / |A|. ASYMMETRIC -- containment(A,B) != containment(B,A)
#                         whenever the two differ in size. Nothing is clustered or plotted from
#                         it, because a dendrogram over an asymmetric matrix is an artefact of
#                         whichever triangle the clustering happened to read.
#   max_containment.csv   max of the two directions. SYMMETRIC, so it CAN be clustered -- and it
#                         is, into its own heatmap and dendrogram. This is the tree to read for a
#                         panel whose members differ in size; the Jaccard tree beside it is the
#                         classic workflow's output and is kept for parity.
#
def write_newick(csv_name, out_name, units):
    """Cluster a symmetric similarity CSV and write the tree as Newick.

    ⛔ THE DENDROGRAM sourmash PUBLISHES IS A PICTURE, AND A PICTURE IS NOT A TREE. `sourmash plot`
    emits PNG or PDF and nothing else -- no branch lengths, no leaf order, nothing another tool can
    read. So the same linkage it draws is also written here, which is what iTOL, FigTree, ggtree and
    ete3 consume. ⚠ THIS IS SHARED BY BOTH MATRICES ON PURPOSE: the similarity tree and the
    max-containment tree are only comparable if they were built the same way, and two copies of this
    logic would be two chances to drift apart -- which is exactly the failure this workflow keeps
    finding elsewhere.

    ⚠ average linkage, and it is a CHOICE rather than a default worth hiding: `sourmash plot` itself
    clusters with scipy's default `single` linkage, which chains through intermediates and on a
    panel with haplotype pairs produces long ladders. `average` (UPGMA) is what a reader expects of
    a distance dendrogram. ⛔ THE PNG AND THE NEWICK ARE THEREFORE NOT THE SAME TREE, for either
    matrix, and the method is named in each output's label so nobody has to guess.

    ⚠ BRANCH LENGTHS ARE A DISSIMILARITY, NOT AN EVOLUTIONARY DISTANCE -- `units` names which one.
    """
    import csv as _csv

    import numpy as _np
    from scipy.cluster.hierarchy import linkage as _linkage
    from scipy.cluster.hierarchy import to_tree as _to_tree
    from scipy.spatial.distance import squareform as _squareform

    with pathlib.Path(csv_name).open(newline="", encoding="utf-8") as _fh:
        _rows = [r for r in _csv.reader(_fh) if r and not r[0].startswith("#")]
    _labels = _rows[0]
    _M = _np.array([[float(v) for v in r] for r in _rows[1:]], dtype=float)
    # ⛔ THE SHAPE IS ASSERTED BEFORE THE TREE IS BUILT. A matrix that is not square against its own
    # labels yields a tree with the wrong leaves and no error anywhere.
    if _M.shape != (len(_labels), len(_labels)):
        sys.exit(f"{csv_name} is {_M.shape}, not square against its {len(_labels)} labels -- "
                 f"refusing to build a tree from a matrix whose shape it cannot trust.")
    # ⛔ SYMMETRISED AND ZERO-DIAGONALLED BEFORE squareform, WHICH REFUSES OTHERWISE. Both matrices
    # are symmetric by construction, but the two triangles are computed separately and differ in the
    # last bits; `squareform` raises on that rather than rounding, and the raise is about "not
    # symmetric" with no hint that the asymmetry is 1e-16.
    _D = 1.0 - _M
    _D = (_D + _D.T) / 2.0
    _np.fill_diagonal(_D, 0.0)
    _Z = _linkage(_squareform(_D, checks=False), method="average")

    def _newick(_node, _parent_height):
        """One clade, with the branch length that reaches it from its parent."""
        _length = _parent_height - _node.dist
        if _node.is_leaf():
            return f"{_labels[_node.id]}:{_length:.6f}"
        _kids = ",".join(_newick(_c, _node.dist) for _c in (_node.get_left(), _node.get_right()))
        return f"({_kids}):{_length:.6f}"

    _root = _to_tree(_Z)
    # ⛔ PURE NEWICK, NO LEADING COMMENT. A bracketed `[...]` comment is legal Newick and was written
    # here first, but acceptance across the tools this file exists for -- iTOL, FigTree, ggtree -- is
    # uneven, and a tree some readers reject defeats the point of emitting it. The method and the
    # units live in the output's label and help text instead, where nothing can choke on them.
    pathlib.Path(out_name).write_text(
        f"({','.join(_newick(_c, _root.dist) for _c in (_root.get_left(), _root.get_right()))});\n",
        encoding="utf-8")
    print(f"wrote {out_name} ({len(_labels)} leaves, average linkage over {units})", file=sys.stderr)


# ⛔ THE SIMILARITY TREE IS WRITTEN UNCONDITIONALLY, BECAUSE similarity.csv ALWAYS EXISTS. The
# max-containment pair below is gated on `--containment`; this one is not, and the asymmetry is the
# whole reason the two were not comparable before: one tree left Galaxy as data and the other only
# as a PNG, so "do the two clusterings agree?" could not be asked of the outputs at all.
write_newick("similarity.csv", "similarity.newick", "1 - Jaccard similarity")

# ⚠ THE TWO TREES CAN DISAGREE, AND THAT IS THE POINT. Where they do, the Jaccard one is the one
# distorted by size. Neither is labelled "the" tree here.
if str(a.containment).lower() in ("true", "1", "yes"):
    subprocess.run(["sourmash", "compare", "--processes", a.processes, "--ksize", a.ksize, "--containment",
                    "--csv", "containment.csv", *sigs], check=True)
    subprocess.run(["sourmash", "compare", "--processes", a.processes, "--ksize", a.ksize, "--max-containment",
                    "-o", "maxc", "--csv", "max_containment.csv", *sigs], check=True)
    # `sourmash plot` writes <prefix>.matrix.png and <prefix>.dendro.png beside the input.
    subprocess.run(["sourmash", "plot", "--labels", "maxc"], check=True)

    # ⛔ THE DENDROGRAM sourmash PUBLISHES IS A PICTURE, AND A PICTURE IS NOT A TREE. `sourmash
    # plot` emits PNG or PDF and nothing else -- no branch lengths, no leaf order, nothing any
    # other tool can read. So the same linkage it draws is also written as NEWICK here, which is
    # what iTOL, FigTree, ggtree and ete3 consume. Same clustering, in a form that leaves Galaxy.
    #
    # ⚠ average linkage, and it is a CHOICE rather than a default worth hiding: `sourmash plot`
    # itself clusters with scipy's default `single` linkage, which chains through intermediates and
    # on a panel with haplotype pairs produces long ladders. `average` (UPGMA) is what a reader
    # expects of a distance dendrogram. ⛔ THE TWO THEREFORE DISAGREE: the PNG and this Newick are
    # not the same tree, and the Newick names its method in a comment so nobody has to guess.
    #
    # ⚠ 1 - max_containment IS A DISSIMILARITY, NOT AN EVOLUTIONARY DISTANCE. Branch lengths here
    # are in that unit. Read the tree as "shares content with", never as descent: on this panel a
    # haplotype pair scores ~0.62 while the 8x size pair scores 0.858, because neither haplotype
    # contains the other while a small genome IS largely inside a big one.
    write_newick("max_containment.csv", "max_containment.newick", "1 - max_containment")
else:
    # ⚠ OFF BY DEFAULT, AND NOTHING IS WRITTEN WHEN OFF -- deliberately, rather than emitting an
    # empty CSV. A zero-row matrix is the silent-success shape this project keeps finding: it
    # arrives as a real dataset, joins cleanly, and describes nothing. The outputs are declared
    # optional so Galaxy simply does not produce them.
    print("containment: not requested (containment=false); writing no matrix", file=sys.stderr)

# ---- publish the per-strain signatures ----------------------------------------------------------
# ⛔ THE CLASSIC WORKFLOW PUBLISHES THESE AND THIS PORT DID NOT -- WF-A's one real parity loss.
# `sourmash_sketch` mapped over the panel leaves one .sig dataset per strain, described in the
# classic's own port list as "BRC-reusable"; collapsing sketch and compare into ONE job (which is
# what the identifier assertion needs) turned them into work-dir files nobody could reach. A
# discovered collection gives them back without splitting the job in two again.
#
# ⛔ AND THE FILENAME IS THE POINT, WHICH IS WHY THIS COPIES RATHER THAN DISCOVERING `stage/`.
# Discovery takes each element identifier from the FILENAME, and `stage/` is deliberately keyed by
# INDEX -- pointing discovery at it hands back a collection keyed `0000`..`000N`, which joins
# against nothing downstream: not the sizes, not the self-pairs, not the relabel map, every one of
# which keys on the strain. A collection whose identifiers are ordinals is worse than no collection
# at all, because it looks like one.
#
# ⚠ SAFE ONLY BECAUSE OF THE IDENTIFIER GUARDS ABOVE, and this is now the second reason they exist:
# a duplicate, a path separator, a dot-alias or a leading dash would collide here or write outside
# the work dir, and all four are refused before anything is sketched.
#
# ⛔ THE INDEX PREFIX IS WHAT PUTS THIS COLLECTION IN PANEL ORDER, AND IT COSTS NOTHING. Discovery
# sorts on the FULL filename (`sort_key: filename`) but takes the element identifier from the
# `designation` REGEX GROUP -- so `0000_cs10.sig` sorts by its index while the element is still
# called `cs10`. Every other WF-A output is in panel order (`sizes`, `busco_summaries` and
# `fasta_index` are map-over collections; `self_pairs` and `relabel_map` are built line-per-strain
# off the identifier list), so an alphabetical `signatures` was the odd one out.
#
# ⚠ AN EARLIER VERSION OF THIS COMMENT SAID PANEL ORDER "CANNOT BE MADE TO MATCH", and that was
# wrong in a way worth recording: it is true only of the `sort_key` FIELD, whose every choice is
# name-derived. The filenames are ours to pick, which settles it -- and the claim had already been
# merged into two PRs before anyone tested the alternative. An impossibility claim needs a
# measurement as much as a number does.
#
# ⚠ MEASURED ON usegalaxy.org 26.1, 2026-09-08 -- history bbd44e69cb8906b51601636ef7c09ef9. A panel
# given as `cs10, PvP01, strain.2`, whose panel order and alphabetical order DIFFER (or the run
# could not have told the two hypotheses apart), came back as a 3-element list keyed
# `cs10, PvP01, strain.2` at indices 0,1,2 -- panel order, with the `0000_` prefix nowhere in the
# identifiers. An earlier run without the prefix (history bbd44e69cb8906b54a726562e523b991) is where
# the rest was measured and still holds: the declared `json` format is honoured, the elements are
# hidden while the collection is visible, `strain.2` survives the pattern with its dot, and the
# three sizes differ -- which is what says discovery found three files, not one file three times.
pathlib.Path("signatures").mkdir(exist_ok=True)
for _i, name in enumerate(ids):
    shutil.copyfile(f"stage/{_i:04d}.sig", f"signatures/{_i:04d}_{name}.sig")
# ⛔ COUNT WHAT LANDED RATHER THAN ASSUMING IT. Discovery publishes whatever it finds, so a short
# collection is a GREEN job with a genome missing -- the same silent loss the duplicate-identifier
# guard exists to prevent, one step further down. If two names ever collapse onto one file by a
# route those guards do not model (a case-insensitive filesystem, a future relaxation), this says so.
_n_sigs = len(list(pathlib.Path("signatures").glob("*.sig")))
if _n_sigs != len(ids):
    sys.exit(f"published {_n_sigs} signature file(s) for {len(ids)} strain(s): the discovered "
             f"`signatures` collection would be short a genome and the job would still succeed.")
print(f"published {_n_sigs} signature(s) to signatures/", file=sys.stderr)
