# libs/ TODO — genome-io + pangenome-helpers

Findings from a full review of both libraries (2026-07-14). Ranked by severity.
Update this file as items are fixed (move to "Done" or delete the line) or as
new issues are found — this is meant to be a living reference, not a one-time
report.

---

## 🔴 Correctness bugs (fix first)

### 1. Minus-strand splice-site detection is broken
**File:** `pangenome-helpers/src/pangenome_helpers/triage.py:255-265` (`_has_splice_issue`)

`_has_splice_issue` always applies the plus-strand coordinate formula
(`exons_sorted[i][1] + 1` for the donor site) and never reverse-complements
the extracted bases, regardless of `gene.strand`.

The original `tools/phase_c2_triage/phase_c2_triage.py::get_splice_sites`
swaps which exon boundary is used for donor/acceptor *and*
reverse-complements when `strand == '-'`. On minus-strand, multi-exon genes
the ported version now reads incorrect genomic positions, which will make
rule **R6** fire spuriously (or fail to fire) unpredictably.

**Why untested:** `tests/test_triage.py` only exercises single-exon genes on
an unspecified/positive strand, so the bug produces no visible test failure.

**Fix:** Port the strand-aware branch from the original `get_splice_sites`
(swap donor/acceptor exon boundary based on strand, reverse-complement the
extracted bases for `-` strand) into `_has_splice_issue`. Add a regression
test with a multi-exon minus-strand gene and a deliberately non-canonical
splice site to confirm R6 fires correctly in both orientations.

### 2. `is_variant_antigen` drops the original label
**Files:** `pangenome-helpers/src/pangenome_helpers/hub.py` and
`pangenome-helpers/src/pangenome_helpers/selection.py` (identical duplicate
function in both)

Original (`tools/build_hub_bb/build_hub_bb.py::is_variant_antigen`) returns
the **full original label** (e.g. `"PIR_beta"`) when a family keyword
substring-matches; the ported version returns just the matched family key
(e.g. `"PIR"`), discarding label detail in the `gene_family` BED12+5 column.

It also iterates a `set` of family keywords — if a label happened to match
more than one keyword, the returned family is hash-order dependent (could
differ between runs/interpreters).

**Fix:** Return `label` (matching original semantics) instead of `fam`.
Iterate a fixed-order sequence (tuple/list) instead of a `set` regardless,
for determinism. Fix in `selection.py` only once `hub.py`'s duplicate copy
is deleted (see #3).

---

## 🟠 Redundancy — duplicate logic across modules

Global rule: avoid duplicate logic; use the same implementation for similar
cases. The following violate that and should be consolidated.

### 3. `bh_fdr` / `qval_to_rgb` / `qval_to_score` / `rgb_to_int` / `is_variant_antigen` defined twice
**Files:** `hub.py` (dead — not exported in `__init__.py`) and `selection.py`
(canonical — exported)

These five functions are byte-for-byte identical in both files. `hub.py`'s
copies have no internal callers and aren't re-exported, so they're pure dead
code left over from before `selection.py` existed.

**Fix:** Delete the "Selection helpers shared with build_hub_bb" section from
`hub.py` entirely; `selection.py` is the canonical home.

### 4. Three different `normalize_gene_id` implementations, three different behaviors
- `genome_io/gff.py:48-62` (`normalize_gene_id`) — strips `_t\d+`, `.\d+`,
  *and* small numeric extra-copy suffixes (`_NN` where `NN` ≤ 2 digits).
- `pangenome_helpers/triage.py:83-89` (private `_normalize_gene_id`) — only
  strips the extra-copy suffix. Missing the `_t\d+`/`.\d+` cases.
- `pangenome_helpers/merge.py:261-267` (private `_normalize_gene_id`) — same
  subset as triage's version, also missing the transcript-suffix cases.

`cds.py` and `consensus.py` correctly import and use
`genome_io.gff.normalize_gene_id`. `triage.py` and `merge.py` instead
maintain private, weaker reimplementations. A Liftoff transcript ID with a
`_t1` suffix will normalize differently depending on which module touches
it — a real cross-module correctness/consistency risk, not just style.

**Fix:** Delete the private `_normalize_gene_id` in both `triage.py` and
`merge.py`; import `normalize_gene_id` from `genome_io.gff` instead. Add a
test asserting `triage.py`/`merge.py` behavior matches `genome_io.gff` for
the `_t\d+` / `.\d+` cases (this is exactly the kind of regression a shared
import prevents).

### 5. Sequence primitives reimplemented in `triage.py`
**File:** `pangenome_helpers/triage.py`

`_revcomp`, `_has_internal_stop`, `STOP_CODONS`, `_slice` /
`_extract_cds_sequence` duplicate `genome_io.sequence.revcomp`,
`has_internal_stop`, `STOP_CODONS`, `extract_sequence`.

`genome_io.sequence.extract_sequence(fasta, chrom, start, end, strand)`
already works fine when `fasta` is a plain `dict[str, str]` (which is what
`triage.py` uses) — `fasta[chrom]` returns a `str`, and `str(...)` on a `str`
is a no-op. There is no functional reason for the private copies.

**Fix:** Import `revcomp`, `has_internal_stop`, `STOP_CODONS`,
`extract_sequence` from `genome_io.sequence` and delete the private copies
in `triage.py`. Note: `_extract_cds_sequence` builds a full CDS from
multiple `(start, end)` segments — check whether `genome_io.sequence`
needs a small addition (e.g. a segments-only CDS builder that doesn't
require the 6-tuple `(chrom, start, end, strand, phase, parent)` shape used
by `extract_cds`) rather than force-fitting the existing signature.

### 6. `genome_io/maf.py` — `parse_blocks` and `iter_maf_blocks` duplicate block-splitting logic
Eager (`parse_blocks`) and lazy (`iter_maf_blocks`) block-splitting share
near-identical logic maintained as two separate implementations.
`tests/test_maf.py::test_parse_blocks_and_iter_blocks_agree` is a good
regression guard, but doesn't remove the duplication/maintenance risk if one
is changed without the other.

**Fix (low priority, more invasive):** Consider implementing `parse_blocks`
as `header, blocks = <extract header lines>, list(iter_maf_blocks(...))` so
there's a single source of truth for block-splitting.

### 7. `genome_io/gff.py` — two independent GFF attribute parsers
`parse_gff_attributes_to_dict` (semicolon/`=`-aware, used almost everywhere)
vs. `parse_gene_id` (ad hoc `rfind("ID=")`, used only by `gff_to_bed_rows`).

**Fix:** Have `gff_to_bed_rows` use `parse_gff_attributes_to_dict(...).get("ID")`
(stripping a leading `gene-` prefix if still needed) instead of the separate
`parse_gene_id` helper. Delete `parse_gene_id` once unused.

---

## 🟡 Minor bugs / dead code / API surface

### 8. Stray docstring-looking string sits after imports in `__init__.py`
**File:** `pangenome_helpers/__init__.py:17`

```python
from .merge import (...)
"""Pangenome workflow orchestration helpers built on genome-io."""   # <-- dead statement, not the module docstring

from .anchors import ...
```

Because it's not the first statement in the file, this is just an inert
expression statement — `pangenome_helpers.__doc__` is `None`, not this
string. Leftover from an import reordering.

**Fix:** Move this string literal to line 1 of the file (before any
imports) so it actually becomes the module docstring.

### 9. `genome_io/io.py::open_maybe_gz` — redundant branch
```python
if "w" in mode and mode.endswith("t"):
    return open(path, mode)
return open(path, mode)
```
Both branches resolve to the identical call. The whole `if "w" in mode...`
block is redundant with the unconditional fallthrough line.

**Fix:** Delete the redundant `if` block; keep just the final
`return open(path, mode)`.

### 10. Public API surface inconsistencies in `genome_io`
- `genome_io.maf.reorder_block` is used cross-package
  (`pangenome_helpers/maf.py` imports the submodule directly:
  `from genome_io import maf as maf_utils`), but unlike its siblings
  (`find_ref_index`, `emit_bed_record`, `species_of`, ...) it is **not**
  re-exported from `genome_io/__init__.py`'s `__all__`. Anyone doing
  `from genome_io import reorder_block` will fail even though it's a
  legitimate, actively-used function.
- `genome_io.orthology.safe_name` is defined but never exported, and never
  imported anywhere under `libs/`. A near-duplicate exists standalone in
  `tools/group_cds_by_og/group_cds_by_og.py` (outside `libs/`, so out of
  scope for this doc, but worth knowing about if `group_cds_by_og` is ever
  ported into `pangenome_helpers.cds`).

**Fix:** Add `reorder_block` to `genome_io/__init__.py` imports/`__all__`
for consistency. Either export `safe_name` too (if it's meant to be public)
or delete it if it's truly unused.

---

## 🧪 Testing gaps

### 11. `selection.py` has zero test coverage
No `tests/test_selection.py` exists despite `selection.py` being a fully new
module (file parsing, BUSTED JSON extraction from dir/tarball, BED-row
builders, BH-FDR stats). This is exactly the kind of logic that regresses
silently without tests.

**Fix:** Add `tests/test_selection.py` covering:
- `load_sizes`, `load_bed12` (including the isoform-suffix-stripping / first-isoform-wins behavior)
- `load_ortholog_table` (including the `ref_column` missing → `ValueError` case, `|`-separated gene lists, `gene_prefix` filtering)
- `extract_busted_pvalues` for both directory and `.tar.gz` sources
- `build_selection_bed_rows` / `build_orthogroup_bed_rows` (sorting, chrom-size clipping)
- `bh_fdr` / `qval_to_rgb` / `qval_to_score` / `rgb_to_int` / `is_variant_antigen` (once #2 is fixed, assert full label is preserved)

### 12. No minus-strand test case for triage splice rule
Directly related to bug #1 — add a multi-exon minus-strand `GeneRecord`
fixture with a non-canonical splice site to `tests/data/triage/` and assert
`R6_splice` fires (and that a canonical minus-strand case does *not* fire).

### 13. `genome-io` test suite is the stronger of the two
For reference/comparison: every `genome-io` module has direct unit tests
including edge cases (empty attrs, trailing semicolons, gzip magic-byte
sniffing, strand-aware CDS extraction). Use it as the bar when adding tests
to `pangenome-helpers`.

---

## ✅ What's working well (don't break these)

- Clean two-layer architecture: `genome-io` = format-parsing primitives,
  `pangenome-helpers` = per-workflow-stage orchestration built on top.
- Most `pangenome-helpers` modules (`anchors`, `manifest`, `overlap`,
  `graph_edges`, `cds`, `consensus`, `orthology`) are thin and correctly
  delegate to `genome-io` rather than reimplementing primitives — `triage.py`
  and `merge.py` are the outliers (see #4, #5).
- Consistent dataclass/result-object conventions (`AnchorPrepResult`,
  `MafProcessResult`, `TriageResult`, `MergeOutputs`, `Phase2Output`) instead
  of bare tuples/dicts almost everywhere.
- `pyproject.toml` dependency declaration
  (`pangenome-helpers` → `genome-io>=0.1.0`) is correct and versioned.

---

## Suggested fix order

1. #1 minus-strand splice bug (correctness, silent failure risk)
2. #2 `is_variant_antigen` label bug (correctness, data quality)
3. #3 dedupe `hub.py` / `selection.py` (quick, mechanical)
4. #4 unify `normalize_gene_id` usage in `triage.py`/`merge.py`
5. #5 unify sequence primitives in `triage.py`
6. #11 add `test_selection.py`
7. #8, #9, #10 mechanical cleanups (batchable together)
8. #6, #7 lower priority, more invasive — revisit if touching `maf.py`/`gff.py` for other reasons anyway
