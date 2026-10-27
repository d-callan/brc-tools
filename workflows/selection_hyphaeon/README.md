# WF-GH `selection_hyphaeon` — HyphAeon selection scan + UCSC/TACC assets (v0)

Replaces the two-stage WF-G (iqtree3 per-OG phylogeny) + WF-H (`hyphy_busted`
mapped over collections) with a single tree-free workflow:

1. **HyphAeon** (`hyphaeon_meme`, `hyphaeon_epistasis`, `hyphaeon_dms`,
   `hyphaeon_splits` — all stock IUC wrappers, `--no-tree`) maps directly
   over the per-orthogroup codon-alignment collection.
2. **`brc-hyphaeon-bundle`** (UDT; source `tools/hyphaeon_bundle/`) projects
   codon coordinates through the reference BED12 and emits two UCSC track
   BEDs, one AutoSql schema, a per-OG summary table, and the TACC
   comparative bundle. All statistics are descriptive (site-level
   HyphAeon p/q-values and site-call counts).
3. `sort` + `ucsc_bedtobigbed` (IUC) convert the two BEDs to BigBeds.

## v0 scope — no recombination handling

This version performs **no recombination scan or partitioning**: every
orthogroup alignment is analysed whole and every collection element is
keyed by the bare `{og}` identifier.

**Scientific caveat:** for highly recombining organisms (Neisseria,
pneumococcus, H. pylori, coronaviruses), uncorrected recombinant tracts can
inflate the site-level selection signal — interpret `selection_sites`
accordingly. For effectively clonal organisms the impact is negligible.

Recombination-aware partition-and-merge is planned **upstream** in the
HyphAeon monorepo (`rhizaeon partition` emitting FASTAs + `manifest.json`,
`hyphaeon --site-offset`, `hyphaeon merge`/`verdict`) — see
`specs/recombination-partition.md` there. When it lands, v1 of this
workflow adds the scan+partition steps and a `recombination_tracts`
track; the v0 bundle deliberately carries no recombination inputs — the
v1 wiring is written against whatever contract upstream settles on.

The same applies to **gene-level statistics**: v0 ships no combined
gene p/q-value. A Cauchy/ACAT omnibus over within-gene site p-values is
unvalidated (site p's are dependent — shared tree, shared alignment,
and in v0 unpartitioned recombinant tracts), and the BUSTED gene-level
signal was previously pulled from upstream prod for exactly this reason.
What a defensible gene-level claim looks like is v1 upstream design work;
`gene_selection.tsv` keeps only ortholog context, tallies of HyphAeon
site calls, and passthrough HyphAeon metrics -- no per-OG extreme-value
statistics (e.g. `min(p)`), which would read as gene-level claims.

The interim Python partitioner stays on the development branch only and
is not part of this change set; it serves as scaffolding for
cross-checking the upstream port.

## Inputs

| Input | Type | Source | Notes |
|---|---|---|---|
| `codon_alignments` | `list` collection of `fasta` | WF-F `pal2nal_codon` | Element identifier = orthogroup id; FASTA headers = ortholog-table column names |
| `ortholog_table` | `tabular` | WF-E `consensus` | `orthogroup_id`, `label`, `n_strains`, `max_copies`, `clique`, `density`, one column per strain |
| `ref_bed12` | `bed` | WF-A `anchor_prep` | Reference gene models; `name` = gene ID |
| `ref_sizes` | `tabular` | WF-A `brc-samtools-faidx` | 2-column `chrom<TAB>length` |
| `ref_column` | text, default `REF` | — | Ortholog-table column of the reference strain; also its FASTA header and the DMS `focal_taxon`. Matches WF-F `ref_name` |
| `bundle_id` | text (no default) | operator | Strain-panel identifier for `manifest.json` (e.g. `panel5`) |
| `bundle_version` | text, default `v1` | operator | Version for `manifest.json` |

Organism-agnostic by construction: reference-gene resolution is set-membership
against the `ref_bed12` names (composite `id1|id2,id3` cells resolved like
`group_cds_by_og`), the `label` column supplies orthogroup classes, and no
gene-ID prefixes or species names appear anywhere.

## Outputs

- `meme_json_collection`, `epistasis_json_collection`, `dms_json_collection`,
  `splits_newick_collection`, `splits_csv_collection` — per-orthogroup raw
  results (`epistasis_graphml_collection` also emitted).
- `selection_sites_bb` (BED9+7), `orthogroup_bb` (BED12) — the two
  comparative tracks consumed by WF-K (`workflows/ucsc_hub`).
- `gene_selection_table` — per-OG joined summary (ortholog + selection +
  epistasis + DMS metrics; descriptive only — no recombination fields
  and no gene-level p/q column in v0).
- `tacc_comparative_bundle` — `manifest.json` + `tables/` + per-OG
  `genes/*.json.gz` payloads for the BRC comparative-genomics UI.
- `qc_report` — MultiQC custom-content table.

## UDTs used (deploy via `scripts/udt_deploy.py`)

| UDT | Renders | Why it exists |
|---|---|---|
| `brc-hyphaeon-bundle` | `tools/hyphaeon_bundle/hyphaeon_bundle.py` | Codon->genome projection, track BEDs, TACC packaging |

It is a `CONVERTED`-class entry in `scripts/check_udt_definitions.py::UDT_RATIONALE`
and is generated verbatim from `tools/` by `scripts/build_selection_udts.py`.

## Consumed by

- **WF-K `ucsc_hub`** takes `selection_sites_bb` and `orthogroup_bb` as
  the inputs to its reference-only `brc_pangenome_select` composite
  (the old `busted_strict`/`busted_relaxed` inputs and
  `build_hub_bb` step are retired; the `selection_genes` and
  `recombination_tracts` subtracks return with v1).
- **TACC bundle** is staged out-of-band to the static hosting tree;
  `tables/gene_selection.json` and `genes/{og}.json.gz` drive the BRC
  organism-page comparative drawer.
