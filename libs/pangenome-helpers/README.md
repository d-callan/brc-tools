# pangenome-helpers

Pipelines for comparative genomics lean on a lot of repeatable orchestration: loading manifests, filtering orthogroups, triaging Liftoff projections, merging consensus tables, stitching UCSC hubs, etc. **pangenome-helpers** collects that orchestration into a reusable Python library built on top of [`genome-io`](../genome-io/README.md).

Goals:

- Keep Galaxy tool wrappers thin (argparse + logging only).
- Provide pytest-covered functions for every pangenome workflow stage.
- Ship on PyPI + Bioconda so workflows outside this repo can reuse the logic.

## Install

From PyPI (future):
```bash
pip install pangenome-helpers
```

From source (development):
```bash
cd libs/pangenome-helpers
pip install -e ".[test]"
```

## Run tests

```bash
cd libs/pangenome-helpers
pytest
```

## Command-line interface

A comprehensive CLI is available for running pangenome orchestration functions directly:

```bash
pangenome-helpers --help
pangenome-helpers triage --help
pangenome-helpers consensus --help
# ... etc
```

See [CLI.md](CLI.md) for complete documentation of all commands and examples.

## Modules (available today)

| Module | Source tool(s) | Purpose |
| --- | --- | --- |
| `pangenome_helpers.manifest` | `group_cds_by_og`, `pansn_rename` | Collection manifest parsing + validation |
| `pangenome_helpers.cds` | `group_cds_by_og` | CDS/protein extraction orchestration (yields per-OG FASTAs) |
| `pangenome_helpers.pansn` | `pansn_rename` | FASTA header renaming orchestration |
| `pangenome_helpers.maf` | `process_maf`, `multiz_order` | MAF filtering + fold-order generation |
| `pangenome_helpers.anchors` | `anchor_prep` | Anchor BED12 filtering + isoform TSV emission |
| `pangenome_helpers.overlap` | `phase_e_rbest_overlap` | Reciprocal-best block projection edges |
| `pangenome_helpers.graph_edges` | `phase_e_graph_edges` | PGGB path co-membership edges |
| `pangenome_helpers.consensus` | `phase_e_consensus` | Orthogroup consensus builder (UnionFind + labels) |
| `pangenome_helpers.triage` | `phase_c2_triage` | R1–R8 rule engine + summaries |
| `pangenome_helpers.orthology` | `group_cds_by_og` | Orthogroup filtering utilities shared across stages |
| `pangenome_helpers.merge` | `phase_c4_merge` | Merge Liftoff + TOGA annotations, emit classification TSV + merged GFF |
| `pangenome_helpers.hub` | `build_genomes_txt`, `build_trackdb`, `build_hub_bb` | UCSC browser hub manifest/track builders |
| `pangenome_helpers.selection` | `build_hub_bb` | Selection BED builders + BUSTED parsing + visualization helpers |
| `pangenome_helpers.phase_c2` | `phase_c2_triage` | Triage orchestration + TSV/BED/GFF/JSON reporting |

## Modules (planned / in progress)

None at this time. All core workflow orchestration helpers are implemented and tested.
