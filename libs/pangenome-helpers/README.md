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

## Modules (planned)

| Module | Source tool(s) | Purpose |
| --- | --- | --- |
| `pangenome_helpers.manifest` | `group_cds_by_og`, `pansn_rename` | Collection manifest parsing + validation |
| `pangenome_helpers.orthology` | `group_cds_by_og` | Orthogroup filtering, min-intact logic |
| `pangenome_helpers.cds` | `group_cds_by_og` | CDS/protein extraction orchestration |
| `pangenome_helpers.triage` | `phase_c2_triage` | R1–R8 rule engine |
| `pangenome_helpers.merge` | `phase_c4_merge` | Merge summaries + fallback decisions |
| `pangenome_helpers.consensus` | `phase_e_consensus` | Orthogroup consensus building |
| `pangenome_helpers.overlap` | `phase_e_rbest_overlap` | Reciprocal-best block projection orchestration |
| `pangenome_helpers.graph_edges` | `phase_e_graph_edges` | Graph edge emission |
| `pangenome_helpers.maf` | `process_maf`, `multiz_order` | MAF folding + order generation |
| `pangenome_helpers.pansn` | `pansn_rename` | FASTA header renaming orchestration |
| `pangenome_helpers.anchors` | `anchor_prep` | Anchor input assembly |
| `pangenome_helpers.hub` | `build_genomes_txt`, `build_trackdb`, `build_hub_bb` | UCSC hub assembly |

Initial release (0.1.0) ships the scaffolding; subsequent minor releases will fill in each module with the extracted logic and tests.
