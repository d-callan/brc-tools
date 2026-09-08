# Run the full Pv4 pangenome A→K on the Galaxy instance (cluster-backed)

Agent instructions for executing the **full 8-strain Pv4** *P. vivax* pangenome
pipeline (Phases A→K) **through the Galaxy instance** — workflows are submitted to
the Galaxy URL via `planemo run` / bioblend, and **Galaxy dispatches the jobs onto
the new cluster** (`nekrut/cluster`) as its compute backend. We do **not** run the
bash impl on the cluster; `nekrut/Pv4-pangenome` (`v3/`) is only the reference
impl + config + ground truth to validate against.

Everything below targets the Galaxy instance (`$GALAXY_URL` + API key); the
cluster is invisible to the workflows — it's just where Galaxy's jobs land.

## 0. Ground rules
- **Each workflow phase runs in its own Galaxy history.** Every `planemo run`
  gets a distinct `--history_name` (e.g. `Pv4-full :: C align_chain_project`).
  Every bioblend-driven phase calls `create_history` once and runs all its tools
  there. **Never append one phase's jobs to another phase's history.**
- **All artifacts on the 4 TB SSD** (the volume Galaxy serves over NFS). Set
  `PV4_SSD` to its mount; use `$PV4_SSD/pv4_full/` as the single artifacts root
  for inputs, `planemo --output_directory`, bioblend download dirs, config files,
  BUSTED archives, the multiz/hub assembly, and all logs. Nothing on head-node
  local disk. Verify Galaxy's `file_path`/`job_working_directory`/`new_file_path`
  resolve onto the SSD/NFS mount before starting.

## 1. Read first (institutional knowledge — avoids re-debugging)
- `execution/galaxy_full_run_status.md` — full A→K test-panel run + every bug/fix.
- `execution/jk_relaxed_status.md` — J / K / relaxed-MSA specifics.
- `ONE-CLICK-READINESS.md` — per-workflow one-click status + gxformat2 limits.
- `workflows/{align_chain_project,msa,vcf_projection,ucsc_hub}/README.md`.
All fixes are already on `brc-tools` main (group_cds ref-gene; vcf_projection
versioned tool_ids; one-click WF-C; helpers removed). `git pull` first.

## 2. Galaxy + tools on the cluster
- Install the `brc-tools` local tools (`local_tool_conf.xml` entries) and ensure
  the IUC/toolshed tools resolve: kegalign, batched_lastz, ucsc_* (axtchain,
  chainsort, chainprenet, chainnet, netchainsubset, chainswap, axttomaf,
  gff3togenepred, genepredtobed, bedtobigbed, fatotwobit, hubcheck), chainStitchId,
  liftoff, mafft (rbc_mafft), pal2nal, trimal, iqtree3, hyphy_busted, bcftools_*,
  crossmap_vcf, pggb, sourmash, and the native collection ops
  (`__CROSS_PRODUCT_FLAT__`, `__FILTER_FROM_FILE__`, `__FILTER_FAILED_DATASETS__`,
  `__RELABEL_FROM_FILE__`).
- **Containers**: the cluster likely uses Singularity/Apptainer — confirm tool
  requirement resolution. (The docker-image issues in the single-node docs are
  irrelevant here.)
- **Scheduler**: jobs go to Slurm/Pulsar per the cluster config — the single-node
  4-worker bottleneck does NOT apply. **KegAlign → GPU partition.**

## 3. Inputs to stage (under `$PV4_SSD/pv4_full/inputs/`)
- `assemblies/{strain}.fa` for all 8 (accessions in `Pv4-pangenome/v3/species.conf`).
- `annotations/plasmodb-68/{anchor}.gff3` for the 3 anchors (PvW1, PAM, PvSY56);
  derive `{anchor}.bed12` + `{anchor}.isoforms.tsv` (gffread).
- **`{anchor}.gene.gff3`** per anchor: copy the anchor GFF3 and rename gene-level
  feature types (`protein_coding_gene`/`ncRNA_gene`/`pseudogene` → `gene`).
  Liftoff's default mode finds **0 genes** in the native types — this is required.
- Cohort VCF for J (real MalariaGEN *P. vivax* if available; else synthesize on
  REF_STRAIN).

## 4. Config
- Production config: `Pv4-pangenome/v3/species.conf` (8 strains; anchors
  PvW1/PAM/PvSY56; `MIN_INTACT_STRICT=7`, `RELAXED=5`). **Confirm `REF_STRAIN`**
  (defaulted to PvW1; the test used PvP01, which is not a full-panel anchor).
- WF-C one-click config files — generate onto the SSD:
  ```
  python3 brc-tools/execution/cluster/gen_wfc_config.py \
    $PV4_SSD/pv4_full/wfc_cfg \
    PvP01,PvW1,PAM,PvSY56,Sal-I,PvT01,PvC01,MHC087 PvW1,PAM,PvSY56
  ```
  → `self_pairs.txt` (8), `anchor_self_pairs.txt` (3), `relabel_map.tsv` (56 rows,
  one per directed non-self chain — the `A_A` diagonal is omitted so the row count
  matches the collection WF-C relabels, which is what `strict: true` requires).

## 5. Per-phase recipe (each its own history; gotchas baked in)
| Ph | Workflow / driver | Notes |
|----|-------------------|-------|
| A | `inventory` (sourmash) | all-vs-all → matrix |
| B | `softmask` | → softmasked FASTAs + `.sizes` |
| C | `align_chain_project` | **one-click** `planemo run` with the §4 config files + `gene.gff3` anchors. 28 KegAlign pairs (GPU) → 56 cleaned + 28 rbest chains; 21 projections |
| D | `pggb` | graph build |
| E | `consensus` | rbest overlap + graph edges + union-find → ortholog_table |
| F/G/H | `msa`/`trees`/`selection` | **drive tool-by-tool via bioblend `{"batch":True,"values":[hdca]}`** — the gxformat2 batch wiring collapses a collection. Run **twice**: strict (mi=7) and relaxed (mi=5), each its own history. Filter failed elements (`__FILTER_FAILED_DATASETS__`) before iqtree/hyphy. **hyphy tree input = IQ-TREE treefile, not the report.** |
| I | `multiz` | per-reference fold from C's axt outputs |
| J | `vcf_projection` | one-click (tool_ids versioned). `cleaned_chains` + `target_fastas` collections keyed by the same target id. Needs a cohort VCF |
| K | `ucsc_hub` | workflow builds tracks; **then assemble the hub tree out-of-band** with `execution/assemble_hub.py` + run `hubCheck`. Multiz s-lines must be `species.chrom`-named **and space-separated** before `maf_to_bigmaf_bed`; the workflow does NOT emit hub.txt/trackDb |

## 6. Recommended order
1. **Smoke-test the 5-strain test panel on the cluster first**
   (`Pv4-pangenome/v3/pipeline/test_data/species.conf` + `make_test_data.sh`) to
   confirm the cluster Galaxy, tool installs, GPU, and containers reproduce the
   known-green results in `execution/galaxy_full_run_status.md`.
2. Then full Pv4 (8 strains), phase by phase, each in its own history.
3. Validate counts/outputs against `Pv4-pangenome` ground truth where present.

## 7. Scale (8 strains)
28 KegAlign pairs → 56 directed chains; 21 anchor×query projections; thousands of
OGs; strict(mi7)+relaxed(mi5) MSA/tree/BUSTED sets; 8 multiz folds. GPU for
KegAlign. Budget shared-SSD space for the full set (≫ the test panel).
