# anchor_prep

Builds the two annotation files TOGA2 needs from one anchor GFF3:

- **BED12** — protein-coding transcripts, `gffread --bed`, transcript id in column 4.
- **isoforms TSV** — `gene<TAB>transcript`, one row per `mRNA` feature, read from
  the `Parent` and `ID` attributes of column 9.

Both come out of a single call on a single GFF3 by design. TOGA2 ties each
projection back to a gene through the isoforms table, so the transcript ids in
the BED12 and in the table have to agree; one step is what guarantees they do.

Used by WF-A (`workflows/inventory/inventory.gxwf.yml`), which maps it over the
`anchor_gene_gff3s` collection and publishes both outputs for WF-C2.

## Verifying the isoforms output by hand

The transformation is small enough to reproduce with stock Galaxy tools, which is
worth doing once on a new annotation source:

| step | tool | setting |
|------|------|---------|
| 1 | `Filter1` | `c3=='mRNA'` |
| 2 | `Cut1` | `c9` |
| 3 | `Convert characters` | `;` → tab |
| 4 | `Convert characters` | `=` → tab |
| 5 | `Cut1` | `c4,c2` |
| 6 | `sort1` | c1 then c2 |

On PlasmoDB-68 this matches `anchor_prep` row for row (PvW1 6,075 / PAM 6,462 /
PvSY56 5,328). It is a check rather than a replacement: it assumes `ID=` precedes
`Parent=` in column 9 and that neither value contains `=` or `;`, none of which
GFF3 guarantees, and it cannot produce the BED12.

Note that Unix `cut -f4,2` emits fields in ascending order and would give you the
columns backwards; Galaxy's `Cut1` honours the order you ask for.

## Deployment note

The tool declares two conda requirements (`gffread`, `python`), so Galaxy resolves
it to a **merged `mulled-v1-<hash>` environment** rather than to the single-package
`__gffread@…` env. On a host where `conda_auto_install` fails to build that merged
env, the job dies with **exit 1 and completely empty stderr** — the failure is
recorded only in `<job_dir>/working/conda_activate.log`:

    EnvironmentLocationNotFound: Not a conda environment: .../envs/mulled-v1-<hash>

Build it explicitly, using the exact env name from the generated `tool_script.sh`:

    conda create -y -p $GALAXY_TOOL_DEPS/_conda/envs/mulled-v1-<hash> \
        -c conda-forge -c bioconda gffread=0.12.7 python=3.12

If compute nodes mount `tool_deps` over NFS, they can take a few minutes to see a
newly created env; until the client attribute cache expires, jobs landing on those
nodes fail with the same message while jobs on the head node succeed.
