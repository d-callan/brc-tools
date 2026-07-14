# pangenome-helpers CLI

Command-line interface for pangenome workflow orchestration helpers.

## Installation

```bash
pip install -e .
```

This installs the `pangenome-helpers` command-line tool.

## Usage

```bash
pangenome-helpers --help
pangenome-helpers <command> --help
```

## Commands

### Phase C.2: Triage

Filter Liftoff genes and flag those needing CESAR2 fallback.

```bash
pangenome-helpers triage \
  <liftoff_gff> \
  <query_fasta> \
  <reference_bed> \
  <output_dir> \
  <query_name> \
  [--family-list FAMILY_LIST] \
  [--core-identity-min 0.95] \
  [--core-coverage-min 0.90] \
  [--family-identity-min 0.85] \
  [--subtelomere-bp 100000]
```

**Arguments:**
- `liftoff_gff`: Liftoff GFF3 file with gene annotations
- `query_fasta`: Query genome FASTA file
- `reference_bed`: Reference genes BED file
- `output_dir`: Output directory for results
- `query_name`: Query strain/genome name

**Options:**
- `--family-list`: TSV file mapping gene IDs to family names
- `--core-identity-min`: Minimum sequence identity for core genes (default: 0.95)
- `--core-coverage-min`: Minimum coverage for core genes (default: 0.90)
- `--family-identity-min`: Minimum identity for family genes (default: 0.85)
- `--subtelomere-bp`: Subtelomeric flank size in bp (default: 100000)

**Outputs:**
- `triage.tsv`: Triage decisions for each gene (R1-R8 rules)
- `needs_cesar2.bed`: BED file of genes flagged for CESAR2 fallback
- `liftoff_clean.gff3`: Filtered GFF3 with only passing genes
- `summary.json`: Summary statistics and thresholds used

**Example:**
```bash
pangenome-helpers triage \
  liftoff.gff3 \
  query.fa \
  reference.bed \
  ./triage_out \
  query_strain \
  --family-list families.tsv
```

---

### Phase C.4: Merge

Merge Liftoff clean annotations with TOGA2 projections.

```bash
pangenome-helpers merge \
  <query> \
  <liftoff_clean_gff> \
  <loss_summary_tsv> \
  <orthology_tsv> \
  <query_annotation_bed> \
  <query_genes_bed> \
  <reference_bed> \
  <output_dir>
```

**Arguments:**
- `query`: Query strain name
- `liftoff_clean_gff`: Liftoff clean GFF3 (from Phase C.2)
- `loss_summary_tsv`: TOGA loss summary TSV
- `orthology_tsv`: TOGA orthology TSV
- `query_annotation_bed`: Query annotation BED (from Liftoff)
- `query_genes_bed`: Query genes BED (from TOGA)
- `reference_bed`: Reference genes BED
- `output_dir`: Output directory

**Outputs:**
- `classification.tsv`: Classification rows combining Liftoff + TOGA
- `merged.gff3`: Merged GFF3 with both sources

---

### Phase E: Consensus

Build consensus orthogroups from Liftoff classifications and edges.

```bash
pangenome-helpers consensus \
  <liftoff_dir> \
  <anchors>... \
  --strains <strains>... \
  --ref-strain <ref_strain> \
  [--rbest-edges RBEST_EDGES] \
  [--graph-edges GRAPH_EDGES] \
  <output_tsv>
```

**Arguments:**
- `liftoff_dir`: Directory containing Liftoff classification TSVs
- `anchors`: Anchor strain names (space-separated)
- `--strains`: All strain names (space-separated)
- `--ref-strain`: Reference strain name
- `output_tsv`: Output consensus orthogroup TSV

**Options:**
- `--rbest-edges`: JSON file of reciprocal-best edges
- `--graph-edges`: JSON file of graph co-membership edges

**Outputs:**
- TSV with orthogroup_id, label (CORE-1:1, CORE-VAR, FAMILY, etc.), and per-strain gene columns

---

### Phase E: Reciprocal-best Edges

Compute reciprocal-best overlap edges from chain projections.

```bash
pangenome-helpers rbest-edges \
  <chains_pattern> \
  <annotations_pattern> \
  [--min-overlap 0.90] \
  <output_json>
```

**Arguments:**
- `chains_pattern`: Glob pattern for `*.chain` files (e.g., `chains/*.rbest.chain`)
- `annotations_pattern`: Glob pattern for `*.bed` files
- `output_json`: Output edges JSON

**Options:**
- `--min-overlap`: Minimum fractional overlap on both projections (default: 0.90)

**Output format:**
```json
[
  {
    "strain_a": "strainA",
    "gene_a": "geneA",
    "strain_b": "strainB",
    "gene_b": "geneB",
    "overlap_a": "0.950",
    "overlap_b": "0.920"
  }
]
```

---

### Phase E: Graph Edges

Compute graph co-membership edges from PGGB path data.

```bash
pangenome-helpers graph-edges \
  <paths_tsv> \
  <annotations_pattern> \
  --strains <strains>... \
  <output_json>
```

**Arguments:**
- `paths_tsv`: ODGI `paths --haplotypes` output TSV
- `annotations_pattern`: Glob pattern for `*.bed` files
- `--strains`: Strain names to include (space-separated)
- `output_json`: Output edges JSON

**Output format:**
```json
[
  {
    "strain_a": "sampleA",
    "gene_a": "geneA",
    "strain_b": "sampleB",
    "gene_b": "geneB",
    "path_id": "contig1",
    "overlap": "1.000"
  }
]
```

---

### Hub Building

Generate UCSC hub manifest and track database files.

```bash
pangenome-helpers hub \
  <genomes_metadata_tsv> \
  <output_dir> \
  [--hub-name HUB_NAME] \
  [--hub-email HUB_EMAIL]
```

**Arguments:**
- `genomes_metadata_tsv`: Genomes metadata TSV
- `output_dir`: Output directory

**Options:**
- `--hub-name`: Hub name (default: from metadata)
- `--hub-email`: Hub contact email (default: from metadata)

**Outputs:**
- `genomes.txt`: UCSC hub genomes manifest
- `trackDb.txt`: UCSC hub track database

---

### PanSN FASTA Renaming

Rename FASTA headers with PanSN (Pangenome Sequence Naming) prefixes.

```bash
pangenome-helpers pansn-rename \
  <input_fasta> \
  <output_fasta> \
  <sample> \
  [--haplotype 1] \
  [--delimiter '#'] \
  [--gzip]
```

**Arguments:**
- `input_fasta`: Input FASTA file
- `output_fasta`: Output FASTA file
- `sample`: Sample name

**Options:**
- `--haplotype`: Haplotype number (default: 1)
- `--delimiter`: PanSN delimiter (default: `#`)
- `--gzip`: Gzip output

**Example:**
```bash
pangenome-helpers pansn-rename input.fa output.fa sample1 --haplotype 2
# >contig1 description -> >sample1#2#contig1 description
```

---

### MAF Processing

Filter, reorder, and sort a MAF (Multiple Alignment Format) file.

```bash
pangenome-helpers process-maf \
  <ref_species> \
  <input_maf> \
  <output_maf>
```

**Arguments:**
- `ref_species`: Reference species name (for reordering)
- `input_maf`: Input MAF file
- `output_maf`: Output MAF file

**Behavior:**
- Filters blocks without reference species
- Reorders blocks to put reference first
- Sorts blocks by reference coordinates

---

### Multiz Ordering

Derive multiz query order from sourmash compare similarity matrix.

```bash
pangenome-helpers multiz-order \
  <compare_csv> \
  <hinge> \
  <queries>... \
  <output_txt>
```

**Arguments:**
- `compare_csv`: Sourmash compare CSV matrix
- `hinge`: Hinge strain name (reference for ordering)
- `queries`: Query strain names (space-separated)
- `output_txt`: Output file (one strain per line)

**Behavior:**
- Sorts queries by descending similarity to hinge
- Missing queries go last

---

### Anchor Prep

Prepare anchor inputs for TOGA2 (filter BED12, emit isoforms).

```bash
pangenome-helpers anchor-prep \
  <gff> \
  <raw_bed> \
  <output_bed> \
  <output_isoforms>
```

**Arguments:**
- `gff`: GFF3 file
- `raw_bed`: Raw gffread BED12 file
- `output_bed`: Output filtered BED
- `output_isoforms`: Output isoforms TSV (gene_id<TAB>transcript_id)

**Behavior:**
- Filters to protein-coding genes
- Rewrites BED name column to gene IDs
- Emits gene-transcript pairs

---

### CDS Grouping

Group CDS sequences by orthogroup.

```bash
pangenome-helpers group-cds \
  <ortho_table> \
  <ref_strain> \
  <ref_gff> \
  <ref_fasta> \
  <query_gff_manifest> \
  <query_fasta_manifest> \
  <output_json> \
  [--min-intact 2]
```

**Arguments:**
- `ortho_table`: Ortholog table TSV
- `ref_strain`: Reference strain name
- `ref_gff`: Reference GFF3
- `ref_fasta`: Reference FASTA
- `query_gff_manifest`: Query GFF manifest (id<TAB>path)
- `query_fasta_manifest`: Query FASTA manifest (id<TAB>path)
- `output_json`: Output orthogroup sequences JSON

**Options:**
- `--min-intact`: Minimum intact strains (default: 2)

**Output format:**
```json
[
  {
    "orthogroup_id": "OG000001",
    "reference_gene_id": "gene1",
    "cds": {
      "ref": "ATGCGTAA...",
      "strainA": "ATGCGTAA...",
      "strainB": "ATGCGTAA..."
    },
    "proteins": {
      "ref": "MRVI...",
      "strainA": "MRVI...",
      "strainB": "MRVI..."
    }
  }
]
```

---

## Exit Codes

- `0`: Success
- `1`: Error (check stderr for details)

## Examples

### Complete Phase C.2 triage workflow

```bash
pangenome-helpers triage \
  liftoff.gff3 \
  query.fa \
  reference.bed \
  ./c2_triage \
  query_strain \
  --family-list families.tsv \
  --core-identity-min 0.95 \
  --core-coverage-min 0.90
```

### Build consensus orthogroups with edges

```bash
# Compute edges
pangenome-helpers rbest-edges \
  chains/*.rbest.chain \
  annotations/*.bed \
  --min-overlap 0.90 \
  rbest_edges.json

pangenome-helpers graph-edges \
  graph_paths.tsv \
  annotations/*.bed \
  --strains ref strainA strainB \
  graph_edges.json

# Build consensus
pangenome-helpers consensus \
  liftoff_dir \
  ref \
  --strains ref strainA strainB \
  --ref-strain ref \
  --rbest-edges rbest_edges.json \
  --graph-edges graph_edges.json \
  consensus.tsv
```

### Prepare for UCSC hub

```bash
# Rename FASTA headers
pangenome-helpers pansn-rename \
  query.fa \
  query.pansn.fa \
  query_strain

# Process MAF
pangenome-helpers process-maf \
  ref_species \
  alignments.maf \
  alignments.sorted.maf

# Derive multiz order
pangenome-helpers multiz-order \
  compare.csv \
  ref_strain \
  strainA strainB strainC \
  multiz_order.txt

# Build hub
pangenome-helpers hub \
  genomes.tsv \
  ./hub_output \
  --hub-name "My Pangenome" \
  --hub-email "contact@example.com"
```

---

## Troubleshooting

### Command not found

Ensure the package is installed:
```bash
pip install -e .
```

### File not found errors

Check that file paths are correct and use absolute paths when possible:
```bash
pangenome-helpers triage \
  /absolute/path/to/liftoff.gff3 \
  /absolute/path/to/query.fa \
  ...
```

### Permission errors

Ensure the output directory is writable:
```bash
mkdir -p output_dir
chmod 755 output_dir
```

### Memory issues with large files

For very large FASTA files, consider processing in chunks or using streaming approaches. File a GitHub issue if you encounter memory problems.

---

## See Also

- [pangenome-helpers README](README.md)
- [genome-io documentation](../genome-io/README.md)
- [BRC Analytics workflows](../../workflows/)
