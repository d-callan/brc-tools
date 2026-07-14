# Libraries

This directory contains reusable Python libraries for genomic data processing, each published independently to PyPI and Bioconda.

## `genome-io`

Lightweight library for parsing and processing genomic data formats (BED, GFF3, FASTA, MAF, orthology tables, etc.).

- **PyPI**: `genome-io`
- **Bioconda**: `genome-io`
- **Status**: Stable (v0.1.0+)
- **Dependencies**: Minimal (pyfaidx optional)

See `genome-io/README.md` for details.

## `pangenome-helpers`

Pangenome workflow orchestration helpers built on top of `genome-io` (manifest loading, orthogroup filtering, CDS grouping, triage/merge stages, UCSC hub builders).

- **PyPI**: `pangenome-helpers` (scaffolding in progress)
- **Bioconda**: `pangenome-helpers`
- **Status**: 0.1.0 scaffold (logic landing in subsequent releases)
- **Dependencies**: `genome-io>=0.1.0`

---

## Development

Each library has its own `pyproject.toml`, tests, and documentation. To work on a specific library:

```bash
cd genome-io
pip install -e ".[test]"
pytest
```

To run all tests:

```bash
for lib in genome-io pangenome-helpers; do
  [ -d "$lib" ] && (cd "$lib" && pytest) || true
done
```
