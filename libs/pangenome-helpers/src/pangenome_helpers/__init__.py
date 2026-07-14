"""Pangenome workflow orchestration helpers built on genome-io."""

from .anchors import AnchorPrepResult, prepare_anchor_inputs
from .cds import OrthogroupSequences, iter_orthogroup_sequences
from .consensus import ConsensusResult, build_consensus_table, summarize_labels
from .graph_edges import compute_graph_edges
from .manifest import (
    ManifestError,
    ensure_matching_collections,
    ensure_reference_not_in_queries,
    load_manifest_map,
)
from .maf import MafProcessResult, derive_multiz_order, process_maf_file
from .pansn import PansnError, rename_fasta
from .overlap import RbestEdge, compute_rbest_edges
from .triage import (
    GeneRecord,
    TriageResult,
    TriageSettings,
    parse_liftoff_gff,
    read_family_list,
    read_reference_bed,
    run_triage,
)

__all__ = [
    "AnchorPrepResult",
    "ConsensusResult",
    "OrthogroupSequences",
    "TriageResult",
    "TriageSettings",
    "GeneRecord",
    "RbestEdge",
    "ManifestError",
    "PansnError",
    "MafProcessResult",
    "load_manifest_map",
    "ensure_matching_collections",
    "ensure_reference_not_in_queries",
    "rename_fasta",
    "prepare_anchor_inputs",
    "process_maf_file",
    "derive_multiz_order",
    "iter_orthogroup_sequences",
    "compute_rbest_edges",
    "compute_graph_edges",
    "build_consensus_table",
    "summarize_labels",
    "parse_liftoff_gff",
    "read_reference_bed",
    "read_family_list",
    "run_triage",
]
