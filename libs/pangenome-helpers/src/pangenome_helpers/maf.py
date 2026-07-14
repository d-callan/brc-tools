"""End-to-end orchestration around genome_io.maf + genome_io.multiz."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

from genome_io import maf as maf_utils
from genome_io import multiz


class MafProcessResult(tuple):
    """Named result for processed MAF stats (kept, dropped)."""

    __slots__ = ()

    @property
    def kept(self) -> int:  # pragma: no cover - property is trivial
        return self[0]

    @property
    def dropped(self) -> int:  # pragma: no cover
        return self[1]


def process_maf_file(
    ref_species: str,
    input_path: str | Path,
    output_path: str | Path,
) -> MafProcessResult:
    """Filter/reorder/sort a MAF (maftoBigMaf pre-processing).

    Returns ``MafProcessResult(kept, dropped_no_ref)``.
    """

    with open(input_path) as fh:
        header, blocks = maf_utils.parse_blocks(fh)
    kept: list[tuple[str, int, int, list[str]]] = []
    dropped = 0
    for block in blocks:
        idx = maf_utils.find_ref_index(block, ref_species)
        if idx is None:
            dropped += 1
            continue
        if idx != 0:
            block = maf_utils.reorder_block(block, idx)
        coords = maf_utils.ref_coords(block, ref_species)
        if coords is None:
            dropped += 1
            continue
        kept.append((coords[0], coords[1], coords[2], block))
    kept.sort(key=lambda t: (t[0], t[1]))
    with open(output_path, "w") as out:
        for line in header:
            out.write(line)
        for _, _, _, block in kept:
            for line in block:
                out.write(line)
            out.write("\n")
    return MafProcessResult((len(kept), dropped))


def derive_multiz_order(
    compare_csv: str | Path,
    hinge: str,
    queries: Iterable[str],
) -> list[str]:
    """Return queries ordered by descending similarity to the hinge."""

    labels, data = multiz.load_matrix(compare_csv)
    sims = multiz.similarities_to_hinge(labels, data, hinge)
    ordered = multiz.order_queries(list(queries), sims)
    if not ordered:
        raise ValueError("no query strains available for ordering")
    return ordered
