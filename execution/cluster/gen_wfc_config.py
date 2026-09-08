#!/usr/bin/env python3
"""Generate the WF-C (align_chain_project) one-click config inputs for a panel.

Usage: gen_wfc_config.py <out_dir> <STRAIN1,STRAIN2,...> <ANCHOR1,ANCHOR2,...>
Emits into <out_dir> (put it on the shared SSD artifacts root):
  self_pairs.txt         one X_X per strain  (FILTER_FROM_FILE remove self pairs)
  anchor_self_pairs.txt  one A_A per anchor  (projection grid self-cells)
  relabel_map.tsv        A_B<TAB>A.B per DISTINCT ordered strain pair (Phase E ids)
These are the panel-specific config the native cross-product/filter/relabel
built-ins need (no native primitive computes the self-cross diagonal; the
join separator can't be '.').

The relabel map omits the A_A diagonal so its row count matches the collection
WF-C relabels; strict mode checks that count first and refuses the run
otherwise. tools/collection_relabel_map/relabel_map.py is the other producer of
this file and must agree; scripts/check_relabel_strict.py checks both.
"""
import itertools
import os
import sys

out, strains_s, anchors_s = sys.argv[1], sys.argv[2], sys.argv[3]
strains = [s for s in strains_s.split(",") if s]
anchors = [a for a in anchors_s.split(",") if a]
os.makedirs(out, exist_ok=True)
with open(f"{out}/self_pairs.txt", "w") as f:
    f.writelines(f"{s}_{s}\n" for s in strains)
with open(f"{out}/anchor_self_pairs.txt", "w") as f:
    f.writelines(f"{a}_{a}\n" for a in anchors)
with open(f"{out}/relabel_map.tsv", "w") as f:
    # No diagonal: WF-C drops the A_A cells before relabelling, and strict mode
    # compares row count against element count.
    for a, b in itertools.product(strains, strains):
        if a != b:
            f.write(f"{a}_{b}\t{a}.{b}\n")
n = len(strains)
print(f"strains={n} anchors={len(anchors)}")
print(f"  self_pairs.txt: {n}")
print(f"  anchor_self_pairs.txt: {len(anchors)}")
print(f"  relabel_map.tsv: {n*(n-1)} rows ({n}x{n} less the diagonal); "
      f"directed non-self chains = {n*(n-1)}")
