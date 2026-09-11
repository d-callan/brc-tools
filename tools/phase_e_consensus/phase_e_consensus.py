#!/usr/bin/env python3
"""
Phase E — consensus ortholog table for the Pk v1 pipeline.

Builds a multigraph from per-anchor Phase C.4 classification.tsv files:
  - Each (strain, gene_id) pair is a graph node.
  - Edges connect orthologous genes; weight by source/intactness.
  - Connected components → orthogroups.
  - Labels: CORE-1:1 | CORE-VAR | FAMILY | PARTIAL | LINEAGE-SPECIFIC

Nodes are NATIVE genes. A projected annotation keeps the reference gene's id --
Liftoff writes PVPAM_130008300 into PvC01's GFF -- so a projection is resolved
back onto the native gene at the same locus (--gene-beds, --alias-overlap) and
contributes an EDGE. One that resolves to nothing is dropped, because adding it
as a node would mean a second node for a locus that already has one.

That is what went wrong before. Without the resolution step the 2026-06-12 run
(invocation cc39af39a106fd9e) put each physical gene in the graph twice: 82% of
populated cells held both a native and a projected id, max_copies was 2 in 83.4%
of groups, and since the labels key on max_copies only 21 of 5,817 orthogroups
came out CORE-1:1 -- 0.4%, for eight conspecific P. vivax strains. With the fix
the same inputs give 3,990 of 5,731, or 69.6%, matching what the rbest edges
alone produce (scripts/rbest_baseline.py: 68.6%). The projections now add real
linkage on top of the chains instead of inflating copy number: CORE-1:1 goes
slightly UP versus rbest alone and PARTIAL goes down.

Two columns warn about over-merging, and they catch different failures.

`clique` scores a group against the mutual support 1:1 evidence can actually
provide: sum over strain pairs of min(copies_a, copies_b). Scoring against every
possible pair instead would be wrong, because each gene can link to at most one
gene per other strain, so a group with m copies in each of k strains tops out at
(k-1)/(m*k-1) -- 0.462 for m=2, k=7, which is exactly where the multi-copy median
sat before this was corrected. Low `clique` means chained: the component hangs
together through a few links rather than mutual agreement.

`density` is the raw edge fraction over all member pairs. It is not comparable
across copy numbers, which is why it is not the threshold, but it is what exposes
a blob: OG000001 on the 2026-06-12 data holds 3,110 genes with 625 copies in one
strain and scores 0.525 clique against 0.0008 density. Sort by density to find
blobs, by clique to find chains.

Both are reported, NOT acted on -- union-find can only merge, never split, so a
wrong edge is permanent and an automatic response would compound it. Treat
FAMILY and CORE-VAR labels on low-scoring groups with suspicion.

Output: work/03_consensus/ortholog_table.tsv
  orthogroup_id, label, n_strains, max_copies, clique, {strain columns...}

Ported and parameterized from
  /media/anton/data/sandbox/Pv4/v3/scripts/phase_e_consensus.py
"""

import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


# TOGA2 loss symbols, ranked best-to-worst by TOGA2 itself in
# cesar_wrapper_constants.NUM_TO_CLASS: FI(8) > I(7) > PI(6) > UL(5) > L(4) >
# M(3) > PM(2) > PG(1) > PP(0). FI is "Fully Intact" and is TOGA2's STRONGEST
# call, not a weak one -- these weights were originally written against TOGA1,
# which has no FI, and it was previously scored 0.20. On the first verified
# TOGA2 pair FI accounted for 3906 of 6092 rescues, so under-scoring it
# discarded most of the evidence the rescue pass exists to produce.
WEIGHTS = {
    ('cesar2', 'FI'): 1.00,
    ('cesar2', 'I'):  1.00,
    ('cesar2', 'PI'): 0.70,
    ('cesar2', 'UL'): 0.40,
    ('cesar2', 'PG'): 0.40,
    ('cesar2', 'L'):  0.00,
    ('cesar2', 'M'):  0.00,
    ('liftoff', 'I'): 0.95,
}


def edge_weight(source, intactness):
    if source == 'liftoff':
        return 0.95
    if source == 'none':
        return 0.0
    return WEIGHTS.get((source, intactness), 0.10)


def normalize_gene_id(gid):
    if not gid or gid == 'None':
        return gid
    # _tN PlasmoDB transcript suffix
    m = re.match(r'^(.+)_t\d+$', gid)
    if m:
        return m.group(1)
    # .N transcript suffix
    m = re.match(r'^(.+)\.\d+$', gid)
    if m:
        return m.group(1)
    # _N Liftoff extra-copy suffix (small integer)
    m = re.match(r'^(.+)_(\d+)$', gid)
    if m and len(m.group(2)) <= 2 and not m.group(1).endswith('_'):
        return m.group(1)
    return gid


def load_gene_beds(spec):
    """{strain: {chrom: [(start, end, gene_id), ...]}} from per-strain gene BEDs.

    The BED's filename stem is the strain name, matching phase_e_rbest_overlap.
    """
    import glob as _glob
    idx: dict = {}
    paths = _glob.glob(spec) if any(c in spec for c in "*?[") else \
        [str(p) for p in Path(spec).iterdir()] if Path(spec).is_dir() else [spec]
    for path in sorted(paths):
        strain = Path(path).stem
        by_chrom: dict = defaultdict(list)
        with open(path) as fh:
            for ln in fh:
                f = ln.rstrip('\n').split('\t')
                if len(f) >= 4:
                    try:
                        by_chrom[f[0]].append((int(f[1]), int(f[2]), f[3]))
                    except ValueError:
                        continue
        for c in by_chrom:
            by_chrom[c].sort()
        if by_chrom:
            idx[strain] = dict(by_chrom)
    return idx


def resolve_to_native(gene_index, strain, chrom, start, end, min_overlap):
    """Best-overlapping native gene id at these coordinates, or None.

    A projected annotation keeps the REFERENCE gene's id -- Liftoff writes
    PVPAM_130008300 into PvC01's GFF -- so the projected id never equals the
    native one. On the 2026-06-12 run, zero of PvC01's 4,718 projected ids
    appear among its 6,769 native ids. Without this lookup each physical gene
    enters the graph twice, every group shows two "copies" per strain, and the
    labels, which key on max_copies, become artifacts: that run reported 21
    CORE-1:1 groups out of 5,817 where the rbest edges alone give 3,979 of 5,804.
    """
    best_id, best_ov = None, 0.0
    for ns, ne, gid in gene_index.get(strain, {}).get(chrom, ()):
        if ne <= start:
            continue
        if ns >= end:
            break
        ov = min(end, ne) - max(start, ns)
        if ov <= 0:
            continue
        r = min(ov / max(1, end - start), ov / max(1, ne - ns))
        if r > best_ov:
            best_ov, best_id = r, gid
    return best_id if best_ov >= min_overlap else None


class UnionFind:
    def __init__(self):
        self.parent: dict = {}

    def find(self, x):
        while self.parent.get(x, x) != x:
            self.parent[x] = self.parent.get(self.parent[x], self.parent[x])
            x = self.parent[x]
        self.parent.setdefault(x, x)
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def load_classifications(liftoff_dir, anchors, all_strains):
    """Yield (anchor, query, ref_gene_canonical, q_gene_canonical, source, intactness, weight)."""
    base = Path(liftoff_dir)
    for anchor in anchors:
        sub = base / f'{anchor}-as-ref'
        if not sub.exists():
            continue
        for q in all_strains:
            if q == anchor:
                continue
            cls_path = sub / f'{q}.classification.tsv'
            if not cls_path.exists():
                continue
            with open(cls_path) as fh:
                r = csv.DictReader(fh, delimiter='\t')
                for row in r:
                    rg = row.get('reference_gene_id', '')
                    qg = row.get('query_gene_id', '')
                    src = row.get('source', '')
                    intact = row.get('intactness', '')
                    if not qg or qg == 'None':
                        continue
                    w = edge_weight(src, intact)
                    if w == 0:
                        continue
                    yield anchor, q, rg, qg, src, intact, w


def reciprocal_overlap(a, b):
    s1, e1, _ = a
    s2, e2, _ = b
    ov = max(0, min(e1, e2) - max(s1, s2))
    if ov == 0:
        return 0.0
    return min(ov / max(1, e1 - s1), ov / max(1, e2 - s2))


def collapse_positions(gene_id_list, strain, node_pos):
    """Group gene IDs by overlapping query position; return list of position-clusters."""
    recs = [(gid, node_pos.get(f'{strain}#{gid}')) for gid in gene_id_list]
    positioned = [(p, gid) for gid, p in recs if p is not None]
    no_pos = [gid for gid, p in recs if p is None]
    if not positioned:
        return [[g] for g in dict.fromkeys(no_pos)]
    by_chr: dict = defaultdict(list)
    for (chrom, s, e), gid in positioned:
        by_chr[chrom].append((s, e, gid))
    clusters = []
    for chrom, lst in by_chr.items():
        lst.sort()
        used = [False] * len(lst)
        for i in range(len(lst)):
            if used[i]:
                continue
            s_i, e_i, g_i = lst[i]
            grp = [g_i]
            used[i] = True
            for j in range(i + 1, len(lst)):
                if used[j]:
                    continue
                s_j, e_j, g_j = lst[j]
                if s_j > e_i:
                    break
                ov = max(0, min(e_i, e_j) - max(s_i, s_j))
                if ov and min(ov / (e_i - s_i + 1), ov / (e_j - s_j + 1)) >= 0.2:
                    grp.append(g_j)
                    used[j] = True
            clusters.append(grp)
    for g in dict.fromkeys(no_pos):
        clusters.append([g])
    return clusters


def main():
    ap = argparse.ArgumentParser(description="Phase E: consensus ortholog table")
    ap.add_argument('--liftoff_dir', required=True, help='work/02d_merged/')
    ap.add_argument('--rbest', required=True, help='work/03_consensus/rbest_edges.tsv')
    ap.add_argument('--graph', default=None,
                    help='graph_edges.tsv from a pggb/odgi graph. Optional: the pipeline runs '
                         'without a graph at all, and on the 2026-06-12 data the graph '
                         'contributed zero edges, so the table was already chains + projections.')
    ap.add_argument('--anchors', required=True, help='Space-separated anchor strain list')
    ap.add_argument('--strains', required=True, help='Space-separated all-strain list')
    ap.add_argument('--ref', required=True, help='Reference strain name')
    ap.add_argument('--output', required=True, help='Output TSV path')
    ap.add_argument('--gene-beds', dest='gene_beds', default=None,
                    help='Per-strain native gene BEDs (dir or glob; filename stem = strain). '
                         'Used to map a projected gene back onto the native gene at the same '
                         'locus. Without this every projection enters the graph as its own node.')
    ap.add_argument('--alias-overlap', dest='alias_overlap', type=float, default=0.5,
                    help='Reciprocal overlap needed to call a projected gene the same gene as a '
                         'native one (default 0.5). On the 2026-06-12 data 36.9%% of projections '
                         'sit at >=1.0 and 16.0%% overlap no native gene at all; the rest form a '
                         'gradient, so this is the knob that decides how much of it is merged.')
    ap.add_argument('--keep-unresolved-projections', dest='edges_only',
                    action='store_false', default=True,
                    help='Add a projection that resolves to no native gene as a node of its own, '
                         'rather than dropping it. This was the old behaviour and it is what made '
                         'the 2026-06-12 table unusable -- an unresolvable projection is a second '
                         'node for a locus that already has one, so max_copies (which the labels '
                         'key on) counts it as an extra copy. The default is to treat projections '
                         'purely as edges between native genes.')
    ap.add_argument('--split-many2many', action='store_true',
                    help="Do not union on a projection whose orthology_class is many2many. "
                         "TOGA2 emits that as a warning that the genes should probably NOT "
                         "collapse into one group; union-find can never undo the merge.")
    args = ap.parse_args()

    anchors = args.anchors.split()
    all_strains = args.strains.split()

    gene_index = load_gene_beds(args.gene_beds) if args.gene_beds else {}
    if gene_index:
        print(f'Native gene BEDs: {len(gene_index)} strains '
              f'({sum(len(v) for v in gene_index.values())} contigs)')
    else:
        print('WARNING: no --gene-beds given. Projected genes keep their reference-derived '
              'ids, so each physical gene enters the graph twice and max_copies -- which the '
              'labels key on -- is inflated.', file=sys.stderr)

    print('Loading per-anchor classifications...', flush=True)
    raw = list(load_classifications(args.liftoff_dir, anchors, all_strains))
    print(f'  raw edges: {len(raw)}')

    uf = UnionFind()
    pos_records: dict = defaultdict(list)
    node_pos: dict = {}
    used_edges: set = set()
    aliased_to_native = 0
    unresolved_projections = 0
    many2many_skipped = 0

    # Seed positions from classification.tsv query_chrom/start/end columns
    base = Path(args.liftoff_dir)
    for anchor in anchors:
        sub = base / f'{anchor}-as-ref'
        if not sub.exists():
            continue
        for q in all_strains:
            if q == anchor:
                continue
            cls_path = sub / f'{q}.classification.tsv'
            if not cls_path.exists():
                continue
            with open(cls_path) as fh:
                r = csv.DictReader(fh, delimiter='\t')
                for row in r:
                    rg = normalize_gene_id(row.get('reference_gene_id', ''))
                    qg = normalize_gene_id(row.get('query_gene_id', ''))
                    src = row.get('source', '')
                    intact = row.get('intactness', '')
                    if not qg or qg == 'None':
                        continue
                    if src == 'none' or edge_weight(src, intact) == 0:
                        continue
                    chrom = row.get('query_chrom', '')
                    start_s = row.get('query_start', '')
                    end_s = row.get('query_end', '')
                    coords = None
                    if chrom and start_s and end_s:
                        try:
                            coords = (chrom, int(start_s), int(end_s))
                        except ValueError:
                            coords = None

                    # A projected gene carries the REFERENCE gene's id, so it can
                    # never equal the native id for the same locus. Resolve it to
                    # the native gene by position; keep the projected id only when
                    # nothing native is there, which is real signal -- the
                    # projection found a gene the annotation missed.
                    if gene_index and coords:
                        native = resolve_to_native(gene_index, q, coords[0], coords[1],
                                                   coords[2], args.alias_overlap)
                        if native:
                            qg = normalize_gene_id(native)
                            aliased_to_native += 1
                        else:
                            unresolved_projections += 1
                            if args.edges_only:
                                continue

                    ocls = (row.get('orthology_class') or '').strip()
                    if args.split_many2many and ocls == 'many2many':
                        many2many_skipped += 1
                        continue

                    a_node = f'{anchor}#{rg}'
                    q_node = f'{q}#{qg}'
                    uf.union(a_node, q_node)
                    used_edges.add(frozenset((a_node, q_node)))
                    if coords:
                        pos_records[(q, coords[0])].append((coords[1], coords[2], q_node))
                        node_pos[q_node] = coords

    # Interval-based aliasing (90% reciprocal overlap → same physical gene)
    aliases_merged = 0
    for key, recs in pos_records.items():
        recs.sort()
        for i in range(len(recs)):
            si, ei, ni = recs[i]
            for j in range(i + 1, len(recs)):
                sj, ej, nj = recs[j]
                if sj > ei:
                    break
                if uf.find(ni) == uf.find(nj):
                    continue
                if reciprocal_overlap(recs[i], recs[j]) >= 0.9:
                    uf.union(ni, nj)
                    aliases_merged += 1
    print(f'  position aliases merged: {aliases_merged}')
    if gene_index:
        tot = aliased_to_native + unresolved_projections
        pct = (100.0 * aliased_to_native / tot) if tot else 0.0
        print(f'  projections resolved to a native gene: {aliased_to_native:,} of {tot:,} '
              f'({pct:.1f}%) at >={args.alias_overlap} reciprocal overlap')
        kept = 'dropped (--projections-as-edges-only)' if args.edges_only \
            else 'kept as their own nodes'
        print(f'  projections with no native gene at that locus: {unresolved_projections:,} '
              f'({kept})')
    if args.split_many2many:
        print(f'  many2many projections skipped: {many2many_skipped:,}')

    # Incorporate rbest chain edges
    rbest_edges_added = 0
    if Path(args.rbest).exists():
        with open(args.rbest) as fh:
            r = csv.DictReader(fh, delimiter='\t')
            for row in r:
                sa = row.get('strain_a', '')
                ga = normalize_gene_id(row.get('gene_a', ''))
                sb = row.get('strain_b', '')
                gb = normalize_gene_id(row.get('gene_b', ''))
                if sa and ga and sb and gb:
                    uf.union(f'{sa}#{ga}', f'{sb}#{gb}')
                    used_edges.add(frozenset((f'{sa}#{ga}', f'{sb}#{gb}')))
                    rbest_edges_added += 1
    print(f'  rbest chain edges: {rbest_edges_added}')

    # Incorporate graph co-membership edges
    graph_edges_added = 0
    if args.graph and Path(args.graph).exists():
        with open(args.graph) as fh:
            r = csv.DictReader(fh, delimiter='\t')
            for row in r:
                sa = row.get('strain_a', '')
                ga = normalize_gene_id(row.get('gene_a', ''))
                sb = row.get('strain_b', '')
                gb = normalize_gene_id(row.get('gene_b', ''))
                if sa and ga and sb and gb:
                    uf.union(f'{sa}#{ga}', f'{sb}#{gb}')
                    used_edges.add(frozenset((f'{sa}#{ga}', f'{sb}#{gb}')))
                    graph_edges_added += 1
    print(f'  graph co-membership edges: {graph_edges_added}')

    # Clique completeness needs the edges themselves, not just the components.
    # A group spanning k strains built from 1:1 evidence should carry k*(k-1)/2
    # undirected edges; far fewer means the component was chained together
    # through a handful of links rather than being mutually supported. Union-find
    # cannot undo a merge, so this is measured and reported -- it is the signal
    # to look at before trusting a FAMILY or CORE-VAR label.
    comp_edges: dict = Counter()
    for e in used_edges:
        comp_edges[uf.find(next(iter(e)))] += 1

    # Connected components
    comps: dict = defaultdict(set)
    for node in uf.parent:
        comps[uf.find(node)].add(node)
    print(f'  orthogroup count (all connected components): {len(comps)}')

    N_ALL = len(all_strains)
    rows_out = []
    for cid, nodes in sorted(comps.items(), key=lambda kv: (-len(kv[1]), min(kv[1]))):
        per_strain: dict = defaultdict(list)
        for n in sorted(nodes):      # `nodes` is a set: sort so runs are diffable
            parts = n.split('#', 1)
            if len(parts) == 2:
                per_strain[parts[0]].append(parts[1])
        if len(per_strain) < 2:
            continue
        strain_clusters = {s: collapse_positions(gs, s, node_pos) for s, gs in per_strain.items()}
        present_strains = [s for s in all_strains if s in per_strain]
        n_strains = len(present_strains)
        max_copies = max(len(c) for c in strain_clusters.values())
        # How much of the mutual support this group COULD have does it actually
        # have? Raw density (edges / n*(n-1)/2) is not comparable across copy
        # numbers: the evidence is 1:1, so each gene can link to at most one gene
        # per other strain, and a group with m copies in each of k strains tops
        # out at (k-1)/(m*k-1) -- 0.462 for m=2, k=7. Scoring against n*(n-1)/2
        # would therefore flag every multi-copy group as chained, which is a false
        # alarm on exactly the FAMILY and CORE-VAR groups the number exists for.
        #
        # The achievable maximum is sum over strain pairs of min(copies_a,
        # copies_b). For single-copy groups that is k*(k-1)/2 and this reduces to
        # plain density; above that it stays a real 0..1 score.
        counts = sorted(len(c) for c in strain_clusters.values())
        expected_edges = sum(min(counts[i], counts[j])
                             for i in range(len(counts)) for j in range(i + 1, len(counts)))
        observed = comp_edges.get(cid, 0)
        clique = round(min(observed / expected_edges, 1.0), 3) if expected_edges else 1.0
        # Raw density too: the normalised score is comparable across groups, which
        # is what makes a threshold usable, but it flatters a blob. OG000001 --
        # 3,110 genes, 625 copies in one strain -- scores 0.525 normalised and
        # 0.001 raw. Both are true; only the pair tells the whole story.
        n_nodes = len(nodes)
        pairs = n_nodes * (n_nodes - 1) // 2
        density = round(observed / pairs, 4) if pairs else 1.0
        if n_strains == N_ALL and max_copies == 1:
            label = 'CORE-1:1'
        elif n_strains == N_ALL and max_copies >= 2:
            label = 'CORE-VAR'
        elif max_copies >= 3:
            label = 'FAMILY'
        elif n_strains <= 2:
            label = 'LINEAGE-SPECIFIC'
        else:
            label = 'PARTIAL'
        row = {
            'orthogroup_id': f'OG{len(rows_out)+1:06d}',
            'label': label,
            'n_strains': n_strains,
            'max_copies': max_copies,
            'clique': clique,
            'density': density,
        }
        for s in all_strains:
            if s in strain_clusters:
                row[s] = ','.join('|'.join(c) for c in strain_clusters[s])
            else:
                row[s] = '-'
        rows_out.append(row)

    print(f'  multi-strain orthogroups: {len(rows_out)}')
    labels = Counter(r['label'] for r in rows_out)
    for k, n in labels.most_common():
        print(f'    {k}: {n}  ({100.0 * n / max(1, len(rows_out)):.1f}%)')

    cl = sorted(r['clique'] for r in rows_out)
    if cl:
        ragged = sum(1 for c in cl if c < 0.9)
        print(f'  clique completeness: median {cl[len(cl) // 2]:.3f}, '
              f'{ragged:,} groups ({100.0 * ragged / len(cl):.1f}%) below 0.9')
        print('    Those were chained together rather than mutually supported. '
              'They are reported, not split -- union-find cannot undo a merge.')
        print('    `clique` is scored against the support 1:1 evidence can actually '
              'provide; `density` is the raw edge fraction. Sort by density to find '
              'blobs, by clique to find chains.')
        print('  Compare against scripts/rbest_baseline.py on the same rbest edges: '
              'a healthy run should land near it, not far above on CORE-VAR.')

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = ['orthogroup_id', 'label', 'n_strains', 'max_copies', 'clique', 'density'] \
        + all_strains
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows_out)
    print(f'Wrote {out}')


if __name__ == '__main__':
    main()
