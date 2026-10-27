#!/usr/bin/env python3
"""Generate a per-assembly UCSC track-hub ``trackDb.txt``.

Authored for the BRC pangenome pipeline (LOCAL.md K.5). Emits, in order:

1. ``track {name}_multiz``           - standalone, ``type bigMaf``
2. ``track brc_pangenome_chains``    - composite, ``type bigChain`` + N sub-tracks
3. ``track brc_pangenome_annot``     - composite, ``type bigBed 12`` + sub-tracks
4. ``track brc_pangenome_select``    - reference strain only, ``type bigBed 12``

bigMaf and bigChain CANNOT share a composite: UCSC composites require all
members to be of a single ``type``. They are therefore emitted as separate
top-level stanzas (the bigMaf standalone, the bigChain in its own composite).

Dependencies: Python standard library only.
"""

import argparse
import sys


def parse_pairs(items):
    """Parse ``ACC=Label`` strings into an ordered list of (acc, label) tuples.

    A bare ``ACC`` (no ``=``) keeps the accession as its own label.
    """
    out = []
    for item in items or []:
        if "=" in item:
            acc, label = item.split("=", 1)
        else:
            acc, label = item, item
        out.append((acc.strip(), label.strip()))
    return out


def label_for(acc, pairs, default=None):
    for a, lab in pairs:
        if a == acc:
            return lab
    return default if default is not None else acc


def maf_stanza(out, *, assembly, strain, maf_url, species_order, labels, group, html):
    """Standalone bigMaf track (NOT in a composite)."""
    name = "{}_multiz".format(strain.replace("-", "_").replace(" ", "_"))
    out.append("track {}".format(name))
    if html:
        out.append("html {}".format(html))
    out.append("shortLabel {} multiz".format(strain))
    out.append("longLabel  {} multi-z alignment".format(strain))
    out.append("type bigMaf")
    out.append("bigDataUrl {}".format(maf_url))
    out.append("group {}".format(group))
    out.append("visibility pack")
    if species_order:
        out.append("speciesOrder {}".format(" ".join(species_order)))
    if labels:
        rendered = " ".join('{}="{}"'.format(a, lab) for a, lab in labels)
        out.append("speciesLabels {}".format(rendered))
    out.append("")


def chains_stanza(out, *, assembly, targets, target_labels, group, html, chains_dir):
    """Composite bigChain track with one sub-track per target assembly."""
    out.append("track brc_pangenome_chains")
    if html:
        out.append("html {}".format(html))
    out.append("compositeTrack on")
    out.append("shortLabel Pangenome chains")
    out.append("longLabel  Pairwise chain alignments ({} targets)".format(len(targets)))
    out.append("type bigChain")
    out.append("group {}".format(group))
    out.append("visibility hide")
    out.append("")
    for tgt in targets:
        lab = label_for(tgt, target_labels)
        out.append("    track chain_to_{}".format(tgt))
        if html:
            out.append("    html {}".format(html))
        out.append("    parent brc_pangenome_chains off")
        out.append("    shortLabel chain to {}".format(lab))
        out.append("    longLabel  Chain alignment from {} to {}".format(assembly, lab))
        out.append("    type bigChain {}".format(tgt))
        out.append(
            "    bigDataUrl {}/{}_to_{}.bigChain.bb".format(chains_dir, assembly, tgt)
        )
        out.append(
            "    linkDataUrl {}/{}_to_{}.bigChain.link.bb".format(
                chains_dir, assembly, tgt
            )
        )
        out.append("    visibility hide")
        out.append("")


def annot_stanza(out, *, anchors, group, html):
    """Composite bigBed 12 annotation track, one sub-track per anchor strain."""
    out.append("track brc_pangenome_annot")
    if html:
        out.append("html {}".format(html))
    out.append("compositeTrack on")
    out.append("shortLabel Pangenome annot")
    out.append("longLabel  Gene projections (Liftoff + TOGA2) from anchor strains")
    out.append("type bigBed 12")
    out.append("group {}".format(group))
    out.append("visibility pack")
    out.append("")
    for acc, lab in anchors:
        out.append("    track annot_from_{}".format(lab))
        if html:
            out.append("    html {}".format(html))
        out.append("    parent brc_pangenome_annot off")
        out.append("    shortLabel annot from {}".format(lab))
        out.append("    longLabel  Genes projected from {} via Liftoff+TOGA2".format(lab))
        out.append("    type bigBed 12")
        out.append("    bigDataUrl annot_from_{}.bb".format(lab))
        out.append("    visibility dense")
        out.append("")


def select_stanza(out, *, group, html):
    """Comparative-genomics composite (reference strain only).

    Sub-tracks fed by WF-GH ``selection_hyphaeon`` BigBeds. The parent
    declares ``type bigBed`` (no field count) because the sub-tracks carry
    heterogeneous schemas; each sub-track declares its own ``type bigBed
    N +`` matching its AutoSql, so UCSC renders the extra columns.

    Names are load-bearing: the out-of-band hub staging materializes each
    ``bigDataUrl`` from the WF-GH output of the same base name -- the file
    name (``bb``) therefore tracks the WF-GH output name, not the subtrack
    label. The ``recombination_tracts`` subtrack is re-added when WF-GH
    regains recombination handling (upstream ``rhizaeon partition``).
    """
    out.append("track brc_pangenome_select")
    if html:
        out.append("html {}".format(html))
    out.append("compositeTrack on")
    out.append("shortLabel Comparative")
    out.append("longLabel  HyphAeon selection and orthogroup membership")
    out.append("type bigBed")
    out.append("group {}".format(group))
    out.append("visibility hide")
    out.append("")
    subs = [
        ("selection_sites", "Selection (sites)",
         "HyphAeon codon-resolution episodic selection sites",
         "bigBed 9 +", "selection_sites"),
        ("orthogroup_membership", "Orthogroups",
         "Orthogroup membership per reference gene",
         "bigBed 12", "orthogroup"),
    ]
    for name, short, long, btype, bb in subs:
        out.append("    track {}".format(name))
        if html:
            out.append("    html {}".format(html))
        out.append("    parent brc_pangenome_select off")
        out.append("    shortLabel {}".format(short))
        out.append("    longLabel  {}".format(long))
        out.append("    type {}".format(btype))
        out.append("    bigDataUrl {}.bb".format(bb))
        out.append("    visibility dense")
        out.append("")


def build(args):
    pairs = parse_pairs(args.strain_label)
    # species order / labels for the bigMaf: all strains in the panel.
    species_order = [s for s, _ in pairs] if pairs else []
    labels = pairs

    # chain / annot targets default to every other assembly in the panel.
    all_accs = [s for s, _ in pairs]
    targets = [a for a in all_accs if a != args.assembly]
    anchors = parse_pairs(args.anchor)

    maf_url = args.maf if args.maf else "{}.multiz.maf.bb".format(args.strain)

    out = []
    maf_stanza(
        out,
        assembly=args.assembly,
        strain=args.strain,
        maf_url=maf_url,
        species_order=species_order,
        labels=labels,
        group=args.group,
        html=args.maf_html,
    )
    chains_stanza(
        out,
        assembly=args.assembly,
        targets=targets,
        target_labels=pairs,
        group=args.group,
        html=args.chains_html,
        chains_dir=args.chains_dir,
    )
    annot_stanza(out, anchors=anchors, group=args.group, html=args.annot_html)
    if args.is_reference:
        select_stanza(out, group=args.group, html=args.select_html)

    text = "\n".join(out).rstrip("\n") + "\n"
    if args.output and args.output != "-":
        with open(args.output, "w") as fh:
            fh.write(text)
    else:
        sys.stdout.write(text)


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Generate a per-assembly UCSC trackDb.txt for the BRC pangenome hub."
    )
    p.add_argument("--assembly", "--accession", required=True,
                   help="Assembly accession of THIS genome (e.g. GCA_000002415.2).")
    p.add_argument("--strain", required=True,
                   help="Strain name for THIS genome (e.g. Sal-I).")
    p.add_argument("--strain-label", "--strain-labels", nargs="+", default=[],
                   metavar="ACC=LABEL",
                   help="Panel members as ACC=LABEL (sets speciesOrder/Labels and "
                        "chain targets). A bare ACC uses itself as the label.")
    p.add_argument("--anchor", "--anchors", nargs="+", default=[], metavar="ACC=LABEL",
                   help="Anchor strains contributing gene projections (annot sub-tracks).")
    p.add_argument("--maf", help="bigDataUrl for the standalone bigMaf track.")
    p.add_argument("--chains-dir", default="chains",
                   help="Sub-directory holding the bigChain .bb pairs (default: chains).")
    p.add_argument("--group", default="brc_pangenome",
                   help="trackDb group name (default: brc_pangenome).")
    p.add_argument("--is-reference", action="store_true",
                   help="Emit the reference-only selection/orthogroup composite.")
    p.add_argument("--maf-html", default="../shared/multiz.html")
    p.add_argument("--chains-html", default="../shared/chains.html")
    p.add_argument("--annot-html", default="../shared/annot.html")
    p.add_argument("--select-html", default="../shared/selection.html")
    p.add_argument("--output", "-o", required=True,
                   help="Output trackDb.txt path (or '-' for stdout).")
    args = p.parse_args(argv)
    build(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
