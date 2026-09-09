#!/usr/bin/env python3
"""Stage WF-A's three input collections on a Galaxy server from a PANEL definition.

    python3 scripts/stage_wfa_panel.py --panel panels/cannabis.panel.yml

⛔ THE PANEL IS DATA, BECAUSE THE PIPELINE IS SHARED BETWEEN TWO ORGANISM SETS. These workflows run
over a Cannabis panel and a Plasmodium one from the same code. Which genomes are in a set, which of
them carry an annotation and which serve as anchors are facts about a dataset; they used to be
module constants here, which meant a second organism could only be supported by editing this file.
Add a file under panels/ instead.

⛔ THE THREE COLLECTIONS ARE INDEPENDENT, AND CONFLATING THEM IS WHAT HELD WF-A TO SIX GENOMES.
They were all built from one dict, so `assemblies` could never be larger than the set with protein
FASTAs -- while WF-B masked nineteen. Nothing in the workflow joins them: `assemblies` feeds
sourmash, faidx and the identifier list, `proteomes` feeds BUSCO alone, `anchor_gene_gff3s` feeds
anchor_prep alone. Three branches, three gates, and only BUSCO needs an annotation.

⚠ NOTHING DOWNSTREAM COMPARES THEM EITHER. `proteomes` and `assemblies` meet only at multiqc, which
merges a BUSCO summary keyed by one collection's identifiers with plots keyed by the other's and
compares nothing -- scripts/check_workflow_ports.py says so in as many words. So the containments
are asserted here, before anything is uploaded: every anchor must have a proteome, and every
proteome must name a member of the resolved assembly set.

⛔ PROTEOMES MUST BE ONE PROTEIN PER GENE, WHICH IS NOT WHAT NCBI SHIPS. RefSeq publishes every
isoform and BUSCO reads its `duplicated` fraction straight off the file, so a proteome carrying
isoforms is reported as duplicated for a reason that is not biology. That matters wherever a panel
holds haplotype pairs, where `duplicated` genuinely means an uncollapsed haplotype.
scripts/primary_proteome.py reduces a proteome to longest-per-gene; the files this expects are
named `*_<key>_protein_primary.faa.gz`, keyed by the panel's `proteomes` values.

⚠ AN ASSEMBLY REACHES GALAXY BY ONE OF THREE ROUTES, AND ONLY ONE OF THEM IS FREE. A member of
`raw_collection` is copied server-side; a `by_url` member is fetched by Galaxy straight from the
NCBI FTP, also server-side; a `from_disk` member is UPLOADED from this machine, which is the only
route that spends tunnel bandwidth. Prefer the first two -- but a genome published only as a
figshare tarball has no third-party URL Galaxy can use, because the deposited object is the
archive, not the FASTA inside it, and Galaxy has no tool to open one. Unpack it here and list it
under `from_disk`.

⚠ THE `_primary` IN THE NAME IS A CLAIM, NOT A CERTIFICATE -- this script cannot verify it, since
doing so needs the full proteome and its annotation, neither of which is staged. Audit the files
themselves before a panel's BUSCO numbers are compared across members:

    primary_proteome.py <full>.faa.gz --check <staged>_protein_primary.faa.gz \
        --gene-map-gff3 <same release>_genomic.gff.gz
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import sys
import time
import urllib.request

import yaml

#: Where the `_primary` proteome and anchor GFF3 files live. ⚠ NOT IN THIS REPOSITORY -- they are
#: several GB of downloads, so this points outside it and $WFA_PROTEOMES overrides.
PROTEOMES = pathlib.Path(os.environ.get("WFA_PROTEOMES", "proteomes")).expanduser()

#: Where the resulting collection ids are written, for the workflow driver to read.
OUT_DIR = pathlib.Path(os.environ.get("WFA_OUT_DIR", ".")).expanduser()

#: Keys a panel MAY define, and the value used when it does not. ⛔ SEPARATE FROM `PANEL_KEYS`
#: BECAUSE THIS REPOSITORY IS SHARED. Panels live beside their datasets, not here, so a key added
#: to the required set breaks every panel written before it existed -- including ones this checkout
#: has never seen. A new key is optional or it is a breaking change.
OPTIONAL_KEYS = {
    "from_disk": ({}, "identifier -> file key, for assemblies UPLOADED from local disk. For "
                      "anything with no URL Galaxy can fetch -- a genome deposited only as an "
                      "archive, or one that exists nowhere but this machine."),
}

#: Every key a panel file must define, and what it means. Validated before anything is staged,
#: because a panel missing a key would otherwise fail partway and leave a half-built history.
PANEL_KEYS = {
    "name": "human-readable label, used for the history name",
    "expected_assemblies": "how many assemblies the finished collection must hold (asserted)",
    "raw_collection": "id of a pre-staged collection to copy from, or null to fetch everything",
    "extra": "identifiers not in raw_collection, staged from `by_url`",
    "by_url": "identifier -> {accession, assembly_name} for anything fetched from NCBI",
    "proteomes": "identifier -> file key, for the members that HAVE a protein FASTA",
    "anchors": "identifier -> file key, for the members whose GFF3 seeds the projections",
}


def load_panel(path: pathlib.Path) -> dict:
    """Read and validate a panel definition.

    ⛔ THE PANEL IS DATA BECAUSE THE PIPELINE IS SHARED. These workflows run over two different
    organism sets -- a Cannabis panel and a Plasmodium one -- from the same code, so which genomes
    are in a panel, which of them carry an annotation and which serve as anchors are all facts
    about a dataset, not about the staging logic. They used to be module constants here, which
    meant a second organism could only be supported by editing this file.

    ⚠ VALIDATED UP FRONT, NOT AS IT GOES. Staging creates a history and uploads several GB before
    it would reach a missing key, and a half-built history is harder to reason about than a refusal.
    """
    try:
        doc = yaml.safe_load(path.read_text(encoding="utf-8"))
    except OSError as e:
        sys.exit(f"cannot read the panel {path}: {e}")
    if not isinstance(doc, dict):
        sys.exit(f"{path} does not parse as a mapping of panel keys")
    absent = [k for k in PANEL_KEYS if k not in doc]
    if absent:
        sys.exit(f"{path} is missing {absent}. Every key is required:\n" +
                 "\n".join(f"  {k}: {v}" for k, v in PANEL_KEYS.items()))
    # ⛔ A PANEL OF ZERO IS NOT A PANEL, AND `0 == 0` WOULD HAVE PASSED. The template ships
    # `expected_assemblies: 0`, and with no reachable collection the resolve check compared 0
    # against 0 and staged an empty history reporting success -- the same silent-zero this repo
    # keeps finding. An unfilled template must fail on the template, not on the run.
    if not isinstance(doc["expected_assemblies"], int) or doc["expected_assemblies"] < 1:
        sys.exit(f"{path}: `expected_assemblies` is {doc['expected_assemblies']!r}; it must be a "
                 f"positive integer. A copy of the template that was never filled in ends up here.")
    if str(doc["name"]).startswith("REPLACE ME"):
        sys.exit(f"{path}: `name` is still the template placeholder, so this panel was copied and "
                 f"not filled in. It names the Galaxy history, which is how a run is found later.")
    for k, (default, _) in OPTIONAL_KEYS.items():
        doc.setdefault(k, default if not isinstance(default, dict) else dict(default))
    for k in ("proteomes", "anchors", "by_url", *OPTIONAL_KEYS):
        if not isinstance(doc[k], dict):
            sys.exit(f"{path}: `{k}` must be a mapping, got {type(doc[k]).__name__}")
    if not isinstance(doc["extra"], list):
        sys.exit(f"{path}: `extra` must be a list, got {type(doc['extra']).__name__}")
    for ident, spec in doc["by_url"].items():
        if not (isinstance(spec, dict) and {"accession", "assembly_name"} <= set(spec)):
            sys.exit(f"{path}: by_url[{ident!r}] needs `accession` and `assembly_name` -- the FTP "
                     f"path is built from both, so a missing one yields a 404 at fetch time.")
    # ⛔ ONE SOURCE PER GENOME. An identifier in both `by_url` and `from_disk` would be staged
    # twice into the same collection under one name -- Galaxy accepts that, and the duplicate is
    # invisible in every downstream count because the collection reports the element once.
    both = sorted(set(doc["by_url"]) & set(doc["from_disk"]))
    if both:
        sys.exit(f"{path}: {both} appear in both `by_url` and `from_disk`. Which copy of the "
                 f"genome is intended is not something to guess at.")

    # ⛔ BOTH CONTAINMENTS, CHECKED AGAINST THE PANEL RATHER THAN THE SERVER, so a malformed panel
    # is caught without a network call. Nothing in WF-A can catch them: `proteomes` and
    # `assemblies` meet only at multiqc, which joins a BUSCO summary keyed by one collection to
    # plots keyed by the other and compares nothing.
    stray = sorted(set(doc["anchors"]) - set(doc["proteomes"]))
    if stray:
        sys.exit(f"{path}: anchors {stray} have no proteome. An anchor's GFF3 and its proteome come "
                 f"from the same annotation release, so one without the other is a staging error.")
    return doc


def creds() -> tuple[str, str]:
    """Server chosen by $WFA_SERVER: "" for usegalaxy.org, "_2" for vgp, "_3" for laila.

    ⚠ `_2` (vgp.usegalaxy.org) IS THE ONE TO REACH FOR ON A BIG PANEL. It has more resources
    dedicated to it than main, and it could not run user-defined tools at all until
    galaxyproject/usegalaxy-playbook#472 (deployed 2026-09-08) fixed two things: the app config
    lacked `enable_beta_tool_formats`, and its TPV had no destination accepting
    `tool_type_user_defined`, so an already-registered UDT was accepted at submit and then errored
    with exit_code None.

    ⚠ AND vgp SHARES MAIN'S DATABASE AND OBJECT STORE, so staging twice is waste: a history staged
    here against usegalaxy.org -- datasets and collections both -- is visible and usable from vgp by
    the same ids. Measured 2026-09-08.

    ⚠ ONE SCRIPT, TWO SERVERS, because the two staging runs must produce the SAME inputs. Forking it
    per server is how the collections quietly drift apart and a difference in the RESULT gets
    blamed on the workflow.
    """
    sfx = os.environ.get("WFA_SERVER", "")
    u = os.environ.get(f"GALAXY_URL{sfx}", "").rstrip("/")
    k = os.environ.get(f"GALAXY_API_KEY{sfx}", "")
    if not (u and k):
        sys.exit(f"set GALAXY_URL{sfx} and GALAXY_API_KEY{sfx}")
    print(f"  server {u}")
    return u, k


URL, KEY = creds()


def api(path, payload=None):
    r = urllib.request.Request(URL + path, method="POST" if payload is not None else "GET",
                               headers={"x-api-key": KEY, "content-type": "application/json"},
                               data=json.dumps(payload).encode() if payload is not None else None)
    with urllib.request.urlopen(r, timeout=900) as f:
        return json.load(f)


def upload(history: str, path: pathlib.Path, name: str, ext: str) -> str:
    """Upload a LOCAL file through /api/tools/fetch as multipart.

    ⛔ NOT bioblend's upload_file, AND NOT `upload1`. Both fail here, for different reasons worth
    recording:

      * bioblend uses the TUS protocol on any Galaxy >= 22.01 with no way to opt out, and the TUS
        endpoint Galaxy hands back is built from ITS OWN configured base URL -- `localhost:8080` on
        laila. Behind an SSH tunnel that address is the sandbox's own localhost, not Galaxy's, so
        every chunk goes nowhere. `TusUploadFailed ... Max retries exceeded`.
      * `upload1` is the classic upload tool and simply is not in laila's panel, which carries 63
        tools. `Tool not found`, 400014.

    /api/tools/fetch takes the file as a multipart part and the plan as a JSON string, and is
    present on both servers.

    ⚠ `trust_env=False` for the tunnel. host.docker.internal must NOT go through the sandbox's
    outbound proxy; usegalaxy.org must.
    """
    import requests
    local = "host.docker.internal" in URL or "localhost" in URL
    s = requests.Session()
    s.headers.update({"x-api-key": KEY})
    s.trust_env = not local
    targets = [{"destination": {"type": "hdas"},
                "elements": [{"src": "files", "name": name, "ext": ext,
                              "to_posix_lines": False, "space_to_tab": False}]}]
    last = None
    for attempt in range(5):
        try:
            with path.open("rb") as fh:
                r = s.post(f"{URL}/api/tools/fetch",
                           data={"history_id": history, "targets": json.dumps(targets)},
                           files={"files_0|file_data": (name, fh)}, timeout=3600)
            r.raise_for_status()
            outs = r.json().get("outputs") or []
            if not outs:
                raise RuntimeError(f"fetch returned no output: {r.text[:200]}")
            return outs[0]["id"]
        except Exception as exc:                                    # noqa: BLE001
            last = exc
            print(f"      upload {name} attempt {attempt+1} failed "
                  f"({type(exc).__name__}: {str(exc)[:90]}); retrying", flush=True)
            time.sleep(10 * (attempt + 1))
    raise SystemExit(f"upload {name}: gave up after 5 attempts -- {type(last).__name__}: {last}")


def ncbi_fasta(acc: str, name: str) -> str:
    prefix, num = acc.split("_")
    num = num.split(".")[0]
    safe = name.replace(" ", "_")
    return (f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{num[0:3]}/{num[3:6]}/{num[6:9]}/"
            f"{acc}_{safe}/{acc}_{safe}_genomic.fna.gz")


def wait(ids, label):
    while True:
        st = {i: api(f"/api/datasets/{i}")["state"] for i in ids}
        bad = [i for i, s in st.items() if s == "error"]
        if bad:
            sys.exit(f"{label}: error on {bad}")
        # ⛔ `empty` IS NOT DONE. This used to accept it alongside `ok`, so a zero-byte upload or a
        # server-side fetch that returned nothing passed staging and went into a collection. At 23
        # members that is 23 chances instead of 6, and a panel one genome short is exactly the kind
        # of thing every guard downstream assumes has already been ruled out.
        bad_empty = [i for i, s in st.items() if s == "empty"]
        if bad_empty:
            sys.exit(f"{label}: {len(bad_empty)} dataset(s) staged EMPTY: {bad_empty}. A zero-byte "
                     f"member is not a staged member -- delete the history and retry rather than "
                     f"building a collection around it.")
        left = [i for i, s in st.items() if s != "ok"]
        if not left:
            return
        print(f"    {label}: {len(ids)-len(left)}/{len(ids)} ready", flush=True)
        time.sleep(20)


def one_file(pattern, ident, what):
    """The single file matching `pattern`, or a refusal that names what was being staged.

    ⛔ THIS WAS `next(PROTEOMES.glob(...))`, WHICH THE DOCSTRING ABOVE ALREADY PROMISED IT WAS NOT.
    A missing file raised an uncaught StopIteration -- a traceback with no mention of the genome or
    the pattern -- and TWO matching files silently took whichever the filesystem yielded first,
    which is how a panel ends up staged against the wrong proteome with nothing to show it.
    """
    hits = sorted(PROTEOMES.glob(pattern))
    if not hits:
        sys.exit(f"no {what} for {ident}: nothing matches {pattern!r} under {PROTEOMES}. Set "
                 f"$WFA_PROTEOMES if the files live elsewhere; a partial panel is worse than none.")
    if len(hits) > 1:
        sys.exit(f"{len(hits)} files match {pattern!r} for {ident}'s {what}: "
                 f"{[h.name for h in hits]}. Which one is intended is not something to guess at.")
    return hits[0]


def collection(history, name, pairs):
    return api(f"/api/histories/{history}/contents", {
        "type": "dataset_collection", "collection_type": "list", "name": name,
        "element_identifiers": [{"src": "hda", "id": i, "name": n} for n, i in pairs]})


def self_test() -> int:
    """Exercise `load_panel` offline, breaking each property separately.

    ⛔ THE HAPPY PATH ALONE PROVES NOTHING. A panel that loads tells you the parser ran, not that
    any guard in it works; every check below is written by making a VALID panel invalid in exactly
    one way, so a deleted guard fails this rather than quietly widening what is accepted.
    """
    import tempfile

    base = {"name": "T", "expected_assemblies": 2, "raw_collection": None,
            "extra": ["a", "b"], "by_url": {"a": {"accession": "GCA_1", "assembly_name": "A"}},
            "proteomes": {}, "anchors": {}}

    def load(doc):
        with tempfile.NamedTemporaryFile("w", suffix=".yml", delete=False) as fh:
            yaml.safe_dump(doc, fh)
            return load_panel(pathlib.Path(fh.name))

    def refuses(doc, needle):
        try:
            load(doc)
        except SystemExit as e:
            assert needle in str(e), f"refused for the wrong reason: {e}"
            return
        raise AssertionError(f"accepted a panel it should refuse ({needle})")

    # ⚠ BACKWARD COMPATIBILITY IS THE FIRST CHECK, not an afterthought: a panel written before
    # `from_disk` existed must still load, or this change breaks every dataset in the wild.
    old_panel = load(dict(base))
    assert old_panel["from_disk"] == {}, "the optional key did not default"

    ok = load({**base, "from_disk": {"b": "bkey"}})
    assert ok["from_disk"] == {"b": "bkey"}

    refuses({**base, "from_disk": ["b"]}, "must be a mapping")
    refuses({**base, "from_disk": {"a": "akey"}}, "both `by_url` and `from_disk`")

    # the pre-existing guards must still bite after the edit
    refuses({**base, "anchors": {"a": "A"}}, "have no proteome")
    refuses({**base, "expected_assemblies": 0}, "positive integer")
    print("self-test: all guards bite")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--self-test", action="store_true",
                    help="validate panel parsing offline and exit; touches no server")
    ap.add_argument("--panel", type=pathlib.Path,
                    help="a panel definition under panels/ -- which genomes, which of them have a "
                         "proteome, and which are anchors. Required: this script stages whatever "
                         "panel it is given and knows nothing about any particular organism.")
    args = ap.parse_args()
    if args.self_test:
        return self_test()
    if args.panel is None:
        ap.error("--panel is required (or pass --self-test)")
    panel = load_panel(args.panel)
    print(f"  panel {args.panel.name}: {panel['name']}")

    # ⚠ NAMED FROM THE PANEL AND ITS ASSERTED COUNT, NOT A LITERAL. The label said "(4 genomes)"
    # for as long as the panel had six in it, and a history whose label disagrees with its contents
    # is read as the contents being wrong. The count is checked against what resolves, below.
    hist = api("/api/histories",
               {"name": f"WF-A UDT — {panel['name']} "
                        f"({panel['expected_assemblies']} assemblies)"})["id"]
    print(f"  history {hist}")

    have = {}
    if panel["raw_collection"]:
        try:
            raw = api(f"/api/dataset_collections/{panel['raw_collection']}?instance_type=history")
            have = {e["element_identifier"]: e["object"]["id"] for e in raw["elements"]}
        # ⚠ A DELIBERATE BROAD BOUNDARY: an unreachable pre-staged collection is not an error here,
        # it is a server that has to fetch. Whether that is survivable is decided by the count check
        # below -- a panel whose members mostly lack URLs will fail it, and should.
        except Exception as exc:                                             # noqa: BLE001
            print(f"  raw collection {panel['raw_collection']} not reachable here "
                  f"({type(exc).__name__}); staging from `by_url` alone")

    # ⛔ THE ASSEMBLY SET IS THE COLLECTION PLUS `extra`, DERIVED RATHER THAN LISTED. A second
    # hardcoded copy of the pre-staged identifiers would be a list that can disagree with the
    # collection it describes, and nothing would notice; the collection is the source of truth for
    # its own members.
    wanted = list(have) + [i for i in panel["extra"] if i not in have]
    wanted = list(dict.fromkeys(wanted + [i for i in panel["from_disk"] if i not in wanted]))
    missing = [k for k in wanted
               if k not in have and k not in panel["by_url"] and k not in panel["from_disk"]]
    if missing:
        sys.exit(f"no source for {missing}: not in the raw collection, and no `by_url` or "
                 f"`from_disk` entry")
    if len(wanted) != panel["expected_assemblies"]:
        sys.exit(f"the panel resolves to {len(wanted)} assemblies, expected "
                 f"{panel['expected_assemblies']} ({len(have)} copied + "
                 f"{len(panel['extra'])} extra). A silently shorter panel stages cleanly and every "
                 f"downstream count is then wrong, so this refuses rather than guessing which "
                 f"members went missing.")

    # ⛔ THE PROTEOME CONTAINMENT NEEDS THE RESOLVED SET, so unlike the anchor check in load_panel()
    # it can only run here. A proteome for a genome the panel does not contain would give BUSCO a
    # completeness number attributed to a strain that is not there, and nothing downstream compares
    # the two collections.
    stray_prot = sorted(set(panel["proteomes"]) - set(wanted))
    if stray_prot:
        sys.exit(f"{args.panel}: proteomes names {stray_prot}, which are not in the assembly set.")

    # ⛔ EVERY LOCAL FILE IS RESOLVED BEFORE THE FIRST BYTE GOES UP. `one_file` exits when a file
    # is absent or ambiguous, and discovering that partway through 192 uploads would leave a
    # history holding most of a panel -- which is harder to reason about than a refusal, and
    # tempting to invoke anyway. Fail on the filesystem, before the network.
    disk = {ident: one_file(f"*_{key}_assembly.fasta.gz", ident, "assembly")
            for ident, key in panel["from_disk"].items() if ident in wanted}

    asm = []
    for ident in wanted:
        if ident in have:
            d = api(f"/api/histories/{hist}/contents",
                    {"source": "hda", "content": have[ident], "type": "dataset"})
            asm.append((ident, d["id"]))
    print(f"  copied {len(asm)} assemblies from the staged panel")

    # ⚠ SERVER-SIDE FETCH, NOT AN UPLOAD. These two are not in the panel collection and their FASTAs
    # are not on this machine; Galaxy pulls them from NCBI directly, which costs no tunnel traffic.
    fetch = [{"src": "url", "url": ncbi_fasta(acc, nm), "name": ident, "ext": "fasta.gz",
              "to_posix_lines": False, "space_to_tab": False}
             for ident, (acc, nm) in ((k, (v["accession"], v["assembly_name"]))
                                      for k, v in panel["by_url"].items())
             if ident in wanted and ident not in have]
    if fetch:
        r = api("/api/tools/fetch", {"history_id": hist,
                                     "targets": [{"destination": {"type": "hdas"},
                                                  "elements": fetch}]})
        for o in r.get("outputs", []):
            asm.append((o["name"], o["id"]))
        print(f"  fetching {len(fetch)} assembly FASTA(s) from NCBI")

    # ⚠ SERIAL, AND THAT IS A DELIBERATE CHOICE RATHER THAN AN OVERSIGHT. This is the only route
    # that spends tunnel bandwidth, and a panel of two hundred genomes is tens of GB; uploading
    # them concurrently does not make the tunnel wider, but it does multiply what has to be redone
    # when one of them fails. The progress line exists because a silent hour reads as a hang.
    for n, (ident, f) in enumerate(sorted(disk.items()), 1):
        asm.append((ident, upload(hist, f, f"{ident}.fasta.gz", "fasta.gz")))
        print(f"    [{n}/{len(disk)}] uploaded assembly {ident} "
              f"({f.stat().st_size / 1048576:.0f} MB)", flush=True)

    prot, anch = [], []
    for ident, key in panel["proteomes"].items():
        f = one_file(f"*_{key}_protein_primary.faa.gz", ident, "proteome")
        prot.append((ident, upload(hist, f, f"{ident}.faa.gz", "fasta.gz")))
        print(f"    uploaded proteome {ident}")
    for ident, key in panel["anchors"].items():
        f = one_file(f"*_{key}_genomic.gff.gz", ident, "anchor GFF3")
        anch.append((ident, upload(hist, f, f"{ident}.gff3", "gff3")))
        print(f"    uploaded anchor   {ident}")

    wait([i for _, i in asm + prot + anch], "staging")

    c_asm = collection(hist, "assemblies", asm)
    c_prot = collection(hist, "proteomes", prot)
    c_anch = collection(hist, "anchor_gene_gff3s", anch)
    out = {"server": URL, "history": hist, "assemblies": c_asm["id"], "proteomes": c_prot["id"],
           "anchor_gene_gff3s": c_anch["id"]}
    # ⛔ THE COLLECTION MUST HOLD WHAT THE PANEL PROMISED. Three routes now feed one collection,
    # and a member lost by any of them yields a shorter panel that stages cleanly; the assertion
    # in `wanted` above cannot see this, because it counts the panel, not what reached Galaxy.
    if len(asm) != panel["expected_assemblies"]:
        sys.exit(f"staged {len(asm)} assemblies but the panel expects "
                 f"{panel['expected_assemblies']} -- {len(have)} copied, {len(fetch)} fetched, "
                 f"{len(disk)} uploaded. The history is left in place to inspect.")

    print(f"\n  assemblies        {c_asm['id']}  ({len(asm)})")
    print(f"  proteomes         {c_prot['id']}  ({len(prot)})")
    print(f"  anchor_gene_gff3s {c_anch['id']}  ({len(anch)})")
    dest = OUT_DIR / f"wfa_inputs{os.environ.get('WFA_SERVER', '')}.json"
    dest.write_text(json.dumps(out, indent=1))
    print(f"  wrote {dest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
