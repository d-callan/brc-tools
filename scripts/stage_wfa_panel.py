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

⚠ A BIG PANEL WANTS A DURABLE UPLOAD HISTORY, OR A FAILURE COSTS THE WHOLE TRANSFER AGAIN. By
default this creates a throwaway history and uploads everything into it, so a run that dies at file
150 of 192 starts over from zero on the next attempt -- it cannot even see what already landed.
Pass `--upload-history-file PATH` instead and the inputs accumulate in one history whose id is
recorded there: every later run skips what is already staged and `ok`, and only the three
collections are rebuilt. The Cannabis panel is 562 files and 65 GiB over the tunnel, where that is
the difference between a resume and a restart.

    python3 scripts/stage_wfa_panel.py --panel panels/cannabis.panel.yml \
        --upload-history-file ../data/wfa_upload_history.id

⚠ THE INPUTS THEN OUTLIVE THE RUN, AND RE-INVOKING IS NEARLY FREE. Staging copies the three
collections into a fresh run history BY REFERENCE -- the element HDAs differ but the underlying
`dataset_id` is the same, so a second invocation over an already-staged panel costs three API calls
and no bytes. `wfa_inputs*.json` names the RUN history and the copies, which is what the workflow
driver must use; `upload_history` in the same file names where the inputs actually live.

⛔ THE RESUME IS KEYED ON THE DATASET NAME, AND GALAXY DOES NOT ENFORCE UNIQUE NAMES. Two datasets
may share one name, so `resume_action` refuses an ambiguous name rather than picking one -- the same
rule `one_file` applies to the filesystem, one layer up. It also skips the HIDDEN copy that building
a collection leaves behind, without which nothing resumes at all; both facts were measured on vgp
rather than assumed, and `--self-test` pins every branch.

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


def api(path, payload=None, method=None):
    """A Galaxy API call. GET with no payload, POST with one, or `method` to override.

    ⚠ `method` EXISTS FOR `PUT`, which is how a dataset is purged. It defaults to the old
    behaviour so every existing call site reads the same.
    """
    r = urllib.request.Request(URL + path,
                               method=method or ("POST" if payload is not None else "GET"),
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


#: Dataset states in which the bytes are present and usable -- the only ones a resume may treat as
#: done. ⛔ `empty` IS DELIBERATELY ABSENT, for the reason `wait()` records at length: a zero-byte
#: upload used to pass staging and go into a collection.
DONE_STATES = frozenset({"ok"})

#: Terminal failures: the dataset exists but holds nothing usable. A resume PURGES these and
#: uploads again -- left in place they would occupy the name forever and stall every later run.
DEAD_STATES = frozenset({"error", "empty", "discarded", "failed_metadata"})


def resume_action(entries: list[dict]) -> tuple[str, str | None]:
    """What a resume must do about ONE upload name, given the datasets already carrying it.

    Returns `("reuse"|"replace", hda_id)`, `("upload", None)`, or `("refuse", reason)`.

    ⛔ KEPT PURE, AND SEPARATE FROM THE NETWORK, SO `--self-test` CAN REACH EVERY BRANCH. Each
    case below cost something to learn, and none of them is testable if the decision is written
    inline in `main()` beside the upload call.

    ⛔ GALAXY ACCEPTS DUPLICATE DATASET NAMES. Measured on vgp 2026-09-10: the same name uploaded
    twice yields two datasets, both `ok`, differing only in `hid`. So a name is NOT a key, and
    taking the first match is the very defect `one_file` exists to prevent, one layer up. Refuse.

    ⚠ A NON-TERMINAL STATE IS REFUSED RATHER THAN WAITED ON. `queued`/`running`/`upload` means
    something may be staging it RIGHT NOW -- plausibly another run of this script, since the whole
    point of a durable upload history is that runs share it. Purging it would corrupt that run and
    uploading alongside it would duplicate the name.

    ⚠ A DELETED OR PURGED DATASET IS TREATED AS ABSENT, not as a name in use. Uploading then leaves
    one deleted dataset and one live one under the same name, which the filter above resolves on
    the next pass -- so a purge-and-retry converges instead of deadlocking.

    ⛔ HIDDEN DATASETS ARE SKIPPED, AND WITHOUT THAT NOTHING RESUMES AT ALL. Building a `list`
    collection from HDAs does not reference them -- Galaxy makes a hidden (`visible: false`) copy
    of each element and leaves the original visible. Measured on vgp 2026-09-10: after one clean
    staging run the upload history held 8 datasets under 4 names, one visible and one hidden each,
    so the very next run saw every name as ambiguous and refused. The hidden twin shares the
    underlying dataset, so it costs nothing; it just must not be mistaken for a second upload.
    `visible` defaults to True when absent, because a plain upload listing may omit it.
    """
    live = [e for e in entries
            if not e.get("deleted") and not e.get("purged") and e.get("visible", True)]
    if not live:
        return "upload", None
    if len(live) > 1:
        return "refuse", (f"{len(live)} live datasets share this name ({[e['id'] for e in live]}). "
                          f"Which one the collection should hold is not something to guess at; "
                          f"purge the extras and re-run.")
    only = live[0]
    state = only.get("state")
    if state in DONE_STATES:
        return "reuse", only["id"]
    if state in DEAD_STATES:
        return "replace", only["id"]
    return "refuse", (f"dataset {only['id']} is in state {state!r}, which is not terminal -- "
                      f"another run may be staging it. Wait for it to settle, or purge it, rather "
                      f"than racing it.")


def history_index(history: str) -> dict[str, list[dict]]:
    """Every DATASET in a history, grouped by name, for `resume_action` to judge.

    ⚠ `history_content_type` IS FILTERED RATHER THAN ASSUMED. The listing returns collections too,
    and an HDCA carries a name just as a dataset does but cannot be reused as an upload -- so
    without this filter a collection named `assemblies` could be mistaken for a staged member.
    """
    listing = api(f"/api/histories/{history}/contents?v=dev"
                  f"&keys=name,id,state,deleted,purged,visible,history_content_type")
    out: dict[str, list[dict]] = {}
    for e in listing:
        if e.get("history_content_type") == "dataset":
            out.setdefault(e["name"], []).append(e)
    return out


def purge_dataset(history: str, hda: str) -> None:
    """Purge one dataset, freeing its name for a fresh upload."""
    api(f"/api/histories/{history}/contents/datasets/{hda}", {"purged": True}, method="PUT")


def copy_collection(dest_history: str, hdca: str) -> dict:
    """Copy a whole collection into another history BY REFERENCE, in one call.

    ⚠ THE BYTES DO NOT MOVE, AND THIS IS THE WHOLE REASON THE UPLOAD HISTORY PAYS OFF. Measured on
    vgp 2026-09-10: the copy's element HDA ids differ from the source's, but each element's
    underlying `dataset_id` is IDENTICAL -- one copy on disk, one charge against quota, however
    many run histories point at it. So re-invoking a workflow over an already-staged panel costs
    three API calls and no transfer.
    """
    return api(f"/api/histories/{dest_history}/contents",
               {"source": "hdca", "content": hdca, "type": "dataset_collection"})


def upload_history(panel: dict, explicit: str | None, id_file: pathlib.Path | None) -> str:
    """The durable history local files accumulate in, found or created.

    ⛔ THE ID IS PERSISTED RATHER THAN LOOKED UP BY NAME. Two histories may share a name -- Galaxy
    does not stop that any more than it stops duplicate dataset names -- so resolving by name
    would silently resume into whichever one the listing happened to return first, and a resume
    that picks the wrong history re-uploads everything while reporting progress.
    """
    if explicit:
        got = api(f"/api/histories/{explicit}")
        print(f"  upload history {explicit} ({got.get('name')!r}) -- reusing")
        return explicit
    if id_file and id_file.exists():
        hid = id_file.read_text(encoding="utf-8").strip()
        if hid:
            got = api(f"/api/histories/{hid}")
            # ⛔ A PURGED HISTORY IS NOT A RESUME TARGET. Its datasets are gone but the id still
            # resolves, so without this the run would "resume" into an empty history and upload
            # everything again -- correct in the end, but silently, and after hours.
            if got.get("purged") or got.get("deleted"):
                sys.exit(f"the upload history recorded in {id_file} ({hid}) is deleted or purged. "
                         f"Remove that file to start a fresh one; nothing here can tell whether "
                         f"you meant to re-upload {panel['expected_assemblies']} assemblies.")
            print(f"  upload history {hid} ({got.get('name')!r}) -- resuming from {id_file.name}")
            return hid
    hid = api("/api/histories",
              {"name": f"WF-A inputs — {panel['name']} "
                       f"({panel['expected_assemblies']} assemblies)"})["id"]
    if id_file:
        id_file.parent.mkdir(parents=True, exist_ok=True)
        id_file.write_text(hid + "\n", encoding="utf-8")
        print(f"  upload history {hid} -- created, recorded in {id_file}")
    else:
        print(f"  upload history {hid} -- created (pass --upload-history {hid} to resume)")
    return hid


def stage_local(history: str, index: dict[str, list[dict]], items: list[tuple[str, pathlib.Path]],
                suffix: str, ext: str, what: str) -> tuple[list[tuple[str, str]], list[str]]:
    """Upload `items` into `history`, skipping whatever is already there and usable.

    Returns `[(identifier, hda_id)]` in `items` order, and the ids of the ones actually uploaded
    (the only ones `wait()` then has to poll).

    ⚠ SERIAL, AND STILL DELIBERATELY SO. Resume changes what a failure COSTS, not what the tunnel
    can carry; uploading concurrently would not widen it, and would multiply what has to be redone.
    """
    out: list[tuple[str, str]] = []
    fresh: list[str] = []
    skipped = 0
    for n, (ident, path) in enumerate(items, 1):
        name = f"{ident}{suffix}"
        action, arg = resume_action(index.get(name, []))
        if action == "refuse":
            sys.exit(f"cannot resume {what} {ident} (dataset name {name!r}): {arg}")
        if action == "reuse":
            out.append((ident, arg))
            skipped += 1
            continue
        if action == "replace":
            print(f"    [{n}/{len(items)}] {what} {ident}: purging unusable {arg} and re-uploading",
                  flush=True)
            purge_dataset(history, arg)
        hda = upload(history, path, name, ext)
        out.append((ident, hda))
        fresh.append(hda)
        print(f"    [{n}/{len(items)}] uploaded {what} {ident} "
              f"({path.stat().st_size / 1048576:.0f} MB)", flush=True)
    print(f"  {what}: {len(out)} total -- {skipped} already staged, {len(fresh)} uploaded")
    return out, fresh


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

    # ── the resume decision table ──────────────────────────────────────────────────────────────
    # ⛔ EVERY BRANCH, INCLUDING THE REFUSALS, and each written by constructing the exact listing
    # Galaxy returns. `resume_action` decides whether hours of upload are repeated or skipped; a
    # test of only the happy path would let "reuse anything with the right name" pass.
    def act(entries):
        return resume_action(entries)[0]

    assert act([]) == "upload", "an unseen name must upload"
    assert act([{"id": "1", "state": "ok"}]) == "reuse", "a finished dataset must be reused"
    assert resume_action([{"id": "abc", "state": "ok"}])[1] == "abc", "reuse must return the id"

    # ⛔ MEASURED, NOT ASSUMED: Galaxy accepts two datasets under one name (vgp, 2026-09-10).
    # Taking the first would stage a collection against an arbitrary one of them.
    assert act([{"id": "1", "state": "ok"}, {"id": "2", "state": "ok"}]) == "refuse", \
        "two live datasets under one name must refuse, not pick one"

    # a terminal failure is redone; `empty` is a failure, per the reason wait() records
    for dead in ("error", "empty", "discarded", "failed_metadata"):
        assert act([{"id": "1", "state": dead}]) == "replace", f"{dead} must be replaced"

    # ⚠ IN-FLIGHT IS REFUSED, NOT WAITED ON: another run may own it, and a durable upload history
    # is precisely the thing two runs share.
    for busy in ("queued", "running", "upload", "new", "paused", "setting_metadata", "deferred"):
        assert act([{"id": "1", "state": busy}]) == "refuse", f"{busy} must refuse"

    # ⚠ DELETED/PURGED READS AS ABSENT, so purge-and-retry converges instead of deadlocking on a
    # name that can never be reused and can never be freed.
    assert act([{"id": "1", "state": "ok", "deleted": True}]) == "upload"
    assert act([{"id": "1", "state": "ok", "purged": True}]) == "upload"
    assert act([{"id": "1", "state": "ok", "deleted": True},
                {"id": "2", "state": "ok"}]) == "reuse", "a deleted twin must not make it ambiguous"

    # ⛔ THE HIDDEN COLLECTION-ELEMENT TWIN, which is what a `list` collection leaves behind. This
    # is not hypothetical: it broke the second run of the rehearsal on vgp before the filter
    # existed, turning every name into an ambiguity refusal.
    assert act([{"id": "1", "state": "ok", "visible": True},
                {"id": "2", "state": "ok", "visible": False}]) == "reuse", \
        "a hidden collection-element copy must not make its visible original ambiguous"
    assert resume_action([{"id": "vis", "state": "ok", "visible": True},
                          {"id": "hid", "state": "ok", "visible": False}])[1] == "vis", \
        "the VISIBLE dataset is the one to reuse"
    assert act([{"id": "1", "state": "ok", "visible": False}]) == "upload", \
        "a hidden dataset alone is not a staged input"

    # ⛔ THE THREE SUFFIXES MUST NOT COLLIDE, because all three types share one upload history and
    # the name is the resume key. Distinct suffixes are what keeps `X.faa.gz` from matching the
    # assembly `X.fasta.gz`.
    suffixes = (".fasta.gz", ".faa.gz", ".gff3")
    assert len(set(suffixes)) == len(suffixes), "the upload-name suffixes are not distinct"
    for a in suffixes:
        for b in suffixes:
            assert a == b or not a.endswith(b), f"{a} ends with {b}: names could be confused"

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
    ap.add_argument("--upload-history", metavar="ID",
                    help="reuse this history as the durable upload target, skipping anything "
                         "already staged in it. Without this (and --upload-history-file) the run "
                         "creates a throwaway history and uploads everything, as it always did.")
    ap.add_argument("--upload-history-file", type=pathlib.Path, metavar="PATH",
                    help="read the upload history id from PATH, creating the history and writing "
                         "the id there on the first run. The resumable form: one flag that is the "
                         "same on every run, with no id to copy by hand.")
    args = ap.parse_args()
    if args.upload_history and args.upload_history_file:
        ap.error("--upload-history and --upload-history-file both name the upload history; "
                 "pass one. The file form is the one to use in a script.")
    if args.self_test:
        return self_test()
    if args.panel is None:
        ap.error("--panel is required (or pass --self-test)")
    panel = load_panel(args.panel)
    print(f"  panel {args.panel.name}: {panel['name']}")

    # ⚠ NAMED FROM THE PANEL AND ITS ASSERTED COUNT, NOT A LITERAL. The label said "(4 genomes)"
    # for as long as the panel had six in it, and a history whose label disagrees with its contents
    # is read as the contents being wrong. The count is checked against what resolves, below.
    # ⚠ TWO HISTORIES, AND THE SPLIT IS THE POINT. `hist` is where the inputs LIVE: durable and
    # reused, so a failed run resumes instead of re-uploading. The run history is created after
    # staging succeeds and receives the three collections by reference, which costs no transfer --
    # so re-invoking a workflow over an already-staged panel is three API calls.
    #
    # ⚠ WITHOUT EITHER FLAG THE OLD BEHAVIOUR IS EXACT: one throwaway history, everything
    # uploaded, no index consulted. That matters because panels live outside this repository and
    # other people's scripts call this one.
    resumable = bool(args.upload_history or args.upload_history_file)
    if resumable:
        hist = upload_history(panel, args.upload_history, args.upload_history_file)
        index = history_index(hist)
        print(f"  {sum(len(v) for v in index.values())} dataset(s) already in it, "
              f"{len(index)} distinct name(s)")
    else:
        hist = api("/api/histories",
                   {"name": f"WF-A UDT — {panel['name']} "
                            f"({panel['expected_assemblies']} assemblies)"})["id"]
        index = {}
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
    copied = 0
    fresh_copy: list[str] = []
    for ident in wanted:
        if ident not in have:
            continue
        # ⛔ A COPY IS RENAMED TO THE SAME KEY AN UPLOAD WOULD USE, so all three assembly routes
        # are indexed alike and a resume can tell a copied member from an absent one. Without the
        # rename the copy keeps the SOURCE dataset's name, `resume_action` never matches it, and
        # every run copies it again -- duplicate names in a durable history, which then refuses.
        name = f"{ident}.fasta.gz"
        action, arg = resume_action(index.get(name, []))
        if action == "refuse":
            sys.exit(f"cannot resume copied assembly {ident} (dataset name {name!r}): {arg}")
        if action == "reuse":
            asm.append((ident, arg))
            continue
        if action == "replace":
            purge_dataset(hist, arg)
        d = api(f"/api/histories/{hist}/contents",
                {"source": "hda", "content": have[ident], "type": "dataset"})
        if resumable:
            api(f"/api/histories/{hist}/contents/datasets/{d['id']}", {"name": name}, method="PUT")
        asm.append((ident, d["id"]))
        fresh_copy.append(d["id"])
        copied += 1
    print(f"  copied {copied} assemblies from the staged panel "
          f"({len(asm) - copied} already present)")

    # ⚠ SERVER-SIDE FETCH, NOT AN UPLOAD. These two are not in the panel collection and their FASTAs
    # are not on this machine; Galaxy pulls them from NCBI directly, which costs no tunnel traffic.
    # ⚠ A FETCHED DATASET IS NAMED `ident`, WITH NO SUFFIX -- unlike an upload's `{ident}.fasta.gz`
    # -- and that difference is preserved rather than tidied away. The name is what an existing
    # history already holds, so changing it would make every panel staged before this change look
    # unstaged and re-fetch all of it. The two forms cannot collide: an identifier would have to
    # literally end in `.fasta.gz`, which the slug rule forbids.
    fetch, refetch = [], []
    for ident, spec in panel["by_url"].items():
        if ident not in wanted or ident in have:
            continue
        action, arg = resume_action(index.get(ident, []))
        if action == "refuse":
            sys.exit(f"cannot resume fetched assembly {ident} (dataset name {ident!r}): {arg}")
        if action == "reuse":
            asm.append((ident, arg))
            continue
        if action == "replace":
            purge_dataset(hist, arg)
        fetch.append({"src": "url", "url": ncbi_fasta(spec["accession"], spec["assembly_name"]),
                      "name": ident, "ext": "fasta.gz",
                      "to_posix_lines": False, "space_to_tab": False})
    if fetch:
        r = api("/api/tools/fetch", {"history_id": hist,
                                     "targets": [{"destination": {"type": "hdas"},
                                                  "elements": fetch}]})
        for o in r.get("outputs", []):
            asm.append((o["name"], o["id"]))
            refetch.append(o["id"])
        print(f"  fetching {len(fetch)} assembly FASTA(s) from NCBI")
    else:
        print(f"  0 assemblies to fetch from NCBI ({len(panel['by_url'])} already staged)")

    # ⚠ SERIAL, AND THAT IS A DELIBERATE CHOICE RATHER THAN AN OVERSIGHT. This is the only route
    # that spends tunnel bandwidth, and a panel of two hundred genomes is tens of GB; uploading
    # them concurrently does not make the tunnel wider, but it does multiply what has to be redone
    # when one of them fails. The progress line exists because a silent hour reads as a hang.
    disk_staged, fresh_disk = stage_local(hist, index, sorted(disk.items()),
                                          ".fasta.gz", "fasta.gz", "assembly")
    asm.extend(disk_staged)

    prot_files = [(ident, one_file(f"*_{key}_protein_primary.faa.gz", ident, "proteome"))
                  for ident, key in panel["proteomes"].items()]
    anch_files = [(ident, one_file(f"*_{key}_genomic.gff.gz", ident, "anchor GFF3"))
                  for ident, key in panel["anchors"].items()]
    prot, fresh_prot = stage_local(hist, index, prot_files, ".faa.gz", "fasta.gz", "proteome")
    anch, fresh_anch = stage_local(hist, index, anch_files, ".gff3", "gff3", "anchor GFF3")

    # ⛔ WAIT ONLY ON WHAT THIS RUN CREATED. Polling all 589 members costs one API call each per
    # cycle, and a reused dataset was already `ok` when the index was read -- `resume_action` will
    # not hand back anything else. In legacy mode every id is fresh, so this is the same set as
    # before.
    if resumable:
        settle = sorted(set(fresh_copy + refetch + fresh_disk + fresh_prot + fresh_anch))
        print(f"  waiting on {len(settle)} newly staged dataset(s) of "
              f"{len(asm) + len(prot) + len(anch)}")
    else:
        settle = [i for _, i in asm + prot + anch]
    wait(settle, "staging")

    # ⛔ THE COLLECTION MUST HOLD WHAT THE PANEL PROMISED, AND THIS IS CHECKED BEFORE ANY
    # COLLECTION EXISTS. Three routes feed one collection, and a member lost by any of them yields
    # a shorter panel that stages cleanly; the assertion in `wanted` above cannot see this, because
    # it counts the panel, not what reached Galaxy. It ran AFTER the build until 2026-09-10, which
    # left a wrong-length collection and a run history sitting in Galaxy to be cleaned up by hand.
    if len(asm) != panel["expected_assemblies"]:
        # ⚠ THE BREAKDOWN COUNTS WHAT THIS RUN DID, and says so, because with a resume most of a
        # correct panel is contributed by neither route -- it was already there. A line reading
        # "0 fetched, 1 uploaded" for a 219-member panel is not a contradiction.
        sys.exit(f"staged {len(asm)} assemblies but the panel expects "
                 f"{panel['expected_assemblies']} -- of which THIS RUN copied {copied}, fetched "
                 f"{len(fetch)} and uploaded {len(fresh_disk)}; the rest were already staged "
                 f"({len(have)} in the raw collection, {len(disk)} resolved on disk). The history "
                 f"is left in place to inspect.")

    c_asm = collection(hist, "assemblies", asm)
    c_prot = collection(hist, "proteomes", prot)
    c_anch = collection(hist, "anchor_gene_gff3s", anch)
    out = {"server": URL, "history": hist, "assemblies": c_asm["id"], "proteomes": c_prot["id"],
           "anchor_gene_gff3s": c_anch["id"]}

    # ⚠ THE INPUTS STAY PUT AND THE RUN GETS A COPY. A workflow writes its outputs into the
    # history it runs in, so invoking inside the upload history would mix results into the thing
    # being reused -- and a later `history_index` would then have to tell an input from an output.
    # The copy is by reference (see `copy_collection`), so this costs no transfer and no quota.
    if resumable:
        run_hist = api("/api/histories",
                       {"name": f"WF-A UDT — {panel['name']} "
                                f"({panel['expected_assemblies']} assemblies)"})["id"]
        out = {"server": URL, "history": run_hist, "upload_history": hist,
               "assemblies": copy_collection(run_hist, c_asm["id"])["id"],
               "proteomes": copy_collection(run_hist, c_prot["id"])["id"],
               "anchor_gene_gff3s": copy_collection(run_hist, c_anch["id"])["id"]}
        print(f"\n  run history {run_hist} -- three collections copied by reference")

    # ⚠ THE IDS PRINTED ARE THE ONES WRITTEN, taken from `out` rather than from the build. In
    # resumable mode those are the RUN history's copies, not the upload history's originals, and an
    # id from the wrong history is not an error that says so -- it is a 404 at invocation time.
    print(f"\n  assemblies        {out['assemblies']}  ({len(asm)})")
    print(f"  proteomes         {out['proteomes']}  ({len(prot)})")
    print(f"  anchor_gene_gff3s {out['anchor_gene_gff3s']}  ({len(anch)})")
    dest = OUT_DIR / f"wfa_inputs{os.environ.get('WFA_SERVER', '')}.json"
    dest.write_text(json.dumps(out, indent=1))
    print(f"  wrote {dest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
