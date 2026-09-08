#!/usr/bin/env python3
"""Register a User-Defined Tool on a Galaxy server and smoke-test it.

WHY A SCRIPT AND NOT A SNIPPET. Porting a wrapper is a LOOP -- convert, register, run, read the
error, fix, repeat -- and each turn of it needs the same four API calls. Retyping them invites the
small differences that make one attempt not comparable with the last.

⚠ FOUR API DETAILS THAT ARE NOT IN THE OBVIOUS PLACE, each of which cost a failed attempt:

  * A UDT is INVOKED BY `tool_uuid`, not by `tool_id`. Posting the id to /api/tools returns
    "Tool not found", which reads like the registration failed when it did not.
  * There is NO UPDATE. Every `create` makes a new tool, so an unchanged `version` leaves the old
    definition sitting beside the new one and `--bump` exists to keep them distinguishable.
  * THE VERSION MUST BE PEP 440, INCLUDING WHATEVER `--bump` INVENTS. `0.1.0-probe` -- the obvious
    thing to type for a throwaway -- is refused with `400 Tool failed lint checks:
    ToolVersionPEP404`, and Galaxy's own linter calls that a WARNING, so nothing in the message
    points at the version. Checked here before the request; see pep440_ok().
  * A tool with NO data inputs fails before reaching a node -- no runner, no stderr, error in about
    a second. Give the smoke test a real input.

    python3 scripts/udt_deploy.py udt/dustmasker.yml --smoke-fasta test.fa
    python3 scripts/udt_deploy.py udt/dustmasker.yml --register-only
"""
from __future__ import annotations

import argparse
import os
import pathlib
import sys
import time

import yaml
from bioblend.galaxy import GalaxyInstance

# ⚠ scripts/ IS NOT A PACKAGE and this is run by path, so the sibling import resolves only once
# its directory is on sys.path. Same insert, same reason, as check_workflow_ports.py.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from check_udt_definitions import pep440_ok      # noqa: I001 -- must follow the path insert


def creds(instance: str) -> tuple[str, str]:
    if instance == "main":
        return os.environ["GALAXY_URL"], os.environ["GALAXY_API_KEY"]
    return os.environ["GALAXY_URL_2"], os.environ["GALAXY_API_KEY_2"]


def wait(gi: GalaxyInstance, job_id: str, timeout: int = 1800) -> dict:
    t0 = time.time()
    while time.time() - t0 < timeout:
        j = gi.jobs.show_job(job_id, full_details=True)
        if j.get("state") in ("ok", "error", "deleted", "paused"):
            return j
        time.sleep(10)
    return gi.jobs.show_job(job_id, full_details=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("definition", type=pathlib.Path)
    ap.add_argument("--instance", default="main")
    ap.add_argument("--register-only", action="store_true")
    ap.add_argument("--smoke-fasta", type=pathlib.Path,
                    help="a small FASTA to run the tool on. ⚠ Without an input the job dies before "
                         "reaching a node, with no stderr to read.")
    ap.add_argument("--bump", help="override the version, since create never updates in place")
    ap.add_argument("--input-name", default="input", help="the tool's data input parameter")
    args = ap.parse_args()

    doc = yaml.safe_load(args.definition.read_text())
    if args.bump:
        doc["version"] = args.bump
    # ⛔ CHECK THE VERSION BEFORE THE NETWORK, AND CHECK IT AFTER `--bump` HAS BEEN APPLIED.
    # `--bump` takes a free-form string, and the obvious thing to type for a throwaway is exactly
    # what the server refuses: `0.1.0-probe`. The refusal is
    # `400 Tool failed lint checks: ToolVersionPEP404`, which names a linter rather than the field,
    # and Galaxy's own linter records this as a WARNING -- so nothing about the message suggests
    # the fix is the version you just invented. Measured on usegalaxy.org 26.1, 2026-09-08.
    if not pep440_ok(str(doc.get("version", ""))):
        src = "--bump" if args.bump else f"{args.definition.name}'s `version`"
        sys.exit(f"{src} is {doc.get('version')!r}, which is not PEP 440, and "
                 f"/api/unprivileged_tools refuses it (400 ToolVersionPEP404). PEP 440 has no bare "
                 f"alphabetic suffix: use a release (`0.2.0`), a dev release (`0.1.0.dev1`) or a "
                 f"local version (`0.1.0+probe1`). `0.1.0-1` also works -- `-<digits>` is an "
                 f"implicit post-release -- but `-probe` is not a version at all.")
    url, key = creds(args.instance)
    gi = GalaxyInstance(url=url, key=key)

    r = gi.make_post_request(url.rstrip("/") + "/api/unprivileged_tools",
                             payload={"representation": doc}, params={"key": key})
    uuid = r["uuid"]
    print(f"registered {doc['id']} v{doc['version']} -> {uuid}")
    if args.register_only or not args.smoke_fasta:
        return 0

    h = gi.histories.create_history(name=f"UDT smoke: {doc['id']} {doc['version']}")
    up = gi.tools.paste_content(args.smoke_fasta.read_text(), h["id"], file_type="fasta")
    ds = up["outputs"][0]["id"]
    for _ in range(60):
        if gi.datasets.show_dataset(ds)["state"] == "ok":
            break
        time.sleep(3)

    j = gi.make_post_request(
        url.rstrip("/") + "/api/tools",
        payload={"history_id": h["id"], "tool_uuid": uuid,          # uuid, NOT tool_id
                 "inputs": {args.input_name: {"src": "hda", "id": ds}}}, params={"key": key})
    job = wait(gi, j["jobs"][0]["id"])
    print(f"  job {job['state']}  exit={job.get('exit_code')}")
    print(f"  history: {url.rstrip('/')}/histories/view?id={h['id']}")
    if job["state"] != "ok":
        # ⚠ Print BOTH streams: a container that fails to pull leaves tool_stderr empty and puts
        # nothing useful in stdout either, which is itself the diagnosis.
        print("  stderr:", (job.get("stderr") or "(empty)")[:600])
        print("  stdout:", (job.get("stdout") or "(empty)")[:300])
        return 1
    out = next(iter(job["outputs"].values()))
    d = gi.datasets.show_dataset(out["id"])
    print(f"  output: {d.get('extension')} {d.get('file_size')} bytes")
    body = gi.datasets.download_dataset(d["id"], use_default_filename=False)
    text = body.decode("utf-8", "replace") if isinstance(body, bytes) else str(body)
    print("  first lines:")
    for line in text.splitlines()[:5]:
        print("   ", line[:100])
    if not text.strip():
        print("  ⛔ EMPTY OUTPUT — the job succeeded and produced nothing, which is a silent failure")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
