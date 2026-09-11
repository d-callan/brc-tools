#!/usr/bin/env python3
"""Re-run WF-A's panel_relabel_map with brc-relabel-map 0.2.0 and check what came back.

WHY THIS EXISTS AS A RUN AND NOT A LOCAL CHECK: the 36-row map WF-A published on vgp is what a
WF-C invocation would pick up, and `strict: true` now refuses it on the row count. The fix is one
job, not a re-run of WF-A -- but the OUTPUT is the thing to verify, so this asserts on it.

⚠ THE DIAGONAL TEST CANNOT SPLIT ON THE DOT. Two panel members carry a `.` in their own name
(T2T_GCA_054642775.1_IMPD), so `A.B`.split(".") finds the wrong boundary and under-counts. The
identifier list is known, so the pairs are reconstructed from it instead of parsed out.

    python3 scripts/regen_relabel_map.py <tool_uuid> [history_id] [identifiers_dataset_id]

Register the tool first (`scripts/udt_deploy.py udt/relabel_map.gxtool.yml --instance vgp
--register-only`) and pass the uuid it prints -- a UDT is invoked by uuid, never by tool id.
"""
import os
import sys
import time

from bioblend.galaxy import GalaxyInstance

URL, KEY = os.environ["GALAXY_URL_2"].rstrip("/"), os.environ["GALAXY_API_KEY_2"]
gi = GalaxyInstance(url=URL, key=KEY)
#: The six-genome WF-A UDT run on vgp, and the `collection_element_identifiers` output it used.
#: Overridable, because the same repair applies to any panel; these are the defaults because they
#: are the run whose map was superseded on 2026-09-08.
UUID = sys.argv[1]
HIST = sys.argv[2] if len(sys.argv) > 2 else "bbd44e69cb8906b5895c6e9f6256bf78"
IDENTIFIERS = sys.argv[3] if len(sys.argv) > 3 else "f9cad7b01a4721354330db38673bdc85"

ids = [x for x in gi.datasets.download_dataset(IDENTIFIERS, use_default_filename=False)
       .decode().splitlines() if x.strip()]
n = len(ids)
print(f"{n} identifiers from the same dataset WF-A used", flush=True)

r = gi.make_post_request(URL + "/api/tools", params={"key": KEY}, payload={
    "history_id": HIST, "tool_uuid": UUID,
    "inputs": {"identifiers": {"src": "hda", "id": IDENTIFIERS}}})
job, out = r["jobs"][0]["id"], r["outputs"][0]["id"]
print(f"submitted job {job}", flush=True)

for _ in range(240):
    st = gi.jobs.show_job(job)["state"]
    if st in ("ok", "error", "paused", "deleted"):
        break
    time.sleep(5)
print(f"job {st}", flush=True)
if st != "ok":
    j = gi.jobs.show_job(job, full_details=True)
    print("  stderr:", (j.get("stderr") or "(empty)")[:500])
    sys.exit(1)

rows = [x for x in gi.datasets.download_dataset(out, use_default_filename=False)
        .decode().splitlines() if x.strip()]
want_pairs = {f"{a}_{b}": f"{a}.{b}" for a in ids for b in ids if a != b}
diagonal = {f"{a}_{a}" for a in ids}

got = {}
bad_cols = 0
for line in rows:
    c = line.split("\t")
    if len(c) != 2:
        bad_cols += 1
        continue
    got[c[0]] = c[1]

print(f"  rows            {len(rows)}   (n**2-n = {n*n-n}; the old map had n**2 = {n*n})")
print(f"  malformed rows  {bad_cols}")
print(f"  diagonal rows   {len(diagonal & set(got))}   (must be 0)")
print(f"  exactly the expected off-diagonal pairs: {got == want_pairs}")

ok = (len(rows) == n * n - n and bad_cols == 0 and not (diagonal & set(got))
      and got == want_pairs)
# The point of the job: this file must now pass the strict row-count check.
print(f"\n  {'PASS' if ok else '⛔ FAIL'} -- new map dataset {out}")
print(f"  history {URL}/histories/view?id={HIST}")
sys.exit(0 if ok else 1)
