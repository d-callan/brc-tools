"""Shared pieces of the softmask UDT scripts.

WHY THIS EXISTS. Four scripts drive the same workflow -- run_softmask_udt (the whole thing),
build_up_softmask (tier by tier), check_softmask_stages (one stage at a time) and
verify_softmask_outputs (the result). They had grown byte-identical copies of `connect`,
`await_dataset` and the invocation-state set, and near-identical copies of the UDT registration
and the FASTA counter.

⛔ THE DUPLICATED COPIES WERE NOT HARMLESS. Both hard-won corrections in this area had to be
applied three times by hand: `await_dataset` learning that `state in ("ok", "error")` is not a
readiness test, and SCHEDULING_IN_PROGRESS gaining `requires_materialization` and `cancelling`
from Galaxy's own enum. A fix that must be repeated N times is a fix that will eventually be
applied N-1 times.

Only things that were ALREADY the same live here. `main`, the renderers and the per-script
assertions differ for real reasons and stay where they are.
"""

from __future__ import annotations

import json
import pathlib
import sys
import time

import yaml
from bioblend.galaxy import GalaxyInstance

ROOT = pathlib.Path(__file__).resolve().parent.parent
UDT_DIR = ROOT / "udt"

# ⚠ scripts/ IS NOT A PACKAGE and these modules are run by path, so the sibling import resolves
# only once its directory is on sys.path -- `python3 scripts/x.py` puts it there, `python3 -m` and
# a symlinked entry point do not. Same insert, same reason, as check_workflow_ports.py.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import galaxy_server                             # noqa: E402, I001 -- must follow the path insert
from check_udt_definitions import pep440_ok      # noqa: E402 -- must follow the path insert
WORKFLOW = ROOT / "workflows/softmask/softmask_udt.gxwf.yml"

#: The UDTs the softmask workflow needs, in dependency order.
UDTS = ("fasta_uppercase", "dustmasker_bed3", "windowmasker_bed3", "tantan_bed3",
        "lc_classify", "samtools_faidx", "fastan_gdb", "fastan_scan", "fastan_bed",
        "masking_row", "masking_header")

#: Invocation states in which Galaxy may still create jobs. Everything else means scheduling is
#: finished, whatever the jobs are doing.
#:
#: ⛔ THE VALUES COME FROM GALAXY'S OWN ENUM, not from watching behaviour. lib/galaxy/schema/
#: invocation.py::InvocationState documents each one, and two of these were missing when this set
#: was written from observation alone:
#:     new                       "Brand new workflow invocation"
#:     ready                     "Workflow ready for another iteration of scheduling."
#:     requires_materialization  "an otherwise NEW or READY workflow that requires inputs to be
#:                                materialized (undeferred)"
#:     cancelling                "invocation scheduler will cancel job in next iteration."
#:
#: ⚠ AND THE SAME FILE SETTLES WHY `completed` CANNOT BE USED AS THE WAIT CONDITION. It defines
#: `scheduled` as "Workflow has been scheduled" and `completed` as "All jobs have reached terminal
#: states" -- so `completed` is the state one WANTS, and it is nevertheless unreliable: measured
#: over 60 invocations on this account, 50 `completed` and 9 `scheduled`, interleaved across the
#: whole timeline, with structurally identical runs landing differently and one sitting `scheduled`
#: for 6.8 days with all ten jobs `ok` and `update_time` frozen at creation. Galaxy records the
#: transition in a separate `workflow_invocation_completion` row (model/__init__.py), so an
#: invocation whose completion hook never fires stays `scheduled` forever. Wait on the JOBS.
SCHEDULING_IN_PROGRESS = ("new", "ready", "requires_materialization", "cancelling")

#: JOB states where the job has not finished and the model may still change -- Galaxy's own
#: `Job.non_ready_states` (model/__init__.py:1801-1808), copied rather than approximated.
#:
#: ⛔ THE APPROXIMATION WAS WRONG IN BOTH DIRECTIONS, and it was written out twice. Two drivers
#: carried `("new", "queued", "running", "paused")` inline:
#:
#:   * it MISSED `waiting`, `resubmitted` and `upload`, so a job in one of those made the pending
#:     set empty, the wait broke early, and the run was reported "⛔ jobs did not all succeed"
#:     while still progressing. `resubmitted` is the live one -- public instances resubmit a job
#:     that exceeded walltime to a larger destination, and WF-C's KegAlign step is exactly that
#:     kind of job.
#:   * it INCLUDED `paused`, which is neither terminal (`ok`/`error`/`deleted`) nor non-ready: a
#:     paused job waits for a person. Polling it to the ceiling turned "this run is paused" into
#:     "⛔ TIMED OUT", which reads like an infrastructure problem rather than a decision waiting to
#:     be made. Leaving `paused` out means the wait ends and the verdict reports it, since the
#:     verdict already fails anything that is not `ok`.
JOBS_UNFINISHED = ("new", "resubmitted", "upload", "waiting", "queued", "running")


def connect() -> GalaxyInstance:
    """A client for the server `$WFA_SERVER` selects. See `galaxy_server.creds`.

    ⛔ THIS READ PLAIN `GALAXY_URL` UNTIL 2026-09-10 AND SO COULD NOT SELECT A SERVER AT ALL, while
    `stage_wfa_panel` selected one by suffix. A 219-genome panel was staged on vgp and then invoked
    on main, because staging and invoking are different scripts and only one of them honoured the
    variable. Six scripts call this, so all six were affected; nothing raised, because vgp shares
    main's database and every id resolved on both.
    """
    return GalaxyInstance(*galaxy_server.creds())


def register_one(gi: GalaxyInstance, name: str) -> tuple[str, str, str]:
    """Create ONE UDT from udt/<name>.gxtool.yml, returning (tool_id, version, uuid)."""
    doc = yaml.safe_load((UDT_DIR / f"{name}.gxtool.yml").read_text(encoding="utf-8"))
    # ⛔ REFUSE A NON-PEP-440 VERSION HERE RATHER THAN LET THE SERVER DO IT. The create comes back
    # `400 Tool failed lint checks: ToolVersionPEP404`, which names a linter and not the field, and
    # arrives after the run has already set up a history. Worse, this helper registers a LIST of
    # tools: one bad version fails the batch partway, leaving the earlier tools registered and the
    # workflow un-runnable. Say which file and which value, before anything is created.
    if not pep440_ok(str(doc.get("version", ""))):
        sys.exit(f"{name}.gxtool.yml: version {doc.get('version')!r} is not PEP 440, and "
                 f"/api/unprivileged_tools refuses it (400 ToolVersionPEP404). Use a release "
                 f"(`0.2.0`), a dev release (`0.1.0.dev1`) or a local version (`0.1.0+probe1`).")
    created = gi.make_post_request(f"{gi.url}/unprivileged_tools",
                                   payload={"representation": doc}, params={"key": gi.key})
    return doc["id"], str(doc["version"]), created["uuid"]


def register_all(gi: GalaxyInstance, names: tuple[str, ...] = UDTS, verbose: bool = True) -> dict[str, str]:
    """Create each UDT, returning {tool_id: uuid}.

    ⚠ THERE IS NO UPDATE. Every create makes a NEW tool, so re-running this leaves the previous
    definition beside the new one. That is Galaxy's behaviour, not a bug here, but it means the
    uuid returned by THIS call is the only one guaranteed to match the YAML on disk -- never reuse
    a uuid recorded by an earlier run after editing a definition.
    """
    mapping = {}
    for name in names:
        tool_id, version, uuid = register_one(gi, name)
        mapping[tool_id] = uuid
        if verbose:
            print(f"  registered {tool_id:24} v{version} -> {uuid}")
    if not verbose:
        print(f"  registered {len(mapping)} UDT(s)")
    return mapping


def await_dataset(gi: GalaxyInstance, dataset_id: str, label: str, tries: int = 1200) -> None:
    """Block until a dataset is ready, and DIE if it errored or never settled.

    ⛔ `state in ("ok", "error")` IS NOT A READINESS TEST. It was used as one here, so an upload
    that failed -- bad format detection, truncated transfer -- became a collection element and the
    whole workflow ran on it. Exhausting the retries fell through just as silently.
    """
    state = "unknown"
    for _ in range(tries):
        state = gi.datasets.show_dataset(dataset_id)["state"]
        if state == "ok":
            return
        if state in ("error", "discarded", "failed_metadata"):
            info = gi.datasets.show_dataset(dataset_id).get("misc_info") or ""
            sys.exit(f"{label} is in state {state!r} and cannot be used: {info[:200]}")
        time.sleep(5)
    sys.exit(f"{label} never became ready after {tries * 5}s (last state {state!r}).")


def fasta_stats(text: str) -> tuple[int, int, int]:
    """(sequences, residues, lowercase residues).

    Callers that want only two of the three unpack and discard; this was two functions differing
    solely in whether they counted headers.
    """
    seqs = res = low = 0
    for line in text.splitlines():
        if line.startswith(">"):
            seqs += 1
            continue
        res += len(line)
        low += sum(1 for c in line if "a" <= c <= "z")
    return seqs, res, low


class UpgradeMessagesRefused(RuntimeError):
    """Galaxy refused an invocation because some step's `state:` leaves a parameter unset."""

    def __init__(self, data: dict) -> None:
        self.data = data
        lines = [f"step {i}: {k}: {v}"
                 for i in sorted(data, key=int) for k, v in sorted(data[i].items())]
        super().__init__("Galaxy refused the invocation over unset parameters:\n  "
                         + "\n  ".join(lines))


def invoke(gi: GalaxyInstance, wf_id: str, inputs: dict, history_id: str,
           use_cached_job: bool = False) -> dict:
    """Invoke a workflow WITHOUT allow_tool_state_corrections, reporting what it would have hidden.

    ⛔ THE FLAG WAS NEVER A FIX. `workflow/modules.py::populate_module_and_state` either raises on a
    step's upgrade messages or, with the flag, calls `log.debug` -- to Galaxy's server log, which no
    response exposes. Passing it does not settle which value a parameter takes; it only removes the
    one place that would have told us the question was open. Every parameter this workflow's steps
    can take is now named in softmask_udt.gxwf.yml, so there is nothing to silence, and a refusal
    here is real news rather than noise to be switched off.

    ⚠ THE REFUSAL IS THE ONLY RELIABLE AUDIT. Comparing a step's `state:` against the tool's
    parameter list misses two shapes, both of which really occurred here: a parameter nested inside
    a repeat's conditional, and an OPTIONAL `data` input, which Galaxy counts as unset exactly like
    a required one -- leaving it unconnected is not the same as naming it null. Galaxy also raises
    on the FIRST offending step only, so one refusal is a floor, not a census.

    ⚠ `use_cached_job` MAKES A RE-INVOKE CHEAP, AND ITS PRECONDITION IS A TRAP. Galaxy reuses any
    prior job with identical tool version, inputs and parameters instead of re-running it, which
    turns a re-invoke over an already-staged panel into minutes instead of tens of core-hours. But
    the cache judges a job by its STATE, not by its outputs -- and a job killed at the scheduler
    level can land `state=ok` with `exit_code=None` and ZERO-BYTE outputs. Measured on vgp
    2026-09-10: `scancel` on a wedged busco left the job green with an empty `failed_metadata`
    dataset, and sourmash green with `max_containment` at 0 bytes. Re-invoking with the cache on
    would have reused both of those broken jobs and reported success.
    ⛔ SO BEFORE RE-INVOKING WITH THE CACHE, DELETE EVERY OUTPUT OF THE JOBS BEING REDONE -- all of
    them, not just the ones the workflow wires into collections. An undeleted output is a cache hit.
    """
    try:
        return gi.workflows.invoke_workflow(wf_id, inputs=inputs, history_id=history_id,
                                            use_cached_job=use_cached_job)
    except Exception as exc:
        body = getattr(exc, "body", None)
        if isinstance(body, str):
            try:
                body = json.loads(body)
            except ValueError:
                body = None
        data = body.get("err_data") if isinstance(body, dict) else None
        if data:
            raise UpgradeMessagesRefused(data) from exc
        raise
