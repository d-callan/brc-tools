#!/usr/bin/env python3
"""Which Galaxy server a script talks to, in ONE place.

    from galaxy_server import creds
    url, key = creds()          # honours $WFA_SERVER

⛔ THIS EXISTS BECAUSE THREE COPIES DISAGREED, AND ONE OF THEM COULD NOT SELECT AT ALL.
On 2026-09-10 a 219-genome panel was staged on vgp (`WFA_SERVER=_2`) and then INVOKED on
usegalaxy.org, because `stage_wfa_panel.creds()` read the suffixed variables while
`softmask_lib.connect()` read plain `GALAXY_URL`/`GALAXY_API_KEY` and ignored `$WFA_SERVER`
entirely. `stage_wfa_panel`'s own docstring claimed "both halves read the same variable so they
cannot disagree about which server the ids belong to" -- they did not read the same variable, and
they disagreed for a whole run.

⚠ THE MISROUTING DOES NOT RAISE, WHICH IS WHY IT SURVIVED. vgp shares main's database and object
store, so a history staged against one is visible from the other by the same ids and every API call
returns 200 either way. The only symptom is that the work lands in the wrong queue -- exactly the
resource decision the server choice was made for. A third server (`_3`, laila) does NOT share that
database, so there the same mistake yields a 404 that at least says something is wrong.

⚠ DEPENDENCY-FREE ON PURPOSE. `stage_wfa_panel` avoids bioblend deliberately (see its `upload`
docstring: bioblend forces TUS, whose endpoint is built from Galaxy's own configured base URL and
so points at the wrong host behind a tunnel). Importing the shared selector must not drag bioblend
back in, so this module imports nothing but the standard library.
"""
from __future__ import annotations

import os
import sys

#: $WFA_SERVER value -> what that server is, for error messages. Not a whitelist: an unknown suffix
#: is allowed as long as its variables exist, so a fourth server needs no edit here.
KNOWN = {
    "": "usegalaxy.org (main)",
    "_2": "vgp.usegalaxy.org -- more dedicated resources; runs UDTs since usegalaxy-playbook#472",
    "_3": "laila (does NOT share main's database)",
}


def suffix() -> str:
    """The configured server suffix, `""` meaning main."""
    return os.environ.get("WFA_SERVER", "")


def creds(sfx: str | None = None, *, announce: bool = True) -> tuple[str, str]:
    """`(url, key)` for the selected server, refusing rather than falling back.

    ⛔ NO FALLBACK TO THE UNSUFFIXED VARIABLES. Reading `GALAXY_URL` when `GALAXY_URL_2` is unset
    is precisely the bug this module was written for: it silently sends the work to main while
    every log line still says `_2`. A missing variable is a refusal.
    """
    sfx = suffix() if sfx is None else sfx
    url = os.environ.get(f"GALAXY_URL{sfx}", "").rstrip("/")
    key = os.environ.get(f"GALAXY_API_KEY{sfx}", "")
    if not (url and key):
        sys.exit(f"set GALAXY_URL{sfx} and GALAXY_API_KEY{sfx} "
                 f"(WFA_SERVER={sfx!r} selects {KNOWN.get(sfx, 'an unlisted server')})")
    if announce:
        # ⚠ PRINTED EVERY TIME, because "which server" is the fact that was wrong for a whole run
        # and invisible in the output. Cheap line, and it is the one a reader checks first.
        print(f"  server {url}  (WFA_SERVER={sfx!r})")
    return url, key


def self_test() -> int:
    """Prove the selector picks by suffix and refuses a half-configured server."""
    saved = {k: os.environ.get(k) for k in
             ("WFA_SERVER", "GALAXY_URL", "GALAXY_API_KEY", "GALAXY_URL_2", "GALAXY_API_KEY_2")}
    try:
        os.environ.update({"GALAXY_URL": "https://main.example", "GALAXY_API_KEY": "kmain",
                           "GALAXY_URL_2": "https://vgp.example", "GALAXY_API_KEY_2": "kvgp"})

        os.environ["WFA_SERVER"] = ""
        assert creds(announce=False) == ("https://main.example", "kmain")

        os.environ["WFA_SERVER"] = "_2"
        assert creds(announce=False) == ("https://vgp.example", "kvgp"), "suffix not honoured"

        # ⛔ THE REGRESSION ITSELF: selecting _2 must NEVER yield main's url.
        os.environ.pop("GALAXY_URL_2")
        try:
            creds(announce=False)
        except SystemExit as e:
            assert "GALAXY_URL_2" in str(e), f"refused for the wrong reason: {e}"
        else:
            raise AssertionError("fell back to the unsuffixed variables instead of refusing")

        # a trailing slash must not survive into the url, or every path becomes a double slash
        os.environ.update({"GALAXY_URL_2": "https://vgp.example/", "WFA_SERVER": "_2"})
        assert creds(announce=False)[0] == "https://vgp.example"
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
    print("galaxy_server self-test: selector honours the suffix and refuses a partial config")
    return 0


if __name__ == "__main__":
    sys.exit(self_test())
