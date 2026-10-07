"""sha256 manifest and run log: every input file and every stage run is recorded."""
import hashlib
import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

from bte import config as C

MANIFEST = C.RESULTS / "MANIFEST.sha256"
RUN_LOG = C.RESULTS / "run_log.json"


def sha256(path: Path, chunk=1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while b := f.read(chunk):
            h.update(b)
    return h.hexdigest()


def git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=C.ROOT, text=True).strip()
    except Exception:
        return "unknown"


def record_file(local: Path, url: str, release: str) -> str:
    """Append (or refresh) one line in the manifest: sha256  path  url  release."""
    digest = sha256(local)
    rel = local.relative_to(C.ROOT)
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    lines = [l for l in (MANIFEST.read_text().splitlines() if MANIFEST.exists() else [])
             if not l.split("  ")[1:2] == [str(rel)]]
    lines.append(f"{digest}  {rel}  {url}  {release}")
    MANIFEST.write_text("\n".join(sorted(lines, key=lambda l: l.split("  ")[1])) + "\n")
    return digest


def log_run(stage: str, params: dict, outputs: list[str], seconds: float):
    RUN_LOG.parent.mkdir(parents=True, exist_ok=True)
    log = json.loads(RUN_LOG.read_text()) if RUN_LOG.exists() else {"runs": []}
    log["runs"] = [r for r in log["runs"] if r["stage"] != stage]
    log["runs"].append({"stage": stage, "finished": datetime.now(UTC).isoformat(timespec="seconds"),
                        "commit": git_commit(), "seconds": round(seconds, 1), "params": params,
                        "outputs": outputs})
    RUN_LOG.write_text(json.dumps(log, indent=1) + "\n")
