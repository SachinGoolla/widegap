"""Reproducibility manifest written next to every benchmark result."""

from __future__ import annotations

import datetime as _dt
import importlib.metadata as _md
import os
import platform
import subprocess
import sys
from pathlib import Path


def _git(*args: str, cwd: Path) -> str | None:
    try:
        return subprocess.run(
            ["git", *args], cwd=cwd, capture_output=True, text=True, check=True
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _cpu_model() -> str:
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or "unknown"


def build_manifest(repo: Path, command: list[str], seeds: list[int]) -> dict:
    packages = {}
    for name in ("widegap", "numpy", "scipy", "torch"):
        try:
            packages[name] = _md.version(name)
        except _md.PackageNotFoundError:
            packages[name] = None
    # Installed metadata goes stale in an editable install after a version bump; record the
    # imported source's own version as the authoritative one.
    from widegap import __version__ as _src_version

    packages["widegap_installed_metadata"] = packages["widegap"]
    packages["widegap"] = _src_version
    try:
        import torch

        torch_threads = torch.get_num_threads()
    except ImportError:
        torch_threads = None
    status = _git("status", "--porcelain", cwd=repo)
    return {
        "timestamp_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "command": command,
        "git_commit": _git("rev-parse", "HEAD", cwd=repo),
        "git_dirty": bool(status) if status is not None else None,
        "python": sys.version,
        "platform": platform.platform(),
        "cpu_model": _cpu_model(),
        "logical_cpus": os.cpu_count(),
        "torch_threads": torch_threads,
        "packages": packages,
        "seeds": seeds,
    }
