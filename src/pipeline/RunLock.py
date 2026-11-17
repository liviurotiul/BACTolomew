"""The run lock (`run.lock.json`): who is running a run, and whether that process still is.

A lock names the process that runs the pipeline (pid and host) and, since 0.3.1, enough to
tell that process apart from a later one that got the same pid: its start time, its command
line and a token of the process instance. A lock is stale when its process is gone on this
host or the pid now belongs to another process; a lock from another host is never judged.
No Qt here: the GUI and the pipeline both use it.
"""
from __future__ import annotations

import os
import socket
import subprocess
import sys
import time
import uuid
from datetime import datetime
from pathlib import Path
from typing import Optional

# One token per process: a lock carrying this process's pid but another token was written
# by an earlier process that had the same pid (an app restarted in a container, say).
INSTANCE_TOKEN = uuid.uuid4().hex

# Seconds two readings of one process's start time may differ by (`ps` has 1 s resolution).
_START_TOLERANCE_S = 2.0


def pid_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True  # it exists but belongs to another user
    except OSError:
        return False
    return True


def _ps_field(pid: int, field: str) -> Optional[str]:
    try:
        out = subprocess.run(
            ["ps", "-o", f"{field}=", "-p", str(pid)],
            capture_output=True, text=True, timeout=5, env={**os.environ, "LC_ALL": "C"},
        )
    except Exception:
        return None
    text = out.stdout.strip()
    return text if out.returncode == 0 and text else None


def process_start_time(pid: int) -> Optional[float]:
    """The process's start time (epoch seconds), or None when it cannot be read."""
    stat = Path(f"/proc/{pid}/stat")
    if stat.exists():
        try:
            fields = stat.read_text().rsplit(")", 1)[1].split()
            ticks = int(fields[19])  # field 22 of stat(5): start time in clock ticks after boot
            btime = next(
                int(line.split()[1]) for line in Path("/proc/stat").read_text().splitlines()
                if line.startswith("btime ")
            )
            return btime + ticks / os.sysconf("SC_CLK_TCK")
        except Exception:
            return None
    text = _ps_field(pid, "lstart")  # macOS and other systems without /proc
    if not text:
        return None
    try:
        return time.mktime(time.strptime(" ".join(text.split()), "%a %b %d %H:%M:%S %Y"))
    except Exception:
        return None


def process_command(pid: int) -> Optional[str]:
    """The process's command line, or None when it cannot be read."""
    cmdline = Path(f"/proc/{pid}/cmdline")
    if cmdline.exists():
        try:
            parts = [p for p in cmdline.read_bytes().decode("utf-8", "replace").split("\0") if p]
        except Exception:
            return None
        return " ".join(parts) or None
    return _ps_field(pid, "command")


def _iso_epoch(value) -> Optional[float]:
    """An ISO time (the lock's `start_time`) as epoch seconds, or None."""
    try:
        return datetime.fromisoformat(str(value)).timestamp()
    except (TypeError, ValueError):
        return None


def lock_payload(*, run_id: str, start_time: str) -> dict:
    """The lock this process writes when it starts a run."""
    pid = os.getpid()
    return {
        "run_id": run_id,
        "pid": pid,
        "host": socket.gethostname(),
        "start_time": start_time,
        "instance": INSTANCE_TOKEN,
        "process_start": process_start_time(pid),
        "process_command": process_command(pid) or " ".join([sys.executable, *sys.argv]),
    }


def stale_reason(data: dict) -> Optional[str]:
    """Why the lock `data` is stale, or None when its run may still be going (or the lock
    is from another host, which is not judged here)."""
    host = str(data.get("host") or "").strip()
    if host and host != socket.gethostname():
        return None
    try:
        pid = int(data.get("pid") or 0)
    except (TypeError, ValueError):
        pid = 0
    if not pid_alive(pid):
        return f"process {pid} is not running" if pid > 0 else "the lock names no process"
    if pid == os.getpid():
        # This process writes its token in every lock it takes.
        if data.get("instance") != INSTANCE_TOKEN:
            return f"process {pid} is this app, started after the run was interrupted"
        return None
    recorded_start = data.get("process_start")
    if isinstance(recorded_start, (int, float)):
        live_start = process_start_time(pid)
        if live_start is not None:
            if abs(live_start - float(recorded_start)) > _START_TOLERANCE_S:
                return f"process {pid} is another program started after the run"
            return None
    else:
        # A lock from before 0.3.1 has no process start: a process that started after the
        # run began cannot be the one running it.
        run_start = _iso_epoch(data.get("start_time"))
        live_start = process_start_time(pid) if run_start is not None else None
        if live_start is not None and live_start > run_start + _START_TOLERANCE_S:
            return f"process {pid} is another program started after the run"
    recorded_command = str(data.get("process_command") or "")
    if recorded_command:
        live_command = process_command(pid)
        if live_command is not None and live_command != recorded_command:
            return f"process {pid} is another program ({live_command[:80]})"
    return None
