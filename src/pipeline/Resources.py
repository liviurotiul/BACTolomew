"""The CPU and memory this process may use (Qt free).

On a workstation that is the machine; on a cluster node or in a container it is less:
the CPU affinity mask (taskset, SLURM's task binding), the cgroup's CPU quota and memory
limit (cgroup v2 `cpu.max` / `memory.max`, v1 `cpu.cfs_quota_us` / `memory.limit_in_bytes`)
and the SLURM allocation (`SLURM_CPUS_PER_TASK`, `SLURM_MEM_PER_NODE`, `SLURM_MEM_PER_CPU`).
Each limit found caps the budget; the smallest one wins and is named in `describe()`.
"""
from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping, Optional

_CGROUP_ROOT = Path("/sys/fs/cgroup")
_UNLIMITED_BYTES = 1 << 60  # cgroup v1 writes a page-rounded LONG_MAX for "no limit"


@dataclass(frozen=True)
class ResourceBudget:
    cores: int
    memory_gb: Optional[float]                      # None: unknown
    core_limits: dict = field(default_factory=dict)  # source -> cores
    memory_limits: dict = field(default_factory=dict)  # source -> GB

    def describe(self) -> str:
        def _src(limits: dict, value) -> str:
            hits = [k for k, v in limits.items() if v == value]
            return hits[0] if hits else "machine"

        cores = f"{self.cores} cores ({_src(self.core_limits, self.cores)})"
        if self.memory_gb is None:
            return f"{cores}, memory unknown"
        return f"{cores}, {self.memory_gb:.1f} GB RAM ({_src(self.memory_limits, self.memory_gb)})"


def _read(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8").strip()
    except OSError:
        return ""


def _cgroup_dirs(proc_cgroup: str, controller: str, root: Path) -> list[Path]:
    """The cgroup folders that can hold `controller`'s limits for this process, own folder
    first then each parent (a parent's limit applies too); v2 when `controller` is ""."""
    rel = None
    for line in proc_cgroup.splitlines():
        parts = line.split(":", 2)
        if len(parts) != 3:
            continue
        if controller == "" and parts[0] == "0" and parts[1] == "":
            rel, base = parts[2], root
        elif controller and controller in parts[1].split(","):
            rel, base = parts[2], root / parts[1]
            if not base.is_dir():
                base = root / controller
        else:
            continue
        break
    if rel is None:
        return []
    out: list[Path] = []
    here = base / rel.lstrip("/")
    while True:
        if here.is_dir():
            out.append(here)
        if here == base or base not in here.parents:
            break
        here = here.parent
    if base not in out and base.is_dir():
        out.append(base)  # inside a container the own path is the mount itself
    return out


def cgroup_cpu_limit(root: Path = _CGROUP_ROOT, proc_cgroup: Optional[str] = None) -> Optional[float]:
    text = proc_cgroup if proc_cgroup is not None else _read(Path("/proc/self/cgroup"))
    best: Optional[float] = None
    for folder in _cgroup_dirs(text, "", root):
        fields = _read(folder / "cpu.max").split()
        if len(fields) == 2 and fields[0] != "max":
            try:
                quota = int(fields[0]) / int(fields[1])
            except (ValueError, ZeroDivisionError):
                continue
            best = quota if best is None else min(best, quota)
    for folder in _cgroup_dirs(text, "cpu", root):
        try:
            quota, period = int(_read(folder / "cpu.cfs_quota_us")), int(_read(folder / "cpu.cfs_period_us"))
        except ValueError:
            continue
        if quota > 0 and period > 0:
            best = quota / period if best is None else min(best, quota / period)
    return best


def cgroup_memory_limit(root: Path = _CGROUP_ROOT, proc_cgroup: Optional[str] = None) -> Optional[int]:
    """The tightest cgroup memory limit in bytes, or None."""
    text = proc_cgroup if proc_cgroup is not None else _read(Path("/proc/self/cgroup"))
    best: Optional[int] = None
    candidates = [(d, "memory.max") for d in _cgroup_dirs(text, "", root)]
    candidates += [(d, "memory.limit_in_bytes") for d in _cgroup_dirs(text, "memory", root)]
    for folder, name in candidates:
        value = _read(folder / name)
        if not value.isdigit():
            continue
        limit = int(value)
        if 0 < limit < _UNLIMITED_BYTES:
            best = limit if best is None else min(best, limit)
    return best


def _affinity_cpus() -> Optional[int]:
    try:
        return len(os.sched_getaffinity(0)) or None
    except (AttributeError, OSError):
        return None


def _physical_memory() -> Optional[int]:
    for line in _read(Path("/proc/meminfo")).splitlines():
        if line.startswith("MemTotal:"):
            try:
                return int(line.split()[1]) * 1024
            except (IndexError, ValueError):
                return None
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (AttributeError, ValueError, OSError):
        return None


def _env_int(env: Mapping[str, str], name: str) -> Optional[int]:
    raw = str(env.get(name) or "").strip().upper()
    # SLURM_MEM_* are in MB and may carry a unit suffix (e.g. "4G").
    scale = 1
    if raw[-1:] in ("K", "M", "G", "T"):
        scale = {"K": 1 / 1024, "M": 1, "G": 1024, "T": 1024 * 1024}[raw[-1]]
        raw = raw[:-1]
    try:
        value = int(float(raw) * scale)
    except ValueError:
        return None
    return value if value > 0 else None


def detect_budget(
    env: Optional[Mapping[str, str]] = None,
    *,
    cgroup_root: Path = _CGROUP_ROOT,
    proc_cgroup: Optional[str] = None,
) -> ResourceBudget:
    env = os.environ if env is None else env
    cores: dict[str, int] = {"machine": os.cpu_count() or 1}
    affinity = _affinity_cpus()
    if affinity:
        cores["CPU affinity"] = affinity
    quota = cgroup_cpu_limit(cgroup_root, proc_cgroup)
    if quota:
        cores["cgroup CPU quota"] = max(1, math.ceil(quota - 1e-9))
    slurm_cpus = _env_int(env, "SLURM_CPUS_PER_TASK")
    if slurm_cpus:
        cores["SLURM_CPUS_PER_TASK"] = slurm_cpus

    gb = float(1 << 30)
    memory: dict[str, float] = {}
    physical = _physical_memory()
    if physical:
        memory["machine"] = physical / gb
    cg_mem = cgroup_memory_limit(cgroup_root, proc_cgroup)
    if cg_mem:
        memory["cgroup memory limit"] = cg_mem / gb
    per_node = _env_int(env, "SLURM_MEM_PER_NODE")
    if per_node:
        memory["SLURM_MEM_PER_NODE"] = per_node / 1024
    per_cpu = _env_int(env, "SLURM_MEM_PER_CPU")
    if per_cpu:
        n = slurm_cpus or _env_int(env, "SLURM_CPUS_ON_NODE") or 1
        memory["SLURM_MEM_PER_CPU"] = per_cpu * n / 1024

    core_budget = max(1, min(cores.values()))
    mem_budget = min(memory.values()) if memory else None
    return ResourceBudget(core_budget, mem_budget, cores, memory)


__all__ = ["ResourceBudget", "detect_budget", "cgroup_cpu_limit", "cgroup_memory_limit"]
