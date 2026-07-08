"""Container runtime detection, image acquisition, and DB bootstrap for BACTolomew.

Non-GUI logic used by first-time setup (and testable on its own). Docker is the
tested path this iteration; Singularity/Apptainer is wired but untested.

The unified tools image (docker/Dockerfile) carries every bioinformatics tool as a
conda env; here we detect the runtime, make sure the image is present (pull from a
registry when configured, else build locally), and download the Bakta DB by running
`bakta_db download` inside that image.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Callable, Optional

from pipeline.ToolEnv import TOOLS_IMAGE, project_root

LogFn = Callable[[str], None]


def _log(log: Optional[LogFn], msg: str) -> None:
    (log or print)(msg)


# --------------------------------------------------------------------------- #
# Detection
# --------------------------------------------------------------------------- #
def singularity_bin() -> Optional[str]:
    return shutil.which("singularity") or shutil.which("apptainer")


def docker_available() -> bool:
    if not shutil.which("docker"):
        return False
    try:
        return subprocess.run(
            ["docker", "info"], capture_output=True, text=True, timeout=20
        ).returncode == 0
    except Exception:
        return False


def singularity_available() -> bool:
    sbin = singularity_bin()
    if not sbin:
        return False
    try:
        return subprocess.run(
            [sbin, "--version"], capture_output=True, text=True, timeout=20
        ).returncode == 0
    except Exception:
        return False


def detect_runtimes() -> dict:
    """Which container runtimes are usable right now."""
    return {
        "docker": docker_available(),
        "singularity": singularity_available(),
        "singularity_bin": singularity_bin() or "",
    }


def suggest_runtime() -> Optional[str]:
    r = detect_runtimes()
    if r["docker"] and not r["singularity"]:
        return "docker"
    if r["singularity"] and not r["docker"]:
        return "singularity"
    if r["docker"] and r["singularity"]:
        return "docker"  # both available -> caller should ask; default docker
    return None


# --------------------------------------------------------------------------- #
# Image presence / acquisition
# --------------------------------------------------------------------------- #
def docker_image_present(image: str) -> bool:
    try:
        return subprocess.run(
            ["docker", "image", "inspect", image],
            capture_output=True, text=True,
        ).returncode == 0
    except Exception:
        return False


def sif_present(sif_path: str) -> bool:
    return bool(sif_path) and Path(sif_path).is_file()


def _stream(cmd: list[str], log: Optional[LogFn], *, cwd: Optional[Path] = None) -> int:
    """Run a command, streaming combined output line-by-line to `log`."""
    _log(log, "$ " + " ".join(cmd))
    proc = subprocess.Popen(
        cmd, cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        _log(log, line.rstrip("\n"))
    return proc.wait()


def build_image_local(image: str = TOOLS_IMAGE, *, repo_root: Optional[Path] = None,
                      log: Optional[LogFn] = None) -> None:
    repo = Path(repo_root) if repo_root else project_root()
    dockerfile = repo / "docker" / "Dockerfile"
    if not dockerfile.is_file():
        raise FileNotFoundError(f"Dockerfile not found: {dockerfile}")
    _log(log, f"Building image {image} from {dockerfile} (this can take a while)...")
    rc = _stream(["docker", "build", "-f", str(dockerfile), "-t", image, "."], log, cwd=repo)
    if rc != 0:
        raise RuntimeError(f"docker build failed (exit {rc})")


def pull_image(image: str, registry: str, log: Optional[LogFn] = None) -> None:
    ref = f"{registry.rstrip('/')}/{image}"
    _log(log, f"Pulling {ref} ...")
    if _stream(["docker", "pull", ref], log) != 0:
        raise RuntimeError(f"docker pull failed: {ref}")
    # Tag to the bare local name so the rest of the app refers to `image`.
    if ref != image:
        _stream(["docker", "tag", ref, image], log)


def ensure_docker_image(image: str = TOOLS_IMAGE, *, registry: str = "",
                        repo_root: Optional[Path] = None, log: Optional[LogFn] = None) -> None:
    if docker_image_present(image):
        _log(log, f"Image already present: {image}")
        return
    if registry:
        try:
            pull_image(image, registry, log)
            return
        except Exception as exc:
            _log(log, f"Pull failed ({exc}); falling back to local build.")
    build_image_local(image, repo_root=repo_root, log=log)


def ensure_singularity_image(sif_path: str, image: str = TOOLS_IMAGE, *, registry: str = "",
                             repo_root: Optional[Path] = None, log: Optional[LogFn] = None) -> None:
    # UNTESTED this iteration (no Singularity available locally).
    sbin = singularity_bin() or "singularity"
    if sif_present(sif_path):
        _log(log, f"SIF already present: {sif_path}")
        return
    Path(sif_path).parent.mkdir(parents=True, exist_ok=True)
    if registry:
        ref = f"docker://{registry.rstrip('/')}/{image}"
        if _stream([sbin, "pull", sif_path, ref], log) == 0:
            return
        _log(log, "SIF pull failed; falling back to building from a local docker image.")
    ensure_docker_image(image, registry=registry, repo_root=repo_root, log=log)
    if _stream([sbin, "build", sif_path, f"docker-daemon://{image}"], log) != 0:
        raise RuntimeError("singularity build failed")


def ensure_tools_image(runtime: str, *, image: str = TOOLS_IMAGE, registry: str = "",
                       sif_path: str = "", repo_root: Optional[Path] = None,
                       log: Optional[LogFn] = None) -> None:
    if runtime == "singularity":
        ensure_singularity_image(sif_path, image, registry=registry, repo_root=repo_root, log=log)
    else:
        ensure_docker_image(image, registry=registry, repo_root=repo_root, log=log)


# --------------------------------------------------------------------------- #
# Database bootstrap
# --------------------------------------------------------------------------- #
def download_bakta_db(dest: Path, *, db_type: str = "light", runtime: str = "docker",
                      image: str = TOOLS_IMAGE, sif_path: str = "",
                      log: Optional[LogFn] = None) -> Path:
    """Run `bakta_db download` inside the tools image; returns the DB dir."""
    dest = Path(dest).expanduser().resolve()
    dest.mkdir(parents=True, exist_ok=True)
    db_type = db_type if db_type in ("light", "full") else "light"
    inner = ["conda", "run", "--no-capture-output", "-n", "bakta",
             "bakta_db", "download", "--output", str(dest), "--type", db_type]
    if runtime == "singularity":
        sbin = singularity_bin() or "singularity"
        cmd = [sbin, "exec", "--cleanenv", "--bind", str(dest), "--pwd", str(dest), sif_path, *inner]
    else:
        cmd = ["docker", "run", "--rm", "--user", f"{os.getuid()}:{os.getgid()}",
               "-e", f"HOME={dest}", "-v", f"{dest}:{dest}", "-w", str(dest), image, *inner]
    _log(log, f"Downloading Bakta DB ({db_type}) to {dest} ...")
    if _stream(cmd, log) != 0:
        raise RuntimeError("bakta_db download failed")
    # bakta_db writes to <dest>/db or <dest>/db-<type>
    for cand in (dest / "db", dest / f"db-{db_type}"):
        if cand.is_dir():
            return cand
    return dest


def smoke_test_tool(env: str, exe: str, *args: str, image: str = TOOLS_IMAGE) -> tuple[int, str]:
    """Run one tool in the image (docker) and return (returncode, output)."""
    cmd = ["docker", "run", "--rm", image, "conda", "run", "--no-capture-output", "-n", env, exe, *args]
    res = subprocess.run(cmd, capture_output=True, text=True)
    return res.returncode, (res.stdout or "") + (res.stderr or "")


__all__ = [
    "detect_runtimes", "suggest_runtime", "singularity_bin",
    "docker_available", "singularity_available",
    "docker_image_present", "sif_present",
    "ensure_tools_image", "ensure_docker_image", "ensure_singularity_image",
    "build_image_local", "pull_image",
    "download_bakta_db", "smoke_test_tool",
]
