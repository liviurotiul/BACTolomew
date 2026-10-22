"""Container runtime detection, image acquisition, and DB bootstrap for BACTolomew.

Non-GUI logic used by first-time setup (and testable on its own). Docker is the
tested path; Podman (rootful or rootless, also behind the podman-docker shim) and
Singularity/Apptainer are wired too.

The unified tools image (docker/Dockerfile) carries every bioinformatics tool as a
conda env; here we detect the runtime, make sure the image is present (pull it when the
image reference names a registry or container.registry is set; build it locally only when
container.build_on_missing allows). Databases that the tools download themselves are
installed by pipeline.DbDownload.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import re
import shutil
import signal
import subprocess
import sys
import threading
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Optional

from pipeline.ToolEnv import (
    OCI_RUNTIMES,
    TOOLS_IMAGE,
    TOOLS_IMAGE_PLATFORMS,
    ContainerContext,
    image_has_repository_path,
    oci_exe,
    podman_image_ref,
    project_root,
    singularity_bin,
)

LogFn = Callable[[str], None]


def _log(log: Optional[LogFn], msg: str) -> None:
    (log or print)(msg)


# --------------------------------------------------------------------------- #
# Detection
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class RuntimeStatus:
    """Whether a container runtime is usable and, when not, why: `problem` is one of
    not_installed, permission, not_running, not_usable; `message` is one line saying what
    to do; `details` is the probe's output."""
    runtime: str
    ok: bool
    problem: str = ""
    message: str = ""
    details: str = ""


_PROBE_TIMEOUT = 20
_DOCKER_INSTALL = "https://docs.docker.com/engine/install/"
_DOCKER_DESKTOP_INSTALL = "https://docs.docker.com/desktop/"


def _docker_desktop() -> bool:
    """macOS and Windows run Docker as Docker Desktop (no systemctl, no docker group)."""
    return sys.platform == "darwin" or sys.platform.startswith("win")
_PODMAN_INSTALL = "https://podman.io/docs/installation"
_APPTAINER_INSTALL = "https://apptainer.org/docs/admin/main/installation.html"


def _probe(cmd: list[str]) -> tuple[Optional[int], str]:
    """(exit status, stderr then stdout) of a short probe; status None when it timed out
    or could not start."""
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=_PROBE_TIMEOUT)
    except subprocess.TimeoutExpired:
        return None, f"{' '.join(cmd)} did not answer within {_PROBE_TIMEOUT} s"
    except OSError as exc:
        return None, str(exc)
    return res.returncode, "\n".join(t.strip() for t in (res.stderr, res.stdout) if t and t.strip())


def _first_line(text: str) -> str:
    for line in text.splitlines():
        line = line.strip()
        if line:
            return line[:200]
    return "no output"


def _podman_problem(name: str, out: str) -> RuntimeStatus:
    low = out.lower()
    if any(k in low for k in ("newuidmap", "subuid", "subgid", "uid_map", "user namespace")):
        return RuntimeStatus(name, False, "permission",
                             "Podman cannot map this user's IDs. Add subordinate IDs (sudo usermod "
                             "--add-subuids 100000-165535 --add-subgids 100000-165535 $USER), run "
                             "`podman system migrate` and retry.", out)
    if "permission denied" in low:
        return RuntimeStatus(name, False, "permission",
                             "Podman is installed, but this user may not use it: " + _first_line(out), out)
    return RuntimeStatus(name, False, "not_usable", "Podman is installed but not usable: " + _first_line(out), out)


def docker_status(exe: str = "docker") -> RuntimeStatus:
    """Docker's state from `docker info`: not installed, no permission on the socket, the
    daemon not running, or another failure, each with what to do."""
    desktop = _docker_desktop()
    if not shutil.which(exe):
        return RuntimeStatus("docker", False, "not_installed",
                             f"Docker is not installed. Install Docker Desktop ({_DOCKER_DESKTOP_INSTALL}), "
                             "start it and retry." if desktop else
                             f"Docker is not installed. Install Docker Engine ({_DOCKER_INSTALL}) and retry.")
    rc, out = _probe([exe, "info"])
    if rc == 0:
        return RuntimeStatus("docker", True)
    if engine_info("docker").podman:  # the podman-docker shim
        return _podman_problem("docker", out)
    low = out.lower()
    if rc is None and "did not answer" in low:
        return RuntimeStatus("docker", False, "not_running",
                             "Docker does not answer. Restart Docker Desktop and retry." if desktop else
                             "Docker does not answer. Restart it (sudo systemctl restart docker) and retry.", out)
    if "permission denied" in low:
        return RuntimeStatus("docker", False, "permission",
                             "Docker is installed, but this user may not use it. Add the user to the docker "
                             "group (sudo usermod -aG docker $USER), then log out and back in.", out)
    if any(k in low for k in ("cannot connect", "failed to connect", "daemon running", "daemon is running",
                              "connection refused", "no such file or directory")):
        return RuntimeStatus("docker", False, "not_running",
                             "Docker is installed but not running. Start Docker Desktop, wait until it says "
                             "Engine running, and retry." if desktop else
                             "Docker is installed but not running. Start it (sudo systemctl start docker) "
                             "and retry.", out)
    return RuntimeStatus("docker", False, "not_usable", "Docker is not usable: " + _first_line(out), out)


def podman_status() -> RuntimeStatus:
    if not shutil.which("podman"):
        return RuntimeStatus("podman", False, "not_installed",
                             f"Podman is not installed. Install it ({_PODMAN_INSTALL}) and retry.")
    rc, out = _probe(["podman", "info"])
    return RuntimeStatus("podman", True) if rc == 0 else _podman_problem("podman", out)


def singularity_status() -> RuntimeStatus:
    sbin = singularity_bin()
    if not sbin:
        return RuntimeStatus("singularity", False, "not_installed",
                             f"Apptainer (or Singularity) is not installed. Install it ({_APPTAINER_INSTALL}) "
                             "and retry.")
    rc, out = _probe([sbin, "--version"])
    if rc == 0:
        return RuntimeStatus("singularity", True)
    return RuntimeStatus("singularity", False, "not_usable",
                         f"{Path(sbin).name} is installed but not usable: {_first_line(out)}", out)


def runtime_status(runtime: str) -> RuntimeStatus:
    """The status of config `runtime` (native needs no container engine)."""
    if runtime == "podman":
        return podman_status()
    if runtime == "singularity":
        return singularity_status()
    if runtime == "docker":
        return docker_status()
    return RuntimeStatus(runtime, True)


def docker_available() -> bool:
    return docker_status().ok


def podman_available() -> bool:
    return podman_status().ok


def singularity_available() -> bool:
    return singularity_status().ok


def detect_runtimes() -> dict:
    """Which container runtimes are usable right now."""
    return {
        "docker": docker_available(),
        "podman": podman_available(),
        "singularity": singularity_available(),
        "singularity_bin": singularity_bin() or "",
    }


def suggest_runtime() -> Optional[str]:
    """The first usable runtime of docker, podman, singularity (a caller with several
    usable may ask)."""
    r = detect_runtimes()
    return next((name for name in ("docker", "podman", "singularity") if r[name]), None)


@dataclass(frozen=True)
class EngineInfo:
    exe: str                # the executable to call
    podman: bool = False    # Podman, also when it answers as `docker` (podman-docker)
    rootless: bool = False
    selinux: bool = False   # SELinux enforcing on the host


def selinux_enforcing() -> bool:
    try:
        return Path("/sys/fs/selinux/enforce").read_text().strip() == "1"
    except OSError:
        return False


_ENGINE_INFO: dict[str, EngineInfo] = {}


def engine_info(runtime: str) -> EngineInfo:
    """What runs `runtime`'s containers: the executable, whether
    it is Podman (runtime podman, or `docker --version` naming podman: the podman-docker
    shim, called as podman from then on), whether it runs
    rootless (`podman info` Host.Security.Rootless; `docker info` SecurityOptions holding
    name=rootless) and whether SELinux is enforcing. A failed probe reads as rootful and
    is probed again next time (the daemon may start later in this process); a successful
    one is kept for the process."""
    if runtime in _ENGINE_INFO:
        return _ENGINE_INFO[runtime]
    selinux = selinux_enforcing()
    if runtime == "singularity":
        return EngineInfo(singularity_bin() or "singularity", selinux=selinux)
    exe = "podman" if runtime == "podman" else "docker"
    if not shutil.which(exe):
        return EngineInfo(exe, podman=runtime == "podman", selinux=selinux)
    podman = runtime == "podman" or "podman" in _probe([exe, "--version"])[1].lower()
    if podman:
        exe = "podman" if shutil.which("podman") else exe
        rc, out = _probe([exe, "info", "--format", "{{.Host.Security.Rootless}}"])
        rootless = out.strip().lower() == "true" if rc == 0 else os.getuid() != 0
    else:
        rc, out = _probe([exe, "info", "--format", "{{json .SecurityOptions}}"])
        rootless = rc == 0 and "name=rootless" in out
    info = EngineInfo(exe, podman=podman, rootless=rootless, selinux=selinux)
    if rc == 0:
        _ENGINE_INFO[runtime] = info
    return info


def engine_exe(runtime: str) -> str:
    """The executable of a docker or podman runtime (docker for anything else)."""
    return engine_info(runtime).exe if runtime in OCI_RUNTIMES else "docker"


# --------------------------------------------------------------------------- #
# Image presence / acquisition
# --------------------------------------------------------------------------- #
def _engine_ref(image: str, engine: str) -> str:
    return podman_image_ref(image) if Path(engine).name == "podman" else image


# The image's platforms (docker/build.sh), by the names engines and `uname -m` report.
_PLATFORMS = {"amd64": "linux/amd64", "x86_64": "linux/amd64",
              "arm64": "linux/arm64", "aarch64": "linux/arm64"}
# fastGEAR (an x86_64 MATLAB runtime) is in the tools image only on this platform.
FASTGEAR_PLATFORM = "linux/amd64"


def engine_platform(engine: str = "docker") -> str:
    """The platform the engine runs natively, linux/amd64 or linux/arm64 (Docker Desktop
    on Apple Silicon: linux/arm64), so a pull takes that variant of a multi-platform tag and
    a build never runs under emulation. Falls back to this machine's architecture, then to
    linux/amd64."""
    probe = (["info", "--format", "{{.Host.Arch}}"] if Path(engine).name == "podman"
             else ["version", "--format", "{{.Server.Arch}}"])
    rc, out = _probe([engine, *probe])
    arch = out.strip().splitlines()[-1].strip().lower() if rc == 0 and out.strip() else ""
    return _PLATFORMS.get(arch) or _PLATFORMS.get(platform.machine().lower(), "linux/amd64")


def docker_image_present(image: str, engine: str = "docker") -> bool:
    try:
        return subprocess.run(
            [engine, "image", "inspect", _engine_ref(image, engine)],
            capture_output=True, text=True,
        ).returncode == 0
    except Exception:
        return False


# --------------------------------------------------------------------------- #
# Digest pinning: container.image may be `name[:tag]@sha256:<hex>`
# --------------------------------------------------------------------------- #
def split_image_digest(image: str) -> tuple[str, str]:
    """`name[:tag]@sha256:x` -> (`name[:tag]`, `sha256:x`); an unpinned ref -> (ref, "")."""
    name, sep, digest = str(image).strip().partition("@")
    return (name, digest) if sep and digest.startswith("sha256:") else (str(image).strip(), "")


def _image_id(raw_id: Any) -> str:
    """An image ID as Docker prints it: Podman prints the bare hex."""
    text = str(raw_id or "")
    return f"sha256:{text}" if re.fullmatch(r"[0-9a-f]{64}", text) else text


def docker_image_info(image: str, engine: str = "docker") -> Optional[dict]:
    """`docker image inspect` (or podman) of one local image: Id, RepoDigests, RepoTags,
    Labels (None if absent). Podman's bare-hex Id gets the sha256: prefix, and its Labels
    may sit at the top level."""
    try:
        res = subprocess.run(
            [engine, "image", "inspect", "--format", "{{json .}}", _engine_ref(image, engine)],
            capture_output=True, text=True, timeout=30,
        )
    except Exception:
        return None
    if res.returncode != 0:
        return None
    try:
        raw = json.loads(res.stdout.strip().splitlines()[0])
    except Exception:
        return None
    config = raw.get("Config") if isinstance(raw.get("Config"), dict) else {}
    return {
        "id": _image_id(raw.get("Id")),
        "repo_digests": list(raw.get("RepoDigests") or []),
        "repo_tags": list(raw.get("RepoTags") or []),
        "created": raw.get("Created") or "",
        "labels": dict(config.get("Labels") or raw.get("Labels") or {}),
        "platform": _PLATFORMS.get(str(raw.get("Architecture") or "").lower(), ""),
    }


def unpublished_platform(image: str, engine: str = "docker") -> str:
    """The engine's platform when `image` is the default TOOLS_IMAGE and that is not
    published for it (Apple Silicon until an arm64 image exists), else "". There the
    local build stands in for the pinned digest."""
    if image != TOOLS_IMAGE:
        return ""
    target = engine_platform(engine)
    return "" if target in TOOLS_IMAGE_PLATFORMS else target


def foreign_platform(image: str, engine: str = "docker") -> str:
    """The local image's platform when it is not the engine's own (an amd64 image on
    Apple Silicon runs emulated), else "" (also when absent or unknown)."""
    info = docker_image_info(image, engine)
    plat = (info or {}).get("platform") or ""
    return plat if plat and plat != engine_platform(engine) else ""


def _matches_digest(info: Optional[dict], digest: str) -> bool:
    if not info or not digest:
        return False
    return info.get("id") == digest or any(
        str(rd).endswith("@" + digest) for rd in info.get("repo_digests") or []
    )


def pinned_digest_mismatch(image: str, engine: str = "docker") -> str:
    """A warning when `image` pins a digest that the local image does not have, else "".

    The local image is the pinned ref itself or, failing that, the same name[:tag]
    without the digest (a local build or an older pull)."""
    name, digest = split_image_digest(image)
    if not digest:
        return ""
    info = docker_image_info(image, engine) or docker_image_info(name, engine)
    if info is None or _matches_digest(info, digest):
        return ""
    have = ", ".join(info.get("repo_digests") or []) or "none (built locally)"
    return (
        f"Local image {name} does not match the pinned digest {digest}: "
        f"its ID is {info.get('id') or '?'} and its repo digests are {have}. "
        f"Pull {image} to get the pinned image."
    )


def resolve_pinned_image(image: str, log: Optional[LogFn] = None, engine: str = "docker") -> str:
    """The ref `docker run` should use. Unpinned refs are returned unchanged. A pinned
    ref is kept when that exact image is present, unless it is for another platform and a
    native name[:tag] is present (a local build on Apple Silicon); when only name[:tag] is
    present with another digest, warn and use it."""
    name, digest = split_image_digest(image)
    if not digest:
        return image
    if docker_image_present(image, engine):
        if foreign_platform(image, engine) and docker_image_present(name, engine) \
                and not foreign_platform(name, engine):
            _log(log, f"Using the local native image {name}, not the emulated {image}.")
            return name
        return image
    if docker_image_present(name, engine):
        warning = "" if unpublished_platform(image, engine) else pinned_digest_mismatch(image, engine)
        if warning:
            _log(log, "WARNING: " + warning + " Using the local image.")
        return name
    return image


def sif_present(sif_path: str) -> bool:
    return bool(sif_path) and Path(sif_path).is_file()


# The command `_stream` runs on each thread, so another thread can stop it (stop_stream).
_RUNNING: dict[int, subprocess.Popen] = {}
_RUNNING_LOCK = threading.Lock()


def _signal_group(proc: subprocess.Popen, sig: int) -> None:
    """Send `sig` to the command and everything it started (its own process group)."""
    try:
        os.killpg(proc.pid, sig)
    except (AttributeError, ProcessLookupError, PermissionError):
        try:
            proc.send_signal(sig)
        except ProcessLookupError:
            pass


def _group_alive(proc: subprocess.Popen) -> bool:
    proc.poll()  # reap the command itself, or its zombie counts as alive
    try:
        os.killpg(proc.pid, 0)
    except (AttributeError, ProcessLookupError):
        return proc.poll() is None
    except PermissionError:
        return True
    return True


def _end_process(proc: subprocess.Popen, grace: float = 10.0) -> None:
    """SIGTERM to the command's process group (a tool under `conda run`, or `docker run`,
    which passes it on to its container), SIGKILL to what still runs `grace` seconds later.
    Returns once the command and the processes it started have ended."""
    if not _group_alive(proc):
        return
    _signal_group(proc, signal.SIGTERM)
    deadline = time.monotonic() + grace
    while _group_alive(proc):
        if time.monotonic() >= deadline:
            _signal_group(proc, signal.SIGKILL)
            break
        time.sleep(0.05)
    proc.wait()


def _stream(cmd: list[str], log: Optional[LogFn], *, cwd: Optional[Path] = None,
            env: Optional[dict[str, str]] = None) -> int:
    """Run a command, streaming combined output line-by-line to `log`. When `log` raises
    (a cancelled job), the command and what it started are ended before the exception goes
    on. The command runs in its own process group so stop_stream reaches all of it."""
    _log(log, "$ " + " ".join(cmd))
    proc = subprocess.Popen(
        cmd, cwd=str(cwd) if cwd else None, env=env,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1,
        start_new_session=True,
    )
    assert proc.stdout is not None
    key = threading.get_ident()
    with _RUNNING_LOCK:
        _RUNNING[key] = proc
    try:
        for line in proc.stdout:
            _log(log, line.rstrip("\n"))
        return proc.wait()
    except BaseException:
        _end_process(proc)
        raise
    finally:
        proc.stdout.close()
        with _RUNNING_LOCK:
            _RUNNING.pop(key, None)


def stop_stream(thread_id: int, grace: float = 10.0) -> bool:
    """End the command `_stream` runs on thread `thread_id` and what it started: SIGTERM
    to its process group, SIGKILL when any of it still runs `grace` seconds later. The
    streaming thread then sees the output end."""
    with _RUNNING_LOCK:
        proc = _RUNNING.get(thread_id)
    if proc is None or not _group_alive(proc):
        return False
    threading.Thread(target=_end_process, args=(proc, grace), daemon=True).start()
    return True


def _verify_archive(archive: Path, log: Optional[LogFn]) -> None:
    """Check `archive` against `<archive>.sha256` (sha256sum format) when that file exists."""
    sums = archive.with_name(archive.name + ".sha256")
    if not sums.is_file():
        return
    want = (sums.read_text(encoding="utf-8").split() or [""])[0].lower()
    have = _sha256_file(archive)
    if want != have:
        raise RuntimeError(f"{archive} does not match {sums.name} (sha256 {have}, expected {want}).")
    _log(log, f"Checked {archive.name} against {sums.name}.")


# Image IDs of digest-pinned images, for a local copy without its registry digest (a
# `docker load` of a saved archive): registry digest -> image ID (config digest).
KNOWN_IMAGE_IDS = {
    # pgcbioinfo/fastgear:latest, pushed 2022-12-01 (DEFAULT_FASTGEAR_IMAGE)
    "sha256:afb299e7fb327129ebf2acc1bcae00e0a54cd6fbe71d114d0a9f1642a2216b79":
        "sha256:ac1746a455a64b9977edfe679b28afbe8ff87d8a6a8c681a5a7844dab48b03fb",
}


def _local_copy_ok(ref: str, image: str, engine: str, log: Optional[LogFn]) -> bool:
    """Whether the local `ref` may stand in for `image`: present and, when `image` pins a
    digest, that image (by repo digest, or by the image ID KNOWN_IMAGE_IDS records)."""
    if not docker_image_present(ref, engine):
        return False
    digest = split_image_digest(image)[1]
    if not digest:
        return True
    info = docker_image_info(ref, engine)
    if _matches_digest(info, digest) or (info or {}).get("id") == KNOWN_IMAGE_IDS.get(digest, "-"):
        return True
    _log(log, f"Local image {ref} (ID {(info or {}).get('id') or '?'}) is not the pinned {image}: not used.")
    return False


def ensure_fastgear_image(image: str, archive: str = "", *, engine: str = "docker",
                          log: Optional[LogFn] = None) -> str:
    """Make the fastGEAR base image present for a local build of the tools image and
    return the local reference to build from: `image` when present (or its name[:tag],
    say after a `docker load`, which drops the registry digest, when it is the pinned
    image: see _local_copy_ok), else loaded from `archive` (config
    container.fastgear_archive, a `docker save` file, checked against `<archive>.sha256`
    when present), else pulled. Raises when none works."""
    name = split_image_digest(image)[0]
    for ref in (image, name):
        if _local_copy_ok(ref, image, engine, log):
            return ref
    if archive:
        path = Path(archive).expanduser()
        if not path.is_file():
            raise RuntimeError(f"container.fastgear_archive {path} does not exist.")
        _verify_archive(path, log)
        _log(log, f"Loading the fastGEAR image from {path} ...")
        lines: list[str] = []
        if _stream([engine, "load", "-i", str(path)], lambda line: (lines.append(line), _log(log, line))) != 0:
            raise RuntimeError(f"{engine} load -i {path} failed")
        loaded = [ln.split(":", 1)[1].strip() for ln in lines if ln.startswith("Loaded image:")]
        for ref in (image, name, *loaded):
            if _local_copy_ok(ref, image, engine, log):
                return ref
        raise RuntimeError(f"{path} loaded, but it does not hold {image}.")
    pull_ref_ = _engine_ref(image, engine)
    _log(log, f"Pulling {pull_ref_} ...")
    if _stream([engine, "pull", "--platform", FASTGEAR_PLATFORM, pull_ref_], log) != 0:
        raise RuntimeError(
            f"The fastGEAR image {image} is not available locally and could not be pulled. Set "
            "container.fastgear_archive to a saved copy (docker save) to build without Docker Hub."
        )
    return image


def build_image_local(image: str = TOOLS_IMAGE, *, repo_root: Optional[Path] = None,
                      log: Optional[LogFn] = None, engine: str = "docker",
                      fastgear_image: str = "", fastgear_archive: str = "") -> None:
    """Build the tools image from docker/Dockerfile for the engine's own platform
    (engine_platform). On linux/amd64 with `fastgear_image` (and `fastgear_archive`) the
    fastGEAR stage is built from that image, made present first; on linux/arm64 the image
    is built without fastGEAR (FastBAPS is its population method), as docker/build.sh does."""
    repo = Path(repo_root) if repo_root else project_root()
    dockerfile = repo / "docker" / "Dockerfile"
    if not dockerfile.is_file():
        raise FileNotFoundError(f"Dockerfile not found: {dockerfile}")
    image = split_image_digest(image)[0]  # a build has no registry digest to tag with
    target = engine_platform(engine)
    build_args: list[str] = []
    if target != FASTGEAR_PLATFORM:
        _log(log, f"Building for {target}: fastGEAR (x86_64 only) is left out.")
        build_args = ["--build-arg", "FASTGEAR=none", "--build-arg", "FASTGEAR_PROVENANCE="]
    elif fastgear_image:
        base = ensure_fastgear_image(fastgear_image, fastgear_archive, engine=engine, log=log)
        # The label records the pinned ref (base was checked against it), or for an
        # unpinned one the image ID it resolved to.
        provenance = fastgear_image
        if not split_image_digest(fastgear_image)[1]:
            provenance += f" (image ID {(docker_image_info(base, engine) or {}).get('id') or '?'})"
        build_args = ["--build-arg", f"FASTGEAR_IMAGE={_engine_ref(base, engine)}",
                      "--build-arg", f"FASTGEAR_PROVENANCE={provenance}"]
    _log(log, f"Building image {image} from {dockerfile} (this can take a while)...")
    rc = _stream([engine, "build", "--platform", target, "-f", str(dockerfile), *build_args,
                  "-t", image, "."], log, cwd=repo)
    if rc != 0:
        raise RuntimeError(f"{engine} build failed (exit {rc})")


def pull_ref(image: str, registry: str = "") -> str:
    """The reference to pull `image` from: the image itself when it names a repository
    path (`ghcr.io/owner/name:tag@sha256:...`), `registry/image` for a bare name when
    `registry` is set, else "" (a bare local name with no registry cannot be pulled)."""
    if image_has_repository_path(image):
        return image
    return f"{registry.rstrip('/')}/{image}" if registry else ""


def pull_image(image: str, registry: str = "", log: Optional[LogFn] = None, engine: str = "docker") -> None:
    ref = _engine_ref(pull_ref(image, registry) or image, engine)
    _log(log, f"Pulling {ref} ...")
    if _stream([engine, "pull", "--platform", engine_platform(engine), ref], log) != 0:
        raise RuntimeError(f"{engine} pull failed: {ref}")
    # Tag to the bare local name so the rest of the app refers to `image` (a digest
    # cannot be a tag: a pinned image is tagged with its name[:tag] part).
    target = split_image_digest(image)[0]
    if ref != target:
        _stream([engine, "tag", ref, target], log)


def _no_build_error(image: str, pull_error: str) -> RuntimeError:
    why = pull_error or f"{image} names no registry and container.registry is empty"
    return RuntimeError(
        f"The tools image {image} is not available: {why}. "
        "container.build_on_missing is false, so it is not built locally. Check the image "
        "reference and the network (a private package needs `docker login`), or set "
        "container.build_on_missing to true to build it from docker/Dockerfile (about an "
        "hour; a local build does not match a pinned digest)."
    )


def ensure_docker_image(image: str = TOOLS_IMAGE, *, registry: str = "",
                        repo_root: Optional[Path] = None, log: Optional[LogFn] = None,
                        build_on_missing: bool = True, engine: str = "docker",
                        fastgear_image: str = "", fastgear_archive: str = "") -> None:
    """Make `image` present (with docker or podman, `engine`): pull it when it names a
    registry (or `registry` is set), else, or when the pull fails, build it from
    docker/Dockerfile, but only when `build_on_missing` (config container.build_on_missing)
    is true; otherwise raise. A build takes its fastGEAR stage from `fastgear_image`,
    loaded from `fastgear_archive` when missing."""
    name, digest = split_image_digest(image)
    if docker_image_present(image, engine):
        _log(log, f"Image already present: {image}")
        plat = foreign_platform(image, engine)
        if plat:
            _log(log, f"Note: it is a {plat} image, so its tools run emulated. Remove it "
                      f"(`{engine} image rm {image}`) and rerun setup to build a native one.")
        return
    ref = pull_ref(image, registry)
    if ref and ref != image and docker_image_present(ref, engine):
        _log(log, f"Image already present: {ref}")
        return
    pull_error = ""
    unpublished = unpublished_platform(image, engine)
    if ref and unpublished:
        # A single-platform digest: the pull would download it all and then fail.
        _log(log, f"The published tools image is for {', '.join(TOOLS_IMAGE_PLATFORMS)} only "
                  f"and this engine runs {unpublished}: not pulling it.")
        pull_error = f"no {unpublished} variant is published"
        ref = ""
    if ref:
        try:
            pull_image(image, registry, log, engine)
            return
        except Exception as exc:
            pull_error = str(exc)
            _log(log, f"Pull failed ({exc}).")
    if not build_on_missing:
        raise _no_build_error(image, pull_error)
    if digest and docker_image_present(name, engine):
        _log(log, f"Using the local image {name}.")
    else:
        _log(log, "container.build_on_missing is true: building the image locally.")
        build_image_local(name, repo_root=repo_root, log=log, engine=engine,
                          fastgear_image=fastgear_image, fastgear_archive=fastgear_archive)
    warning = "" if unpublished else pinned_digest_mismatch(image, engine)
    if warning:
        _log(log, "WARNING: " + warning)


def ensure_singularity_image(sif_path: str, image: str = TOOLS_IMAGE, *, registry: str = "",
                             repo_root: Optional[Path] = None, log: Optional[LogFn] = None,
                             build_on_missing: bool = True, fastgear_image: str = "",
                             fastgear_archive: str = "") -> None:
    """Make the SIF at `sif_path` present: `singularity pull docker://<ref>` when the image
    names a registry (no Docker needed, the HPC case), else convert a local Docker image,
    which needs a Docker daemon and is built only when `build_on_missing` is true."""
    # UNTESTED this iteration (no Singularity available locally).
    sbin = singularity_bin() or "singularity"
    if sif_present(sif_path):
        _log(log, f"SIF already present: {sif_path}")
        return
    if not sif_path:
        raise RuntimeError("container.sif_path is empty: set where the .sif file should be written.")
    Path(sif_path).parent.mkdir(parents=True, exist_ok=True)
    ref = pull_ref(image, registry)
    pull_error = ""
    if ref:
        if _stream([sbin, "pull", sif_path, f"docker://{ref}"], log) == 0:
            return
        pull_error = f"singularity pull docker://{ref} failed"
        _log(log, pull_error + ".")
    if not build_on_missing:
        raise _no_build_error(image, pull_error)
    docker = docker_status()
    if not docker.ok:
        raise RuntimeError(
            (pull_error + "; " if pull_error else f"{image} names no registry; ")
            + "building the SIF from a local Docker image needs a Docker daemon, which is not "
            f"usable ({docker.message}). Pull the image on a machine with network access "
            "(`singularity pull <file>.sif docker://<image>`) and set container.sif_path."
        )
    _log(log, "Converting a local Docker image to a SIF.")
    ensure_docker_image(image, registry=registry, repo_root=repo_root, log=log,
                        build_on_missing=build_on_missing, fastgear_image=fastgear_image,
                        fastgear_archive=fastgear_archive)
    if _stream([sbin, "build", sif_path, f"docker-daemon://{resolve_pinned_image(image, log)}"], log) != 0:
        raise RuntimeError("singularity build failed")


def ensure_tools_image(runtime: str, *, image: str = TOOLS_IMAGE, registry: str = "",
                       sif_path: str = "", repo_root: Optional[Path] = None,
                       log: Optional[LogFn] = None, build_on_missing: bool = True,
                       fastgear_image: str = "", fastgear_archive: str = "") -> None:
    if runtime == "singularity":
        ensure_singularity_image(sif_path, image, registry=registry, repo_root=repo_root, log=log,
                                 build_on_missing=build_on_missing, fastgear_image=fastgear_image,
                                 fastgear_archive=fastgear_archive)
    else:
        ensure_docker_image(image, registry=registry, repo_root=repo_root, log=log,
                            build_on_missing=build_on_missing, engine=engine_exe(runtime),
                            fastgear_image=fastgear_image, fastgear_archive=fastgear_archive)


def smoke_test_tool(env: str, exe: str, *args: str, image: str = TOOLS_IMAGE,
                    engine: str = "docker") -> tuple[int, str]:
    """Run one tool in the image (docker or podman) and return (returncode, output)."""
    cmd = [engine, "run", "--rm", _engine_ref(image, engine), "conda", "run", "--no-capture-output", "-n", env, exe, *args]
    res = subprocess.run(cmd, capture_output=True, text=True)
    return res.returncode, (res.stdout or "") + (res.stderr or "")


# --------------------------------------------------------------------------- #
# Run provenance: which image and tool versions a run used
# --------------------------------------------------------------------------- #
IMAGE_TOOL_VERSIONS_PATH = "/opt/bactolomew/tool_versions.json"
# Written at build time by docker/verify_image.sh (images after 0.2.0): the in-image
# mob_suite database's date and checksums.
IMAGE_CONTENTS_PATH = "/opt/bactolomew/image_contents.json"
# Older images have only mob_init's status line ("Download date: YYYY-MM-DD. ...").
IMAGE_MOBSUITE_STATUS_GLOB = "/opt/conda/envs/mobsuite/lib/python3.*/site-packages/mob_suite/databases/status.txt"
_SECTION = "@@bactolomew-section@@"
_IMAGE_READ_SCRIPT = (
    f"cat {IMAGE_TOOL_VERSIONS_PATH} 2>/dev/null; echo; echo {_SECTION}; "
    f"cat {IMAGE_CONTENTS_PATH} 2>/dev/null; echo; echo {_SECTION}; "
    f"cat {IMAGE_MOBSUITE_STATUS_GLOB} 2>/dev/null | head -n 1"
)
_OCI_LABELS = ("org.opencontainers.image.version", "org.opencontainers.image.revision",
               "org.opencontainers.image.created", "org.opencontainers.image.source")
# The fastGEAR image the tools image's fastGEAR stage was built from (docker/Dockerfile).
FASTGEAR_BASE_LABEL = "org.bactolomew.fastgear.image"
_DATE_RE = re.compile(r"\d{4}-\d{2}-\d{2}")


def _compact_tool_versions(raw: Any) -> Optional[dict]:
    """tool_versions.json -> {env: {tool: "version build"}} plus the fastGEAR identity."""
    if not isinstance(raw, dict) or not isinstance(raw.get("envs"), dict):
        return None
    envs: dict[str, dict[str, str]] = {}
    for env, entry in raw["envs"].items():
        tools = entry.get("tools") if isinstance(entry, dict) else None
        if isinstance(tools, dict):
            envs[env] = {
                name: " ".join(str(t.get(k)) for k in ("version", "build") if t.get(k))
                for name, t in tools.items() if isinstance(t, dict)
            }
    out: dict[str, Any] = {"source": IMAGE_TOOL_VERSIONS_PATH, "envs": envs}
    fastgear = raw.get("fastgear")
    if isinstance(fastgear, dict) and fastgear.get("fastGEAR_sha256"):
        out["fastgear_sha256"] = fastgear["fastGEAR_sha256"]
    return out


def _json_or_none(text: str) -> Any:
    try:
        return json.loads(text) if text.strip() else None
    except Exception:
        return None


def _image_mobsuite_db(contents: Any, status_line: str) -> Optional[dict]:
    """The in-image mob_suite database: image_contents.json when the image has it, else
    the date in mob_init's status.txt."""
    if isinstance(contents, dict) and isinstance(contents.get("mobsuite_db"), dict):
        return {"location": "image", "source": IMAGE_CONTENTS_PATH, **contents["mobsuite_db"]}
    m = _DATE_RE.search(status_line or "")
    if m:
        return {"location": "image", "download_date": m.group(0),
                "source": "mob_suite databases/status.txt (no checksums recorded in this image)"}
    return None


def read_image_metadata(ctx: ContainerContext) -> dict:
    """The image's tool_versions.json (compacted) and in-image database identities, read
    in one short container start. {"tool_versions": dict | None, "databases": dict}."""
    out: dict[str, Any] = {"tool_versions": None, "databases": {}}
    if ctx.runtime == "singularity":
        if not sif_present(ctx.sif_path):
            return out
        cmd = [singularity_bin() or "singularity", "exec", "--cleanenv", ctx.sif_path,
               "sh", "-c", _IMAGE_READ_SCRIPT]
    elif ctx.runtime in OCI_RUNTIMES and ctx.image:
        cmd = [oci_exe(ctx), "run", "--rm", "--network", "none",
               podman_image_ref(ctx.image) if ctx.podman else ctx.image, "sh", "-c", _IMAGE_READ_SCRIPT]
    else:
        return out
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except Exception:
        return out
    if res.returncode != 0:
        return out
    parts = (res.stdout.split(_SECTION) + ["", "", ""])[:3]
    out["tool_versions"] = _compact_tool_versions(_json_or_none(parts[0]))
    mobsuite = _image_mobsuite_db(_json_or_none(parts[1]), parts[2].strip())
    if mobsuite:
        out["databases"]["mobsuite"] = mobsuite
    return out


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def _sif_hash_cache_path() -> Path:
    base = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
    return Path(base) / "bactolomew" / "sif_sha256.json"


def sif_sha256(sif_path: str) -> str:
    """sha256 of a SIF file, cached per path, size and mtime (hashing a 10 GB SIF takes
    about a minute, so it is done once per file version)."""
    path = Path(sif_path).resolve()
    st = path.stat()
    key = str(path)
    stamp = {"size": st.st_size, "mtime_ns": st.st_mtime_ns}
    cache_file = _sif_hash_cache_path()
    try:
        cache = json.loads(cache_file.read_text(encoding="utf-8"))
        cache = cache if isinstance(cache, dict) else {}
    except Exception:
        cache = {}
    hit = cache.get(key)
    if isinstance(hit, dict) and all(hit.get(k) == v for k, v in stamp.items()) and hit.get("sha256"):
        return str(hit["sha256"])
    digest = _sha256_file(path)
    cache[key] = {**stamp, "sha256": digest}
    try:
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        tmp = cache_file.with_suffix(".tmp")
        tmp.write_text(json.dumps(cache, indent=1), encoding="utf-8")
        tmp.replace(cache_file)
    except Exception:
        pass
    return digest


def _find_labels(node: Any) -> dict:
    if isinstance(node, dict):
        if isinstance(node.get("labels"), dict):
            return node["labels"]
        for value in node.values():
            found = _find_labels(value)
            if found:
                return found
    return {}


def sif_labels(sif_path: str) -> dict:
    """The SIF's labels (`singularity inspect --labels --json`); a SIF pulled from
    docker:// keeps the image's OCI labels. {} when unavailable."""
    try:
        res = subprocess.run(
            [singularity_bin() or "singularity", "inspect", "--labels", "--json", sif_path],
            capture_output=True, text=True, timeout=60,
        )
        return _find_labels(json.loads(res.stdout)) if res.returncode == 0 else {}
    except Exception:
        return {}


def _short_labels(labels: dict) -> dict:
    return {k.rsplit(".", 1)[-1]: v for k, v in labels.items() if k in _OCI_LABELS}


def database_provenance(*, bakta_db: Any = "", checkm2_db: Any = "", mobsuite_db_dir: Any = "") -> dict:
    """Identities of the databases outside the image: the Bakta DB's version.json (and the
    AMRFinderPlus DB release it bundles), the CheckM2 DIAMOND file (path, size, mtime; too
    large to hash per run) and an external mob_suite database's status.txt. Never raises."""
    out: dict[str, Any] = {}
    if bakta_db:
        base = Path(str(bakta_db)).expanduser()
        entry: dict[str, Any] = {"path": str(base)}
        for cand in (base, base / "db"):
            version_file = cand / "version.json"
            if version_file.is_file():
                raw = _json_or_none(version_file.read_text(encoding="utf-8", errors="replace"))
                if isinstance(raw, dict):
                    entry.update({k: raw[k] for k in ("date", "major", "minor", "type", "doi") if k in raw})
                amr = cand / "amrfinderplus-db" / "latest" / "version.txt"
                if amr.is_file():
                    entry["amrfinderplus_db"] = amr.read_text(encoding="utf-8", errors="replace").strip()
                break
        else:
            entry["version"] = None
        out["bakta"] = entry
    if checkm2_db:
        path = Path(str(checkm2_db)).expanduser()
        entry = {"path": str(path)}
        try:
            st = path.stat()
            entry["size"] = st.st_size
            entry["mtime"] = datetime.fromtimestamp(st.st_mtime, timezone.utc).isoformat(timespec="seconds")
        except OSError:
            entry["missing"] = True
        out["checkm2"] = entry
    if mobsuite_db_dir:
        path = Path(str(mobsuite_db_dir)).expanduser()
        entry = {"location": "external", "path": str(path)}
        try:
            status = (path / "status.txt").read_text(encoding="utf-8", errors="replace").strip()
            m = _DATE_RE.search(status)
            entry["download_date"] = m.group(0) if m else status
        except OSError:
            entry["download_date"] = None
        out["mobsuite"] = entry
    return out


def environment_provenance(ctx: ContainerContext, *, configured_image: str = "") -> dict:
    """The run's tool environment for run_config.json / run_report.json: runtime, image
    reference, local image ID, repo digests, OCI labels, pinned-digest check (for
    Singularity: the SIF's sha256 and labels), the image's tool versions and its in-image
    databases. Never raises; a missing piece is recorded as missing."""
    env: dict[str, Any] = {"runtime": ctx.runtime}
    if ctx.runtime in OCI_RUNTIMES and (ctx.podman or ctx.rootless):
        env["engine"] = {"exe": oci_exe(ctx), "podman": ctx.podman, "rootless": ctx.rootless}
    native_flag = os.environ.get("BACTOLOMEW_NATIVE")
    if native_flag:
        env["BACTOLOMEW_NATIVE"] = native_flag
    if ctx.runtime == "native":
        env["tool_versions"] = None
        env["tool_versions_source"] = "native conda envs: see run_report.json tools.versions"
        return env
    all_labels: dict = {}
    try:
        if ctx.runtime == "singularity":
            env["sif_path"] = ctx.sif_path
            env["image"] = ctx.image
            if sif_present(ctx.sif_path):
                env["sif_sha256"] = sif_sha256(ctx.sif_path)
                all_labels = sif_labels(ctx.sif_path)
                labels = _short_labels(all_labels)
                if labels:
                    env["image_labels"] = labels
            else:
                env["sif_sha256"] = None
        else:
            env["image"] = ctx.image
            if configured_image and configured_image != ctx.image:
                env["configured_image"] = configured_image
            pinned = split_image_digest(configured_image or ctx.image)[1]
            info = docker_image_info(ctx.image, oci_exe(ctx))
            if info is None:
                env["image_id"] = None
            else:
                env["image_id"] = info["id"]
                env["repo_digests"] = info["repo_digests"]
                env["image_created"] = info["created"]
                all_labels = info["labels"]
                labels = _short_labels(all_labels)
                if labels:
                    env["image_labels"] = labels
            if pinned:
                env["pinned_digest"] = pinned
                env["pinned_digest_matches"] = _matches_digest(info, pinned)
    except Exception as exc:
        env["error"] = f"image inspect failed: {exc}"
    # fastGEAR runs from the tools image; its stage came from this fastGEAR image (images
    # built before the label: the digest docker/Dockerfile pinned, DEFAULT_FASTGEAR_IMAGE).
    env["fastgear_image"] = all_labels.get(FASTGEAR_BASE_LABEL) or None
    meta = read_image_metadata(ctx)
    env["tool_versions"] = meta["tool_versions"]
    if meta["tool_versions"] is None:
        env["tool_versions_source"] = (
            f"{IMAGE_TOOL_VERSIONS_PATH} not in the image: see run_report.json tools.versions"
        )
    if meta["databases"]:
        env["databases"] = meta["databases"]
    return env


__all__ = [
    "detect_runtimes", "suggest_runtime", "singularity_bin",
    "RuntimeStatus", "runtime_status", "docker_status", "podman_status", "singularity_status",
    "EngineInfo", "engine_info", "engine_exe", "selinux_enforcing",
    "docker_available", "podman_available", "singularity_available", "ensure_fastgear_image",
    "docker_image_present", "sif_present",
    "split_image_digest", "docker_image_info", "pinned_digest_mismatch", "resolve_pinned_image",
    "read_image_metadata", "environment_provenance", "database_provenance",
    "sif_sha256", "sif_labels", "pull_ref",
    "ensure_tools_image", "ensure_docker_image", "ensure_singularity_image",
    "build_image_local", "pull_image", "engine_platform", "foreign_platform",
    "unpublished_platform",
    "smoke_test_tool",
]
