"""BACTolomew's single entry point.

    bactolomew [gui] [--dev]          the desktop app (the default)
    bactolomew setup --yes [...]      first-time setup without a window
    bactolomew setup --viewer         a viewer-only install: open shared runs, no downloads
    bactolomew run --genomes DIR --out DIR [--threads N] [--resume]
    bactolomew check                  what this machine has and lacks
    bactolomew db --list | --add FASTA | --remove NAME   marker and compound databases
    bactolomew --version

Only `gui` imports Qt; setup, run, check and db work on a headless node.
"""
from __future__ import annotations

import argparse
import os
import platform
import subprocess
import sys
import threading
from pathlib import Path
from typing import Callable, Optional, Sequence

_SRC = Path(__file__).resolve().parent
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

RUNTIME_CHOICES = ("docker", "podman", "singularity", "apptainer")


def _version_text() -> str:
    from pipeline import __version__

    text = f"BACTolomew {__version__}"
    try:
        rev = subprocess.run(["git", "-C", str(_SRC.parent), "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, timeout=5)
        if rev.returncode == 0 and rev.stdout.strip():
            text += f" (git {rev.stdout.strip()})"
    except (OSError, subprocess.SubprocessError):
        pass
    return text


# --------------------------------------------------------------------------- #
# gui
# --------------------------------------------------------------------------- #
def display_problem(env: Optional[dict] = None) -> str:
    """Why the GUI cannot open a window here, or "" when it can."""
    env = os.environ if env is None else env
    if not sys.platform.startswith("linux"):
        return ""
    if env.get("DISPLAY") or env.get("WAYLAND_DISPLAY") or env.get("QT_QPA_PLATFORM") == "offscreen":
        return ""
    return ("No display: DISPLAY and WAYLAND_DISPLAY are not set, so the BACTolomew window cannot open.\n"
            "On a server or cluster node use `bactolomew run` (pipeline) or `bactolomew check`; "
            "for the window, connect with a display (ssh -X, a remote desktop) or run on your workstation.")


def cmd_gui(args: argparse.Namespace) -> int:
    problem = display_problem()
    if problem:
        print(problem, file=sys.stderr)
        return 1
    try:
        from annotation_explorer_polars.app import main as app_main
    except ImportError as exc:
        print(f"The BACTolomew window cannot start: {exc}\n"
              "Run `bactolomew check` to see what is missing.", file=sys.stderr)
        return 1
    sys.argv = [sys.argv[0]]  # QApplication reads sys.argv
    return int(app_main(dev_mode=bool(args.dev)) or 0)


# --------------------------------------------------------------------------- #
# setup
# --------------------------------------------------------------------------- #
NO_RUNTIME = ("No usable container runtime: install Docker, Podman or Apptainer (Singularity), "
              "or pass --runtime. `bactolomew check` says what was found.")


def default_sif_path(root: Path, config: dict) -> Path:
    """Where setup writes the tools SIF when the config names none: the install folder,
    named after the image (bactolomew-0.3.0.sif)."""
    from pipeline.ToolEnv import TOOLS_IMAGE

    container = config.get("container") if isinstance(config.get("container"), dict) else {}
    image = str(container.get("image") or "").strip() or TOOLS_IMAGE
    return root / (image.split("@")[0].rsplit("/", 1)[-1].replace(":", "-") + ".sif")


def _setup_config(args: argparse.Namespace, current: dict, root: Path) -> tuple[dict, dict, str]:
    """(config the setup runs with, keys to write, problem). --runtime / --sif win; a config
    without a runtime gets the first usable one found here (docker, for the plan of a dry
    run, when none is). The problem says why the image step cannot work ("" when it can)."""
    from pipeline import ContainerRuntime as CR

    updates: dict = {}
    problem = ""
    runtime = args.runtime
    if runtime == "apptainer":
        runtime = "singularity"
    if not runtime and args.sif:
        runtime = "singularity"
    if not runtime and not current.get("runtime"):
        runtime = CR.suggest_runtime()
        if runtime:
            print(f"Container runtime: {runtime} (detected)")
        else:
            runtime, problem = "docker", NO_RUNTIME
    if runtime and runtime != current.get("runtime") and problem != NO_RUNTIME:
        updates["runtime"] = runtime
    effective = runtime or str(current.get("runtime") or "").strip().lower()
    effective = "singularity" if effective == "apptainer" else effective
    container = current.get("container") if isinstance(current.get("container"), dict) else {}
    if args.sif:
        updates["container"] = {"sif_path": str(Path(args.sif).expanduser().resolve())}
    elif effective == "singularity" and not str(container.get("sif_path") or "").strip():
        updates["container"] = {"sif_path": str(default_sif_path(root, current))}
    if not problem and effective in ("docker", "podman", "singularity"):
        status = CR.runtime_status(effective)
        if not status.ok:
            problem = f"{effective}: {status.message}"
    merged = {**current, "runtime": effective} if effective else dict(current)
    for key, value in updates.items():
        merged[key] = {**(merged.get(key) or {}), **value} if isinstance(value, dict) else value
    return merged, updates, problem


class _TextProgress:
    """Setup steps as text: one line per state change; on a terminal the running step's
    progress is rewritten in place, elsewhere every tenth percent."""

    def __init__(self, total: int, verbose: bool) -> None:
        self.total = total
        self.verbose = verbose
        self.tty = sys.stdout.isatty()
        self.index: dict[str, int] = {}
        self.last: dict[str, tuple[str, str]] = {}
        self.tail: list[str] = []
        self._inline = False
        self._lock = threading.Lock()

    def _emit(self, text: str, inline: bool = False) -> None:
        if self._inline and not inline:
            sys.stdout.write("\n")
        sys.stdout.write(("\r\033[K" if inline else "") + text + ("" if inline else "\n"))
        sys.stdout.flush()
        self._inline = inline

    def step(self, step) -> None:
        with self._lock:
            n = self.index.setdefault(step.key, len(self.index) + 1)
            head = f"[{n}/{self.total}] {step.label}"
            before = self.last.get(step.key, ("", ""))
            self.last[step.key] = (step.state, step.text)
            if (step.state, step.text) == before:
                return
            if step.state == "running" and before[0] == "running":
                pct = step.text.rstrip("%")
                if self.tty:
                    self._emit(f"{head}: {step.text}", inline=True)
                elif not pct.isdigit() or int(pct) % 10 == 0:
                    self._emit(f"{head}: {step.text}")
                return
            line = f"{head}: {step.text}"
            if step.state == "failed" and step.note:
                line += f" ({step.note})"
            elif step.state == "done" and step.detail:
                line += f" ({step.detail})"
            self._emit(line)

    def log(self, line: str) -> None:
        with self._lock:
            self.tail = (self.tail + [line])[-40:]
            if self.verbose:
                self._emit("    " + line)


def cmd_setup(args: argparse.Namespace) -> int:
    from pipeline.Setup import SetupJob, _read_config, write_viewer_config
    from pipeline.ToolEnv import project_root, save_app_config

    root = project_root()
    if getattr(args, "viewer", False):
        if args.dry_run:
            print(f"[dry-run] would set mode \"viewer\" in {root / 'bactolomew_config.json'} (no downloads)")
            return 0
        print(f"Viewer install: {write_viewer_config(root)} says mode \"viewer\".")
        print("Open shared runs with Import Run... in the run list. To run predictions later: "
              "`bactolomew setup`.")
        return 0
    config, updates, problem = _setup_config(args, _read_config(root), root)
    if args.dry_run:
        print(f"[dry-run] project folder: {root}")
        if updates:
            print(f"[dry-run] would set in bactolomew_config.json: {updates}")
        job = SetupJob(root, config=config, log=lambda _line: None)
        for line in job.dry_run():
            print(f"[dry-run] {line}")
        if problem:
            print(f"[dry-run] setup cannot install the tool image yet: {problem}", file=sys.stderr)
            return 1
        return 0
    if problem:
        print(problem, file=sys.stderr)
        return 1
    if not args.yes:
        if not sys.stdin.isatty():
            print("Setup downloads the tool image and databases (over 40 GB; about 125 GB free needed). "
                  "Pass --yes to go ahead, or --dry-run to see the plan.", file=sys.stderr)
            return 2
        if input("Download the tool image and databases (over 40 GB) now? [y/N] ").strip().lower() not in ("y", "yes"):
            return 2
    if updates:
        save_app_config(updates)
        print(f"Set in bactolomew_config.json: {updates}")

    progress: _TextProgress
    job = SetupJob(root, on_step=lambda s: progress.step(s), log=lambda l: progress.log(l))
    progress = _TextProgress(len(job.steps), args.verbose)
    outcome: list[str] = []

    def work() -> None:
        job.thread_id = threading.get_ident()
        try:
            outcome.append(job.run())
        except Exception as exc:  # noqa: BLE001 - a bug in a step: shown as a failure
            progress.log(f"ERROR: {exc}")
            outcome.append("failed")
        finally:
            job.thread_id = None

    worker = threading.Thread(target=work, name="bactolomew-setup", daemon=True)
    worker.start()
    try:
        while worker.is_alive():
            worker.join(0.5)
    except KeyboardInterrupt:
        print("\nStopping the setup (the download in progress is discarded)...", file=sys.stderr)
        job.cancel()
        worker.join()
    result = outcome[0] if outcome else "cancelled"
    if result == "done":
        print("BACTolomew is ready.")
        return 0
    if result == "cancelled":
        print("Setup stopped. Run it again to continue.", file=sys.stderr)
        return 130
    if not args.verbose and progress.tail:
        print("Last output:\n" + "\n".join("    " + line for line in progress.tail), file=sys.stderr)
    print("Setup could not finish. Fix the cause above and run it again; finished steps are kept.",
          file=sys.stderr)
    return 1


# --------------------------------------------------------------------------- #
# run
# --------------------------------------------------------------------------- #
def cmd_run(args: argparse.Namespace) -> int:
    from pipeline.CLI import workflow_cfg
    from pipeline.Runner import logger, run_workflow

    if not args.output:
        print("run: --out is required", file=sys.stderr)
        return 2
    if not args.genomes and not args.resume:
        print("run: --genomes is required (or --resume to continue the run in --out)", file=sys.stderr)
        return 2
    try:
        cfg = workflow_cfg(args, logger)
        run_workflow(cfg)
    except SystemExit as exc:
        return int(exc.code or 0) if isinstance(exc.code, int) else 1
    except KeyboardInterrupt:
        logger.error("Stopped. `bactolomew run --out <same folder> --resume` continues the run.")
        return 130
    except Exception as exc:  # noqa: BLE001 - one line for the user, the log has the rest
        logger.error(f"Run failed: {exc}")
        return 1
    return 0


# --------------------------------------------------------------------------- #
# check
# --------------------------------------------------------------------------- #
_QT_PROBE = (
    "import PyQt6.QtCore as c; print('PyQt6 ' + c.PYQT_VERSION_STR + ', Qt ' + c.QT_VERSION_STR)\n"
    "import PyQt6.QtWebEngineWidgets\n"
    "import PyQt6Ads\n"
)


def _qt_report() -> tuple[bool, str]:
    """Qt imports in a child process, so a missing system library cannot stop the check."""
    try:
        res = subprocess.run([sys.executable, "-c", _QT_PROBE], capture_output=True, text=True, timeout=60,
                             env={**os.environ, "QT_QPA_PLATFORM": "offscreen"})
    except (OSError, subprocess.SubprocessError) as exc:
        return False, str(exc)
    if res.returncode == 0:
        return True, res.stdout.strip() + ", WebEngine and PyQt6Ads import"
    err = [line for line in res.stderr.strip().splitlines() if line.strip()]
    return False, (err[-1] if err else f"exit {res.returncode}")


# polars build flags -> /proc/cpuinfo names, where they differ
_CPUINFO_FLAG = {"sse3": "pni", "sse4.1": "sse4_1", "sse4.2": "sse4_2", "cmpxchg16b": "cx16", "lzcnt": "abm"}


def polars_cpu_problem(cpuinfo: str = "/proc/cpuinfo") -> str:
    """The CPU features the installed polars build needs and this CPU lacks ("" when none
    or unknown). Without them polars dies with an illegal instruction; polars[rtcompat]
    runs anywhere."""
    import importlib
    import importlib.util

    for name in ("_polars_runtime_compat", "_polars_runtime_64", "_polars_runtime_32"):
        if importlib.util.find_spec(name) is not None:
            break
    else:
        return ""
    flags = str(getattr(importlib.import_module(name), "BUILD_FEATURE_FLAGS", "") or "")
    try:
        with open(cpuinfo, encoding="utf-8") as fh:
            line = next((ln for ln in fh if ln.startswith("flags")), "")
    except OSError:
        return ""
    have = set(line.partition(":")[2].split())
    if not have:
        return ""
    need = [f[1:] for f in flags.split(",") if f.startswith("+") and f != "+ctr-static"]
    missing = [f for f in need if _CPUINFO_FLAG.get(f, f) not in have]
    return ", ".join(missing)


def _database_rows(row: Callable[[str, Optional[bool], str], None]) -> None:
    """The marker and compound databases, installed or not. They are optional (a run
    searches the installed ones), so a missing one is not a problem."""
    from pipeline.DbCatalog import COMPOUND_DB_CATALOG_PATH, MARKER_DB_CATALOG_PATH, catalog_databases

    for label, catalog in (("Marker DBs", MARKER_DB_CATALOG_PATH), ("Compound DBs", COMPOUND_DB_CATALOG_PATH)):
        try:
            databases = catalog_databases(catalog)
        except Exception as exc:  # noqa: BLE001 - the DB folder row above says what to do
            row(label, None, f"catalog not readable: {exc}")
            continue
        installed = sum(db.installed for db in databases)
        row(label, None, f"{installed} of {len(databases)} installed (optional; runs search the installed ones)")
        for db in databases:
            category = f" ({db.category})" if db.category else ""
            print(f"  {db.name + category:<52} {'installed' if db.installed else 'not installed'}")


def cmd_db(args: argparse.Namespace) -> int:
    from pipeline.UserDatabases import run_cli

    return run_cli(args)


def cmd_check(_args: argparse.Namespace) -> int:
    from pipeline import ContainerRuntime as CR
    from pipeline.Setup import SetupJob, _project_db_layout_ready, _project_db_root, _read_config
    from pipeline.ToolEnv import APP_MODE_VIEWER, _resolve_runtime, configured_app_mode, project_root, \
        resource_budget_text

    problems = 0

    def row(label: str, ok: Optional[bool], text: str) -> None:
        nonlocal problems
        mark = {True: "ok", False: "MISSING", None: "info"}[ok]
        if ok is False:
            problems += 1
        print(f"{label:<18} {mark:<8} {text}")

    root = project_root()
    print(_version_text())
    row("Python", sys.version_info >= (3, 10), f"{platform.python_version()} ({sys.executable})")
    libc, libc_version = platform.libc_ver()
    row("C library", None, f"{libc} {libc_version}".strip() or "unknown")
    qt_ok, qt_text = _qt_report()
    print(f"{'Qt (window)':<18} {'ok' if qt_ok else 'no':<8} {qt_text}")
    problem = display_problem()
    print(f"{'Display':<18} {'ok' if not problem else 'no':<8} "
          f"{'available' if not problem else 'none (the window cannot open; setup, run and check work)'}")
    config_path = root / "bactolomew_config.json"
    config = _read_config(root)
    row("Config", config_path.is_file(), str(config_path) + ("" if config_path.is_file() else " (run `setup`)"))
    viewer = config_path.is_file() and configured_app_mode(config) == APP_MODE_VIEWER
    row("Mode", None, "viewer (opens shared runs; `setup` adds the pipeline)" if viewer
        else "full (runs the pipeline)")
    row("Resources", None, resource_budget_text())
    missing_cpu = polars_cpu_problem()
    if missing_cpu:
        row("CPU", False, f"lacks {missing_cpu}, which the installed polars build needs; "
                          f"install the compatible one: {sys.executable} -m pip install 'polars[rtcompat]'")

    usable = [name for name, ok in CR.detect_runtimes().items() if ok is True]
    if viewer:
        # The tools image and databases are not needed to view runs: reported, never counted.
        row("Usable runtimes", None, ", ".join(usable) or "none")
        for label in ("Tool image", "Bakta database", "CheckM2 database", "DB folder"):
            row(label, None, "not needed (viewer)")
        print()
        print("Ready to view runs." if not problems else f"{problems} item(s) missing.")
        return 0 if not problems else 1
    if not config.get("runtime") and not os.environ.get("BACTOLOMEW_NATIVE"):
        # No runtime saved yet: report the one setup will pick, as setup does.
        if usable:
            config = {**config, "runtime": usable[0]}
            row("Runtime", None, f"{usable[0]} (not saved yet; setup will use it)")
        else:
            row("Runtime", False, "none configured or usable: install Docker, Podman or Apptainer (Singularity)")
            config = {**config, "runtime": "docker"}
    else:
        runtime = _resolve_runtime(config)
        if runtime == "native":
            row("Runtime", None, "native (tools from conda environments under envs/)")
        else:
            status = CR.runtime_status(runtime)
            row("Runtime", status.ok, runtime + ("" if status.ok else f": {status.message}"))
    row("Usable runtimes", None, ", ".join(usable) or "none")

    for step, found in SetupJob(root, config=config, log=lambda _line: None).status():
        if step.key != "check":
            row(step.label, found is not None, str(found) if found is not None else "not installed (run `setup`)")
    db_ready = _project_db_layout_ready(root)
    row("DB folder", db_ready, str(_project_db_root(root)) + ("" if db_ready else " (run `setup`)"))
    _database_rows(row)
    print()
    print("Ready to run." if not problems else f"{problems} item(s) missing: `bactolomew setup --yes` installs them.")
    return 0 if not problems else 1


# --------------------------------------------------------------------------- #
# Arguments
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    from pipeline.CLI import add_run_arguments
    from pipeline.UserDatabases import add_cli_arguments as add_db_arguments

    parser = argparse.ArgumentParser(
        prog="bactolomew",
        description="BACTolomew: bacterial comparative genomics. Without a command it opens the desktop app.",
    )
    parser.add_argument("--version", action="store_true", help="print the version and exit")
    sub = parser.add_subparsers(dest="command", metavar="{gui,setup,run,check,db}")

    gui = sub.add_parser("gui", help="open the desktop app (the default)")
    gui.add_argument("--dev", action="store_true",
                     help="Skip Bakta env/DB setup and use current hardcoded dev Bakta env+DB.")

    setup = sub.add_parser("setup", help="first-time setup without a window: tool image, databases, config")
    setup.add_argument("--yes", "-y", action="store_true", help="go ahead without asking (downloads over 40 GB)")
    setup.add_argument("--runtime", choices=RUNTIME_CHOICES,
                       help="container runtime to use and save (default: the configured one, else the first usable)")
    setup.add_argument("--sif", metavar="PATH",
                       help="Apptainer/Singularity image file of the tools (implies --runtime singularity)")
    setup.add_argument("--bakta", choices=("full",), default="full",
                       help="Bakta database type: full, the only one (a 39.4 GB download that needs about "
                            "121 GB free while it unpacks; a light database installed earlier is kept and used)")
    setup.add_argument("--viewer", action="store_true",
                       help="viewer-only install: open runs shared by others, download nothing "
                            "(plain `setup` later adds the pipeline)")
    setup.add_argument("--dry-run", action="store_true", help="show what would be installed, change nothing")
    setup.add_argument("--verbose", "-v", action="store_true", help="print every output line of the tools")

    run = sub.add_parser("run", help="run the pipeline on a genome folder (databases from the config)")
    add_run_arguments(run)

    sub.add_parser("check", help="report Python, Qt, the container runtime, the image and the databases")

    db = sub.add_parser("db", help="list, add or remove marker and compound databases (searched after Bakta)")
    add_db_arguments(db)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()
    if not argv or argv[0].startswith("-") and argv[0] not in ("-h", "--help", "--version"):
        argv = ["gui", *argv]  # no command: the desktop app (bactolomew --dev)
    args = parser.parse_args(argv)
    if args.version:
        print(_version_text())
        return 0
    handlers = {"gui": cmd_gui, "setup": cmd_setup, "run": cmd_run, "check": cmd_check, "db": cmd_db}
    return handlers[args.command or "gui"](args)


if __name__ == "__main__":
    sys.exit(main())
