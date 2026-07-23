"""Run a command and record process-tree resource use for later review.

The monitor is deliberately side-band: it does not determine scientific DAG
dependencies.  Thresholds create structured warnings but never terminate the
wrapped command.  CPU concurrency must still be limited by the command itself.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable
from uuid import uuid4

import psutil


SCHEMA_VERSION = "0.1.0"
GIB = 1024**3


@dataclass(frozen=True)
class Snapshot:
    timestamp_utc: str
    elapsed_seconds: float
    process_count: int
    cpu_percent_machine_capacity: float
    rss_gib: float
    vms_gib: float
    io_read_gib: float
    io_write_gib: float
    project_size_gib: float | None
    filesystem_free_gib: float
    filesystem_used_percent: float


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _portable_path(path: str | Path, project_root: Path) -> str:
    resolved = Path(path).resolve()
    try:
        relative = resolved.relative_to(project_root)
    except ValueError:
        return f"<external>/{resolved.name}"
    return "." if not relative.parts else relative.as_posix()


def _portable_command(command: Iterable[str], project_root: Path) -> list[str]:
    """Redact absolute paths while retaining project-relative provenance."""

    portable: list[str] = []
    for raw_argument in command:
        argument = str(raw_argument)
        prefix, separator, value = argument.partition("=")
        candidate = value if separator else argument
        if Path(candidate).is_absolute():
            replacement = _portable_path(candidate, project_root)
            portable.append(
                f"{prefix}={replacement}" if separator else replacement
            )
        else:
            portable.append(argument)
    return portable


def _atomic_json(path: str | Path, payload: dict[str, Any]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.parent / f".{output.name}.{uuid4().hex}.tmp"
    try:
        temporary.write_text(
            json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        temporary.replace(output)
    finally:
        if temporary.exists():
            temporary.unlink()


def _project_size_bytes(path: Path) -> int:
    completed = subprocess.run(
        ["du", "-sb", "--", str(path)],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return int(completed.stdout.split(maxsplit=1)[0])


def _process_tree(process: psutil.Process) -> list[psutil.Process]:
    try:
        descendants = process.children(recursive=True)
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        descendants = []
    candidates = [process, *descendants]
    alive: list[psutil.Process] = []
    seen: set[tuple[int, float]] = set()
    for candidate in candidates:
        try:
            identity = (candidate.pid, candidate.create_time())
            if candidate.is_running() and identity not in seen:
                alive.append(candidate)
                seen.add(identity)
        except (psutil.NoSuchProcess, psutil.ZombieProcess, psutil.AccessDenied):
            continue
    return alive


def _tree_metrics(
    processes: Iterable[psutil.Process],
) -> tuple[float, int, int, int, int, int]:
    cpu_seconds = 0.0
    rss = 0
    vms = 0
    read_bytes = 0
    write_bytes = 0
    count = 0
    for process in processes:
        try:
            cpu = process.cpu_times()
            memory = process.memory_info()
            cpu_seconds += float(cpu.user + cpu.system)
            rss += int(memory.rss)
            vms += int(memory.vms)
            try:
                io = process.io_counters()
            except (psutil.AccessDenied, AttributeError):
                io = None
            if io is not None:
                read_bytes += int(io.read_bytes)
                write_bytes += int(io.write_bytes)
            count += 1
        except (psutil.NoSuchProcess, psutil.ZombieProcess, psutil.AccessDenied):
            continue
    return cpu_seconds, rss, vms, read_bytes, write_bytes, count


def _warning(
    warnings: dict[str, dict[str, Any]],
    *,
    code: str,
    value: float,
    threshold: float,
    message: str,
    timestamp: str,
) -> None:
    record = warnings.get(code)
    if record is None:
        warnings[code] = {
            "code": code,
            "first_timestamp_utc": timestamp,
            "threshold": threshold,
            "peak_observed": value,
            "message": message,
        }
    else:
        record["peak_observed"] = max(float(record["peak_observed"]), value)


def monitor_command(
    *,
    command: list[str],
    cwd: str | Path,
    project_root: str | Path,
    series_output: str | Path,
    summary_output: str | Path,
    command_log: str | Path,
    interval_seconds: float = 10.0,
    disk_interval_seconds: float = 60.0,
    cpu_warn_percent: float = 40.0,
    project_warn_gib: float = 10.0,
    project_critical_gib: float = 20.0,
    filesystem_free_warn_gib: float = 20.0,
    environment: dict[str, str] | None = None,
) -> dict[str, Any]:
    if not command:
        raise ValueError("command must not be empty")
    if interval_seconds <= 0 or disk_interval_seconds <= 0:
        raise ValueError("sampling intervals must be positive")
    thresholds = (
        cpu_warn_percent,
        project_warn_gib,
        project_critical_gib,
        filesystem_free_warn_gib,
    )
    if any(value < 0 for value in thresholds):
        raise ValueError("warning thresholds must be non-negative")
    if project_critical_gib < project_warn_gib:
        raise ValueError("project critical threshold must be >= warning threshold")

    working_directory = Path(cwd).resolve()
    project = Path(project_root).resolve()
    if not working_directory.is_dir():
        raise NotADirectoryError(working_directory)
    if not project.is_dir():
        raise NotADirectoryError(project)
    series_path = Path(series_output)
    log_path = Path(command_log)
    series_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    started_utc = _utc_now()
    started_monotonic = time.monotonic()
    logical_cpus = psutil.cpu_count(logical=True) or 1
    warnings: dict[str, dict[str, Any]] = {}
    snapshots: list[Snapshot] = []
    last_project_size: int | None = None
    next_disk_sample = 0.0
    previous_cpu_seconds: float | None = None
    previous_cpu_timepoint: float | None = None

    with log_path.open("wb") as command_handle, series_path.open(
        "w", encoding="utf-8", newline=""
    ) as series_handle:
        writer = csv.DictWriter(
            series_handle,
            fieldnames=list(Snapshot.__dataclass_fields__),
            delimiter="\t",
        )
        writer.writeheader()
        series_handle.flush()
        process = subprocess.Popen(
            command,
            cwd=working_directory,
            env=environment,
            stdout=command_handle,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        root_process = psutil.Process(process.pid)
        interrupted = False
        try:
            while True:
                now_monotonic = time.monotonic()
                elapsed = now_monotonic - started_monotonic
                tree = _process_tree(root_process)
                cpu_seconds, rss, vms, read_bytes, write_bytes, count = _tree_metrics(tree)
                if previous_cpu_seconds is None or previous_cpu_timepoint is None:
                    cpu_percent = 0.0
                else:
                    delta_cpu = max(0.0, cpu_seconds - previous_cpu_seconds)
                    delta_wall = max(1e-9, now_monotonic - previous_cpu_timepoint)
                    cpu_percent = 100.0 * delta_cpu / (delta_wall * logical_cpus)
                previous_cpu_seconds = cpu_seconds
                previous_cpu_timepoint = now_monotonic

                if elapsed >= next_disk_sample or last_project_size is None:
                    last_project_size = _project_size_bytes(project)
                    next_disk_sample = elapsed + disk_interval_seconds
                disk = shutil.disk_usage(project)
                project_gib = (
                    float(last_project_size / GIB)
                    if last_project_size is not None
                    else None
                )
                timestamp = _utc_now()
                snapshot = Snapshot(
                    timestamp_utc=timestamp,
                    elapsed_seconds=round(elapsed, 6),
                    process_count=count,
                    cpu_percent_machine_capacity=cpu_percent,
                    rss_gib=rss / GIB,
                    vms_gib=vms / GIB,
                    io_read_gib=read_bytes / GIB,
                    io_write_gib=write_bytes / GIB,
                    project_size_gib=project_gib,
                    filesystem_free_gib=disk.free / GIB,
                    filesystem_used_percent=100.0 * disk.used / disk.total,
                )
                snapshots.append(snapshot)
                writer.writerow(asdict(snapshot))
                series_handle.flush()

                if cpu_percent > cpu_warn_percent:
                    _warning(
                        warnings,
                        code="CPU_CAPACITY_ABOVE_WARNING",
                        value=cpu_percent,
                        threshold=cpu_warn_percent,
                        message=(
                            "Observed process-tree CPU exceeded the configured share "
                            "of whole-machine logical CPU capacity."
                        ),
                        timestamp=timestamp,
                    )
                if project_gib is not None and project_gib > project_warn_gib:
                    _warning(
                        warnings,
                        code="PROJECT_SIZE_ABOVE_WARNING",
                        value=project_gib,
                        threshold=project_warn_gib,
                        message="Project output size exceeded the warning threshold.",
                        timestamp=timestamp,
                    )
                if project_gib is not None and project_gib > project_critical_gib:
                    _warning(
                        warnings,
                        code="PROJECT_SIZE_ABOVE_CRITICAL",
                        value=project_gib,
                        threshold=project_critical_gib,
                        message="Project output size exceeded the critical threshold.",
                        timestamp=timestamp,
                    )
                free_gib = disk.free / GIB
                if free_gib < filesystem_free_warn_gib:
                    _warning(
                        warnings,
                        code="FILESYSTEM_FREE_BELOW_WARNING",
                        value=free_gib,
                        threshold=filesystem_free_warn_gib,
                        message="Filesystem free space fell below the warning threshold.",
                        timestamp=timestamp,
                    )

                exit_code = process.poll()
                if exit_code is not None:
                    break
                try:
                    process.wait(timeout=interval_seconds)
                except subprocess.TimeoutExpired:
                    pass
        except (KeyboardInterrupt, SystemExit):
            interrupted = True
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                process.wait(timeout=15)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                process.wait()
            raise
        finally:
            command_handle.flush()

    finished_utc = _utc_now()
    wall_seconds = time.monotonic() - started_monotonic
    exit_code = int(process.returncode)
    peak = lambda field: max((float(getattr(row, field)) for row in snapshots), default=0.0)
    final_project_size = _project_size_bytes(project) / GIB
    summary: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "status": "interrupted" if interrupted else ("success" if exit_code == 0 else "failed"),
        "exit_code": exit_code,
        "started_utc": started_utc,
        "finished_utc": finished_utc,
        "wall_seconds": wall_seconds,
        "command": _portable_command(command, project),
        "cwd": _portable_path(working_directory, project),
        "project_root": ".",
        "logical_cpu_count": logical_cpus,
        "sampling": {
            "interval_seconds": interval_seconds,
            "disk_interval_seconds": disk_interval_seconds,
            "n_snapshots": len(snapshots),
        },
        "thresholds": {
            "cpu_warn_percent_machine_capacity": cpu_warn_percent,
            "project_warn_gib": project_warn_gib,
            "project_critical_gib": project_critical_gib,
            "filesystem_free_warn_gib": filesystem_free_warn_gib,
            "threshold_action": "warn_only",
        },
        "peaks": {
            "cpu_percent_machine_capacity": peak("cpu_percent_machine_capacity"),
            "rss_gib": peak("rss_gib"),
            "vms_gib": peak("vms_gib"),
            "io_read_gib_observed": peak("io_read_gib"),
            "io_write_gib_observed": peak("io_write_gib"),
            "project_size_gib": max(
                final_project_size, peak("project_size_gib")
            ),
            "filesystem_used_percent": peak("filesystem_used_percent"),
        },
        "final_project_size_gib": final_project_size,
        "warnings": sorted(warnings.values(), key=lambda row: row["code"]),
        "outputs": {
            "series_tsv": _portable_path(series_path, project),
            "command_log": _portable_path(log_path, project),
        },
        "measurement_notes": [
            "CPU is normalized to total logical-CPU machine capacity.",
            "RSS and I/O are sampled across the live descendant process tree; very short-lived descendants can be undercounted.",
            "Thresholds do not terminate the command.",
            "Absolute paths are redacted or converted to project-relative paths.",
        ],
    }
    _atomic_json(summary_output, summary)
    return summary


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cwd", default=".")
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--series-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--command-log", required=True)
    parser.add_argument("--interval-seconds", type=float, default=10.0)
    parser.add_argument("--disk-interval-seconds", type=float, default=60.0)
    parser.add_argument("--cpu-warn-percent", type=float, default=40.0)
    parser.add_argument("--project-warn-gib", type=float, default=10.0)
    parser.add_argument("--project-critical-gib", type=float, default=20.0)
    parser.add_argument("--filesystem-free-warn-gib", type=float, default=20.0)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    return parser


def main() -> None:
    arguments = _parser().parse_args()
    command = list(arguments.command)
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        raise SystemExit("A command is required after --")
    summary = monitor_command(
        command=command,
        cwd=arguments.cwd,
        project_root=arguments.project_root,
        series_output=arguments.series_output,
        summary_output=arguments.summary_output,
        command_log=arguments.command_log,
        interval_seconds=arguments.interval_seconds,
        disk_interval_seconds=arguments.disk_interval_seconds,
        cpu_warn_percent=arguments.cpu_warn_percent,
        project_warn_gib=arguments.project_warn_gib,
        project_critical_gib=arguments.project_critical_gib,
        filesystem_free_warn_gib=arguments.filesystem_free_warn_gib,
        environment=os.environ.copy(),
    )
    raise SystemExit(int(summary["exit_code"]))


if __name__ == "__main__":
    main()
