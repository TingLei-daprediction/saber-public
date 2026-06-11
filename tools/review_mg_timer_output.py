#!/usr/bin/env python
"""Review and plot mg_timer_output timer files.

Supported formats:
1. Matrix/table format used by mg_timer_output, e.g.
      mype, init, upsend, ..., multiply, ..., icount
      0,    0.0000, ...
      1,    0.0000, ...
2. OOPS_STATS lines, e.g.
      OOPS_STATS label : total calls avg
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt


NUM_RE = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")
OOPS_STATS_RE = re.compile(
    r"^\s*OOPS_STATS\s+(?P<label>.+?)\s*:\s*"
    r"(?P<total>[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+"
    r"(?P<calls>\d+)\s+"
    r"(?P<avg>[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*$"
)


@dataclass
class TimerRow:
    label: str
    total_time: float
    calls: int | None
    avg_time: float | None
    raw_line: str
    mean_time: float | None = None
    max_time: float | None = None
    min_time: float | None = None
    num_ranks: int | None = None


def _clean_label(label: str) -> str:
    label = label.strip(" :-|,\t")
    label = re.sub(r"\s+", " ", label)
    return label


def _safe_float(value: str) -> float | None:
    try:
        x = float(value)
    except ValueError:
        return None
    if not math.isfinite(x):
        return None
    return x


def parse_oops_stats(lines: list[str]) -> list[TimerRow]:
    rows: list[TimerRow] = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        match = OOPS_STATS_RE.match(stripped)
        if not match:
            continue
        rows.append(
            TimerRow(
                label=_clean_label(match.group("label")),
                total_time=float(match.group("total")),
                calls=int(match.group("calls")),
                avg_time=float(match.group("avg")),
                raw_line=stripped,
            )
        )
    rows.sort(key=lambda row: row.total_time, reverse=True)
    return rows


def parse_matrix(lines: list[str]) -> list[TimerRow]:
    header_idx = None
    header = None
    for idx, line in enumerate(lines):
        if "," not in line:
            continue
        parts = [_clean_label(part) for part in line.split(",")]
        lowered = [part.lower() for part in parts]
        if "mype" in lowered and "icount" in lowered and len(parts) > 5:
            header_idx = idx
            header = parts
            break
    if header_idx is None or header is None:
        return []

    data_rows: list[list[float]] = []
    for line in lines[header_idx + 1 :]:
        if not line.strip():
            continue
        if "," not in line:
            continue
        parts = [part.strip() for part in line.split(",")]
        if len(parts) != len(header):
            continue
        values: list[float] = []
        ok = True
        for part in parts:
            value = _safe_float(part)
            if value is None:
                ok = False
                break
            values.append(value)
        if ok:
            data_rows.append(values)

    if not data_rows:
        return []

    name_to_idx = {name.lower(): idx for idx, name in enumerate(header)}
    mype_idx = name_to_idx.get("mype")
    icount_idx = name_to_idx.get("icount")
    num_ranks = len(data_rows)

    rows: list[TimerRow] = []
    for idx, name in enumerate(header):
        lname = name.lower()
        if idx == mype_idx or idx == icount_idx:
            continue

        values = [row[idx] for row in data_rows]
        total_time = sum(values)
        mean_time = total_time / num_ranks
        max_time = max(values)
        min_time = min(values)

        calls = None
        avg_time = None
        if icount_idx is not None:
            icount_values = [row[icount_idx] for row in data_rows]
            if all(abs(v - round(v)) < 1e-9 for v in icount_values):
                unique_counts = sorted({int(round(v)) for v in icount_values})
                if len(unique_counts) == 1 and unique_counts[0] > 0:
                    calls = unique_counts[0]
                    avg_time = mean_time / calls

        rows.append(
            TimerRow(
                label=name,
                total_time=total_time,
                calls=calls,
                avg_time=avg_time,
                raw_line="matrix aggregate",
                mean_time=mean_time,
                max_time=max_time,
                min_time=min_time,
                num_ranks=num_ranks,
            )
        )

    rows.sort(key=lambda row: row.total_time, reverse=True)
    return rows


def load_rows(path: Path) -> tuple[str, list[TimerRow]]:
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()

    matrix_rows = parse_matrix(lines)
    if matrix_rows:
        return "matrix", matrix_rows

    oops_rows = parse_oops_stats(lines)
    if oops_rows:
        return "oops_stats", oops_rows

    return "unknown", []


def write_csv(rows: list[TimerRow], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "label",
                "total_time",
                "calls",
                "avg_time",
                "mean_time",
                "max_time",
                "min_time",
                "num_ranks",
                "raw_line",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.label,
                    row.total_time,
                    row.calls,
                    row.avg_time,
                    row.mean_time,
                    row.max_time,
                    row.min_time,
                    row.num_ranks,
                    row.raw_line,
                ]
            )


def write_text_summary(rows: list[TimerRow], path: Path, top_n: int, format_name: str) -> None:
    total = sum(row.total_time for row in rows)
    with path.open("w", encoding="utf-8") as f:
        f.write(f"Detected format: {format_name}\n")
        f.write(f"Parsed timers: {len(rows)}\n")
        f.write(f"Summed total time: {total:.6f}\n\n")
        f.write(f"Top {min(top_n, len(rows))} timers by total time:\n")
        for idx, row in enumerate(rows[:top_n], start=1):
            frac = (row.total_time / total * 100.0) if total > 0 else 0.0
            calls = row.calls if row.calls is not None else "-"
            avg = f"{row.avg_time:.6f}" if row.avg_time is not None else "-"
            mean_time = f"{row.mean_time:.6f}" if row.mean_time is not None else "-"
            max_time = f"{row.max_time:.6f}" if row.max_time is not None else "-"
            min_time = f"{row.min_time:.6f}" if row.min_time is not None else "-"
            ranks = row.num_ranks if row.num_ranks is not None else "-"
            f.write(
                f"{idx:>2}. {row.label}\n"
                f"    total_time={row.total_time:.6f}  calls={calls}  avg_time={avg}  frac={frac:.2f}%\n"
                f"    mean_time={mean_time}  max_time={max_time}  min_time={min_time}  num_ranks={ranks}\n"
            )


def plot_top_timers(rows: list[TimerRow], path: Path, top_n: int) -> None:
    top = rows[:top_n]
    labels = [row.label for row in top][::-1]
    values = [row.total_time for row in top][::-1]

    fig, ax = plt.subplots(figsize=(12, max(6, 0.35 * len(top))))
    ax.barh(labels, values, color="#4472c4")
    ax.set_xlabel("Summed time across parsed rows")
    ax.set_ylabel("Timer")
    ax.set_title(f"Top {len(top)} Timers by Summed Time")
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_cumulative(rows: list[TimerRow], path: Path, top_n: int) -> None:
    top = rows[:top_n]
    total = sum(row.total_time for row in rows)
    cumulative = []
    running = 0.0
    for row in top:
        running += row.total_time
        cumulative.append(running / total * 100.0 if total > 0 else 0.0)

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(range(1, len(top) + 1), cumulative, marker="o", color="#c0504d")
    ax.set_xlabel("Top-N timers")
    ax.set_ylabel("Cumulative share of summed time (%)")
    ax.set_title("Cumulative Time Share")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_calls_vs_time(rows: list[TimerRow], path: Path, top_n: int) -> None:
    top = [row for row in rows if row.calls is not None][:top_n]
    if not top:
        return

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter([row.calls for row in top], [row.total_time for row in top], color="#70ad47")
    for row in top:
        ax.annotate(row.label, (row.calls, row.total_time), fontsize=8, alpha=0.8)
    ax.set_xlabel("Calls")
    ax.set_ylabel("Summed time")
    ax.set_title("Calls vs Summed Time")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Review mg_timer_output-style timer text files.")
    parser.add_argument("timer_file", type=Path, help="Path to mg_timer_output text file")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory for CSV/text/plots. Defaults to <timer_file stem>_review beside the input file.",
    )
    parser.add_argument("--top-n", type=int, default=25, help="Number of top timers to plot/summarize")
    args = parser.parse_args()

    timer_file = args.timer_file
    output_dir = args.output_dir or timer_file.with_name(f"{timer_file.stem}_review")
    output_dir.mkdir(parents=True, exist_ok=True)

    format_name, rows = load_rows(timer_file)
    if not rows:
        raise SystemExit(f"No timer rows could be parsed from {timer_file}")

    write_csv(rows, output_dir / "timer_summary.csv")
    write_text_summary(rows, output_dir / "timer_summary.txt", args.top_n, format_name)
    plot_top_timers(rows, output_dir / "top_timers.png", args.top_n)
    plot_cumulative(rows, output_dir / "cumulative_time_share.png", args.top_n)
    plot_calls_vs_time(rows, output_dir / "calls_vs_total_time.png", args.top_n)

    print(f"Detected format: {format_name}")
    print(f"Parsed timers: {len(rows)}")
    print(f"Output directory: {output_dir}")
    print("Wrote:")
    print(f"  {output_dir / 'timer_summary.csv'}")
    print(f"  {output_dir / 'timer_summary.txt'}")
    print(f"  {output_dir / 'top_timers.png'}")
    print(f"  {output_dir / 'cumulative_time_share.png'}")
    calls_plot = output_dir / "calls_vs_total_time.png"
    if calls_plot.exists():
        print(f"  {calls_plot}")


if __name__ == "__main__":
    main()
