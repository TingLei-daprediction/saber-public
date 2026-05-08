#!/usr/bin/env python
"""Review and plot mg_timer_output-style timer text files.

This script is intentionally tolerant because timer outputs vary across builds.
It scans each line, tries to extract a timer label plus one or more numeric
columns, and uses the last numeric value as the total time unless a better
pattern is obvious.
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


@dataclass
class TimerRow:
    label: str
    total_time: float
    calls: int | None
    avg_time: float | None
    raw_line: str


def _clean_label(label: str) -> str:
    label = label.strip(" :-|,\t")
    label = re.sub(r"\s+", " ", label)
    return label


def _looks_like_label(label: str) -> bool:
    if not label:
        return False
    if len(label) < 3:
        return False
    if not any(ch.isalpha() for ch in label):
        return False
    if label.lower().startswith(("rank ", "thread ", "time ", "total ")):
        return True
    return True


def parse_timer_line(line: str) -> TimerRow | None:
    stripped = line.strip()
    if not stripped:
        return None
    if stripped.startswith(("#", "=", "-", "*")):
        return None

    numbers = list(NUM_RE.finditer(stripped))
    if not numbers:
        return None

    label = _clean_label(stripped[: numbers[0].start()])
    if not _looks_like_label(label):
        return None

    values = [float(match.group(0)) for match in numbers]
    if not values:
        return None

    total_time = values[-1]
    if not math.isfinite(total_time):
        return None
    if total_time < 0:
        return None

    calls = None
    avg_time = None

    # Common loose heuristic:
    #   label ... <calls> <time>
    # or
    #   label ... <calls> <avg> <time>
    if len(values) >= 2:
        maybe_calls = values[-2]
        if float(maybe_calls).is_integer() and maybe_calls > 0:
            calls = int(maybe_calls)
            if len(values) >= 3:
                maybe_avg = values[-3]
                if maybe_avg >= 0 and total_time >= maybe_avg:
                    avg_time = maybe_avg

    if calls is not None and avg_time is None and calls > 0:
        avg_time = total_time / calls

    return TimerRow(
        label=label,
        total_time=total_time,
        calls=calls,
        avg_time=avg_time,
        raw_line=stripped,
    )


def load_rows(path: Path) -> list[TimerRow]:
    rows: list[TimerRow] = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            row = parse_timer_line(line)
            if row is not None:
                rows.append(row)
    return rows


def aggregate_rows(rows: list[TimerRow]) -> list[TimerRow]:
    by_label: dict[str, TimerRow] = {}
    for row in rows:
        if row.label not in by_label:
            by_label[row.label] = TimerRow(
                label=row.label,
                total_time=row.total_time,
                calls=row.calls,
                avg_time=row.avg_time,
                raw_line=row.raw_line,
            )
            continue

        existing = by_label[row.label]
        existing.total_time += row.total_time
        if existing.calls is not None and row.calls is not None:
            existing.calls += row.calls
        else:
            existing.calls = existing.calls or row.calls
        if existing.calls:
            existing.avg_time = existing.total_time / existing.calls
    out = list(by_label.values())
    out.sort(key=lambda row: row.total_time, reverse=True)
    return out


def write_csv(rows: list[TimerRow], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["label", "total_time", "calls", "avg_time", "raw_line"])
        for row in rows:
            writer.writerow([row.label, row.total_time, row.calls, row.avg_time, row.raw_line])


def write_text_summary(rows: list[TimerRow], path: Path, top_n: int) -> None:
    total = sum(row.total_time for row in rows)
    with path.open("w", encoding="utf-8") as f:
        f.write(f"Parsed timers: {len(rows)}\n")
        f.write(f"Summed total time: {total:.6f}\n\n")
        f.write(f"Top {min(top_n, len(rows))} timers by total time:\n")
        for idx, row in enumerate(rows[:top_n], start=1):
            frac = (row.total_time / total * 100.0) if total > 0 else 0.0
            calls = row.calls if row.calls is not None else "-"
            avg = f"{row.avg_time:.6f}" if row.avg_time is not None else "-"
            f.write(
                f"{idx:>2}. {row.label}\n"
                f"    total_time={row.total_time:.6f}  calls={calls}  avg_time={avg}  frac={frac:.2f}%\n"
            )


def plot_top_timers(rows: list[TimerRow], path: Path, top_n: int) -> None:
    top = rows[:top_n]
    labels = [row.label for row in top][::-1]
    values = [row.total_time for row in top][::-1]

    fig, ax = plt.subplots(figsize=(12, max(6, 0.35 * len(top))))
    ax.barh(labels, values, color="#4472c4")
    ax.set_xlabel("Total time")
    ax.set_ylabel("Timer")
    ax.set_title(f"Top {len(top)} Timers by Total Time")
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
    ax.set_ylabel("Cumulative share of total time (%)")
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
    ax.set_ylabel("Total time")
    ax.set_title("Calls vs Total Time")
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

    rows = load_rows(timer_file)
    rows = aggregate_rows(rows)
    if not rows:
        raise SystemExit(f"No timer rows could be parsed from {timer_file}")

    write_csv(rows, output_dir / "timer_summary.csv")
    write_text_summary(rows, output_dir / "timer_summary.txt", args.top_n)
    plot_top_timers(rows, output_dir / "top_timers.png", args.top_n)
    plot_cumulative(rows, output_dir / "cumulative_time_share.png", args.top_n)
    plot_calls_vs_time(rows, output_dir / "calls_vs_total_time.png", args.top_n)

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
