#!/usr/bin/env python3
"""Run the finite-element solver benchmark matrix and save tidy CSV data."""

import argparse
import csv
import datetime
import platform
import re
import subprocess
from pathlib import Path


KEY_VALUE = re.compile(r"([a-z0-9_]+)=([^ ]+)")


def parse_line(line):
    return dict(KEY_VALUE.findall(line.strip()))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--levels", nargs="+", type=int, default=range(2, 9))
    parser.add_argument("--output", type=Path, default=Path(__file__).parent / "generated")
    args = parser.parse_args()

    project = Path(__file__).resolve().parents[1]
    executable = args.build.resolve() / "bin" / "fe_benchmark_mg"
    args.output.mkdir(parents=True, exist_ok=True)

    configurations = [
        ("mg-v", 40),
        ("mg-w", 40),
        ("gs", 2000),
        ("pcg-gs", 3000),
        ("pcg-mg-v", 100),
    ]
    results = []
    histories = []
    raw_log = []

    cases = [
        (level, method, max_iterations)
        for level in args.levels
        for method, max_iterations in configurations
    ]
    if 1 not in args.levels:
        cases.append((1, "gs", 2000))

    for level, method, max_iterations in cases:
        command = [
            str(executable), "--method", method, "--fine", str(level),
            "--max-iterations", str(max_iterations), "--tolerance", "1e-10",
            "--smoothing", "3", "--coarse-smoothing", "500", "--history",
        ]
        if method == "mg-v" and level <= 4:
            field_path = (args.output / f"fields_level{level}.csv").resolve()
            command.extend(["--field-output", str(field_path)])
        print(f"running level={level} method={method}", flush=True)
        completed = subprocess.run(
            command, cwd=args.build.resolve(), text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
        )
        log_command = [part.replace(str(project), "<project>") for part in command]
        raw_log.append(f"$ {' '.join(log_command)}\n{completed.stdout}")
        result_lines = [line for line in completed.stdout.splitlines() if line.startswith("RESULT ")]
        if len(result_lines) != 1:
            raise RuntimeError(f"No unique result for {method}, level {level}:\n{completed.stdout}")
        result = parse_line(result_lines[0])
        result["exit_code"] = str(completed.returncode)
        results.append(result)
        for line in completed.stdout.splitlines():
            if line.startswith("HISTORY "):
                histories.append(parse_line(line))

    def write_csv(path, rows):
        fields = list(rows[0])
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    write_csv(args.output / "results.csv", results)
    write_csv(args.output / "history.csv", histories)
    (args.output / "raw.log").write_text("\n\n".join(raw_log))
    (args.output / "metadata.txt").write_text(
        f"generated_utc={datetime.datetime.now(datetime.timezone.utc).isoformat()}\n"
        f"platform={platform.platform()}\n"
        f"python={platform.python_version()}\n"
        f"git_commit={subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=project, text=True).strip()}\n"
    )


if __name__ == "__main__":
    main()
