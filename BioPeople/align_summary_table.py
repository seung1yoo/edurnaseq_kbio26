#!/usr/bin/env python3

import argparse
import csv
import glob
import re
import sys
from pathlib import Path


DEFAULT_PATTERN = "BioPeople/02.Align/*/align_summary.txt"

FIELDS = [
    "Sample",
    "Input Reads",
    "Mapped Reads",
    "Mapped Rate (%)",
    "Multi-mapped Reads",
    "Multi-mapped Rate (%)",
    "Concordant Pair Rate (%)",
]


def parse_int(value):
    return int(value.replace(",", ""))


def parse_section(text, section_name):
    pattern = re.compile(
        rf"{section_name} reads:\s*"
        r"Input\s*:\s*(?P<input>\d+)\s*"
        r"Mapped\s*:\s*(?P<mapped>\d+)\s*\((?P<mapped_pct>[\d.]+)% of input\)\s*"
        r"of these:\s*(?P<multiple>\d+)\s*\(\s*(?P<multiple_pct>[\d.]+)%\) "
        r"have multiple alignments \((?P<gt20>\d+) have >20\)",
        re.MULTILINE,
    )
    match = pattern.search(text)
    if not match:
        raise ValueError(f"Could not parse {section_name} reads section")

    return {
        "input": parse_int(match.group("input")),
        "mapped": parse_int(match.group("mapped")),
        "mapped_pct": float(match.group("mapped_pct")),
        "multiple": parse_int(match.group("multiple")),
        "multiple_pct": float(match.group("multiple_pct")),
        "multiple_gt20": parse_int(match.group("gt20")),
    }


def parse_align_summary(path):
    text = path.read_text()
    sample = path.parent.name

    left = parse_section(text, "Left")
    right = parse_section(text, "Right")

    concordant_match = re.search(r"(?P<rate>[\d.]+)% concordant pair alignment rate\.", text)

    if concordant_match is None:
        raise ValueError("Could not parse concordant pair alignment rate")

    input_reads = left["input"] + right["input"]
    mapped_reads = left["mapped"] + right["mapped"]
    multi_mapped_reads = left["multiple"] + right["multiple"]
    mapped_rate = mapped_reads / input_reads * 100
    multi_mapped_rate = multi_mapped_reads / mapped_reads * 100

    return {
        "Sample": sample,
        "Input Reads": f"{input_reads:,}",
        "Mapped Reads": f"{mapped_reads:,}",
        "Mapped Rate (%)": f"{mapped_rate:.1f}%",
        "Multi-mapped Reads": f"{multi_mapped_reads:,}",
        "Multi-mapped Rate (%)": f"{multi_mapped_rate:.1f}%",
        "Concordant Pair Rate (%)": f"{float(concordant_match.group('rate')):.1f}%",
    }


def expand_inputs(inputs):
    paths = []
    for value in inputs:
        matches = sorted(glob.glob(value))
        if matches:
            paths.extend(Path(match) for match in matches)
        else:
            paths.append(Path(value))
    return paths


def main():
    parser = argparse.ArgumentParser(
        description="Summarize TopHat2 align_summary.txt files into one TSV row per sample."
    )
    parser.add_argument(
        "inputs",
        nargs="*",
        default=[DEFAULT_PATTERN],
        help=f"Input align_summary.txt files or glob patterns. Default: {DEFAULT_PATTERN}",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output TSV path. Default: stdout",
    )
    args = parser.parse_args()

    paths = expand_inputs(args.inputs)
    if not paths:
        print("[ERROR] No input files found", file=sys.stderr)
        return 1

    rows = []
    for path in paths:
        if not path.is_file():
            print(f"[ERROR] Input file not found: {path}", file=sys.stderr)
            return 1
        try:
            rows.append(parse_align_summary(path))
        except ValueError as error:
            print(f"[ERROR] {path}: {error}", file=sys.stderr)
            return 1

    output_handle = open(args.output, "w", newline="") if args.output else sys.stdout
    try:
        writer = csv.DictWriter(output_handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    finally:
        if args.output:
            output_handle.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
