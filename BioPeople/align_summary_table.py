#!/usr/bin/env python3

import argparse
import csv
import glob
import re
import sys
from pathlib import Path


DEFAULT_PATTERN = "BioPeople/02.Align/*/align_summary.txt"

FIELDS = [
    "sample",
    "left_input",
    "left_mapped",
    "left_mapped_pct",
    "left_multiple",
    "left_multiple_pct",
    "left_multiple_gt20",
    "right_input",
    "right_mapped",
    "right_mapped_pct",
    "right_multiple",
    "right_multiple_pct",
    "right_multiple_gt20",
    "overall_mapping_rate",
    "aligned_pairs",
    "pair_multiple",
    "pair_multiple_pct",
    "discordant",
    "discordant_pct",
    "concordant_pair_alignment_rate",
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

    overall_match = re.search(r"(?P<rate>[\d.]+)% overall read mapping rate\.", text)
    aligned_pairs_match = re.search(r"Aligned pairs:\s*(?P<count>\d+)", text)
    pair_multiple_match = re.search(
        r"Aligned pairs:\s*\d+\s*"
        r"of these:\s*(?P<count>\d+)\s*\(\s*(?P<pct>[\d.]+)%\) have multiple alignments",
        text,
        re.MULTILINE,
    )
    discordant_match = re.search(
        r"(?P<count>\d+)\s*\(\s*(?P<pct>[\d.]+)%\) are discordant alignments",
        text,
    )
    concordant_match = re.search(r"(?P<rate>[\d.]+)% concordant pair alignment rate\.", text)

    required = {
        "overall mapping rate": overall_match,
        "aligned pairs": aligned_pairs_match,
        "pair multiple alignments": pair_multiple_match,
        "discordant alignments": discordant_match,
        "concordant pair alignment rate": concordant_match,
    }
    missing = [name for name, match in required.items() if match is None]
    if missing:
        raise ValueError(f"Could not parse {', '.join(missing)}")

    return {
        "sample": sample,
        "left_input": left["input"],
        "left_mapped": left["mapped"],
        "left_mapped_pct": left["mapped_pct"],
        "left_multiple": left["multiple"],
        "left_multiple_pct": left["multiple_pct"],
        "left_multiple_gt20": left["multiple_gt20"],
        "right_input": right["input"],
        "right_mapped": right["mapped"],
        "right_mapped_pct": right["mapped_pct"],
        "right_multiple": right["multiple"],
        "right_multiple_pct": right["multiple_pct"],
        "right_multiple_gt20": right["multiple_gt20"],
        "overall_mapping_rate": float(overall_match.group("rate")),
        "aligned_pairs": parse_int(aligned_pairs_match.group("count")),
        "pair_multiple": parse_int(pair_multiple_match.group("count")),
        "pair_multiple_pct": float(pair_multiple_match.group("pct")),
        "discordant": parse_int(discordant_match.group("count")),
        "discordant_pct": float(discordant_match.group("pct")),
        "concordant_pair_alignment_rate": float(concordant_match.group("rate")),
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
        description="Convert TopHat align_summary.txt files into a TSV table."
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
