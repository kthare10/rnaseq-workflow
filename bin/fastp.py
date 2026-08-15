#!/usr/bin/env python3

"""Wrapper for fastp QC and adapter trimming.

Handles both paired-end and single-end reads.
Corresponds to the FASTP process in the Nextflow pipeline.
"""

import argparse
import os
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description="fastp QC and trimming")
    parser.add_argument("--read1", required=True, help="Read 1 FASTQ file")
    parser.add_argument("--read2", default=None, help="Read 2 FASTQ file (paired-end)")
    parser.add_argument("--out1", required=True, help="Trimmed read 1 output")
    parser.add_argument("--out2", default=None, help="Trimmed read 2 output (paired-end)")
    parser.add_argument("--html", required=True, help="fastp HTML report output")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    args = parser.parse_args()

    # Create output directories if needed
    for outfile in [args.out1, args.out2, args.html]:
        if outfile:
            out_dir = os.path.dirname(outfile)
            if out_dir:
                os.makedirs(out_dir, exist_ok=True)

    # Build fastp command
    cmd = ["fastp", "-i", args.read1, "-o", args.out1]

    if args.read2 and args.out2:
        cmd.extend(["-I", args.read2, "-O", args.out2])

    cmd.extend([
        "--thread", str(args.threads),
        "--html", args.html,
    ])

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    print(f"Output: {args.out1}")
    if args.out2:
        print(f"Output: {args.out2}")
    print(f"Report: {args.html}")


if __name__ == "__main__":
    main()
