#!/usr/bin/env python3

"""Wrapper for BAM to BigWig conversion using deeptools bamCoverage.

Corresponds to the BAM2BIGWIG process in the Nextflow pipeline.
"""

import argparse
import os
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description="BAM to BigWig conversion")
    parser.add_argument("--input-bam", required=True, help="Input BAM file")
    parser.add_argument("--output", required=True, help="Output BigWig file")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("--bin-size", type=int, default=5, help="Bin size in bp (default: 5)")
    args = parser.parse_args()

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    cmd = [
        "bamCoverage",
        "-b", args.input_bam,
        "-o", args.output,
        "--binSize", str(args.bin_size),
        "-p", str(args.threads),
    ]

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
