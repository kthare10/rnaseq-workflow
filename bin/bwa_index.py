#!/usr/bin/env python3

"""Wrapper for BWA index construction.

Builds a BWA index from a reference genome FASTA file.
Corresponds to the MAKE_BWA_INDEX process in the Nextflow pipeline.
"""

import argparse
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description="Build BWA index")
    parser.add_argument("--reference", required=True, help="Reference genome FASTA file")
    args = parser.parse_args()

    # Build BWA index with prefix "ref_idx"
    cmd = ["bwa", "index", "-p", "ref_idx", args.reference]

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    print("BWA index created: ref_idx.{amb,ann,bwt,pac,sa}")


if __name__ == "__main__":
    main()
