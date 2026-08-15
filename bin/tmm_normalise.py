#!/usr/bin/env python3

"""Wrapper for TMM normalization of read counts using edgeR.

Calls TMM_normalise_counts.R which:
- Reads gene_counts.tsv and ref_gene_df.tsv from working directory
- Filters out rRNA genes
- Applies TMM normalization (edgeR)
- Computes CPM and RPKM values
- Optionally log-transforms

Corresponds to the TMM_NORMALISE_COUNTS process in the Nextflow pipeline.
"""

import argparse
import os
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description="TMM normalization of read counts")
    parser.add_argument(
        "--log-transform", type=str, default="TRUE",
        help="Log-transform counts (TRUE/FALSE, default: TRUE)",
    )
    args = parser.parse_args()

    # Create output subdirectory for organized staging
    os.makedirs("read_counts", exist_ok=True)

    # Symlink hierarchical inputs to flat names expected by R script
    for src in ["read_counts/gene_counts.tsv", "read_counts/ref_gene_df.tsv"]:
        flat = os.path.basename(src)
        if os.path.exists(src) and not os.path.exists(flat):
            os.symlink(src, flat)

    # The R script is staged into the working directory by Pegasus
    r_script = os.path.join(os.getcwd(), "TMM_normalise_counts.R")

    cmd = [
        "Rscript", r_script,
        "-t", args.log_transform,
        "-o", "read_counts/",
    ]

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    for f in ["read_counts/cpm_counts.tsv", "read_counts/rpkm_counts.tsv"]:
        if os.path.exists(f):
            print(f"Output: {f}")
        else:
            print(f"Warning: expected output not found: {f}", file=sys.stderr)


if __name__ == "__main__":
    main()
