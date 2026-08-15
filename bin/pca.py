#!/usr/bin/env python3

"""Wrapper for PCA analysis of normalized expression values.

Calls pca.R which:
- Reads cpm_counts.tsv and sample_metadata.tsv from working directory
- Performs PCA on transposed count matrix
- Generates ggbiplot visualization colored by experimental group

Corresponds to the PCA_SAMPLES process in the Nextflow pipeline.
"""

import os
import subprocess
import sys


def main():
    # Create output subdirectory for organized staging
    os.makedirs("PCA_samples", exist_ok=True)

    # Symlink hierarchical input to flat name expected by R script
    src = "read_counts/cpm_counts.tsv"
    if os.path.exists(src) and not os.path.exists("cpm_counts.tsv"):
        os.symlink(src, "cpm_counts.tsv")

    # The R script is staged into the working directory by Pegasus
    r_script = os.path.join(os.getcwd(), "pca.R")

    cmd = ["Rscript", r_script, "-o", "PCA_samples/"]

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    for f in ["PCA_samples/pca.rds", "PCA_samples/pca_coords.tsv", "PCA_samples/pca_grouped.png"]:
        if os.path.exists(f):
            print(f"Output: {f}")
        else:
            print(f"Warning: expected output not found: {f}", file=sys.stderr)


if __name__ == "__main__":
    main()
