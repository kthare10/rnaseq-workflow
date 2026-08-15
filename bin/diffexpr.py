#!/usr/bin/env python3

"""Wrapper for differential gene expression analysis using DESeq2.

Calls diffexpr.R which:
- Reads gene_counts.tsv, sample_metadata.tsv, contrast_table.tsv from working dir
- Builds DESeqDataSet and runs DESeq2
- Performs pairwise contrasts for each row in the contrast table
- Generates EnhancedVolcano plots

Corresponds to the DIFF_EXPRESSION process in the Nextflow pipeline.
"""

import argparse
import os
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description="Differential gene expression (DESeq2)")
    parser.add_argument(
        "--p-threshold", type=float, default=0.05,
        help="Adjusted p-value threshold (default: 0.05)",
    )
    parser.add_argument(
        "--l2fc-threshold", type=float, default=1.0,
        help="Log2 fold change threshold (default: 1.0)",
    )
    args = parser.parse_args()

    # Create output subdirectory for organized staging
    os.makedirs("diff_expr", exist_ok=True)

    # Symlink hierarchical input to flat name expected by R script
    src = "read_counts/gene_counts.tsv"
    if os.path.exists(src) and not os.path.exists("gene_counts.tsv"):
        os.symlink(src, "gene_counts.tsv")

    # The R script is staged into the working directory by Pegasus
    r_script = os.path.join(os.getcwd(), "diffexpr.R")

    cmd = [
        "Rscript", r_script,
        "-p", str(args.p_threshold),
        "-l", str(args.l2fc_threshold),
        "-o", "diff_expr/",
    ]

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    # Report outputs found
    diff_dir = "diff_expr"
    if os.path.isdir(diff_dir):
        for f in os.listdir(diff_dir):
            if f.startswith("DGE_") or f.startswith("volcano_plot_"):
                print(f"Output: {diff_dir}/{f}")


if __name__ == "__main__":
    main()
