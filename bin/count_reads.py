#!/usr/bin/env python3

"""Wrapper for gene-level read quantification using Rsubread featureCounts.

Calls count_reads.R which:
- Reads BAM files and counts mapped reads from .counts files
- Parses GFF annotation to extract gene info
- Runs Rsubread::featureCounts for gene-level quantification
- Generates library composition plots

Corresponds to the COUNT_READS process in the Nextflow pipeline.
"""

import argparse
import os
import shutil
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description="Gene-level read quantification")
    parser.add_argument("--metadata", required=True, help="Sample metadata TSV file")
    parser.add_argument("--gff", required=True, help="GFF annotation file")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("--bam", action="append", default=[], help="BAM file (repeatable)")
    parser.add_argument("--bai", action="append", default=[], help="BAI file (repeatable)")
    args = parser.parse_args()

    # Symlink hierarchical BAM/BAI paths to flat names in CWD.
    # The R script reads BAMs as paste0(sample, ".bam") (flat names).
    for bam_path in args.bam:
        flat_name = os.path.basename(bam_path)
        if bam_path != flat_name and not os.path.exists(flat_name):
            os.symlink(bam_path, flat_name)
    for bai_path in args.bai:
        flat_name = os.path.basename(bai_path)
        if bai_path != flat_name and not os.path.exists(flat_name):
            os.symlink(bai_path, flat_name)

    # The R script is staged into the working directory by Pegasus
    r_script = os.path.join(os.getcwd(), "count_reads.R")

    cmd = [
        "Rscript", r_script,
        "-m", args.metadata,
        "-g", args.gff,
        "-t", str(args.threads),
    ]

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    # Copy R outputs to read_counts/ subdirectory for organized staging
    os.makedirs("read_counts", exist_ok=True)
    expected = [
        "gene_counts.tsv",
        "gene_counts_pc.tsv",
        "counts_summary.tsv",
        "ref_gene_df.tsv",
        "library_composition.png",
        "library_composition_proportions.png",
    ]
    for f in expected:
        if os.path.exists(f):
            shutil.copy2(f, os.path.join("read_counts", f))
            print(f"Output: {f} -> read_counts/{f}")
        else:
            print(f"Warning: expected output not found: {f}", file=sys.stderr)


if __name__ == "__main__":
    main()
