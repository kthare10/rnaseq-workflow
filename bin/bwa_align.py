#!/usr/bin/env python3

"""Wrapper for BWA-MEM read alignment.

Aligns reads to a reference genome, sorts to BAM, creates index,
and counts uniquely mapped reads.
Handles both paired-end and single-end reads.
Corresponds to the BWA_ALIGN process in the Nextflow pipeline.
"""

import argparse
import os
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description="BWA-MEM alignment")
    parser.add_argument("--read1", required=True, help="Read 1 FASTQ file")
    parser.add_argument("--read2", default=None, help="Read 2 FASTQ file (paired-end)")
    parser.add_argument("--output-bam", required=True, help="Output sorted BAM file")
    parser.add_argument("--output-bai", required=True, help="Output BAM index file")
    parser.add_argument("--output-counts", required=True, help="Output unique mapped read count")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads")
    args = parser.parse_args()

    # Create output directories if needed
    for outfile in [args.output_bam, args.output_bai, args.output_counts]:
        out_dir = os.path.dirname(outfile)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)

    sort_threads = max(1, args.threads - 1)

    # BWA MEM + samtools sort pipeline
    if args.read2:
        bwa_cmd = (
            f"bwa mem -t {args.threads} ref_idx {args.read1} {args.read2} | "
            f"samtools sort -@ {sort_threads} -O bam - > {args.output_bam}"
        )
    else:
        bwa_cmd = (
            f"bwa mem -t {args.threads} ref_idx {args.read1} | "
            f"samtools sort -@ {sort_threads} -O bam - > {args.output_bam}"
        )

    print(f"Running: {bwa_cmd}")
    result = subprocess.run(bwa_cmd, shell=True, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    # Index the BAM file
    idx_cmd = f"samtools index -@ {args.threads} {args.output_bam}"
    print(f"Running: {idx_cmd}")
    result = subprocess.run(idx_cmd, shell=True, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    # Count uniquely mapped reads
    count_cmd = (
        f"samtools view -F 0x4 {args.output_bam} | "
        f"cut -f 1 | sort | uniq | wc -l > {args.output_counts}"
    )
    print(f"Running: {count_cmd}")
    result = subprocess.run(count_cmd, shell=True, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    print(f"Output BAM: {args.output_bam}")
    print(f"Output BAI: {args.output_bai}")
    print(f"Output counts: {args.output_counts}")


if __name__ == "__main__":
    main()
