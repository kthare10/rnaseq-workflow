#!/usr/bin/env python3

"""
Pegasus workflow generator for bacterial RNA-Seq analysis.

Converted from the chienlab-rnaseq Nextflow DSL2 pipeline.

Pipeline steps:
1. fastp        - QC and adapter trimming (per-sample, parallel)
2. bwa_index    - Build BWA index from reference genome
3. bwa_align    - Align reads with BWA-MEM (per-sample, parallel)
4. bam2bigwig   - Generate BigWig files for visualization (per-sample, parallel)
5. count_reads  - Gene-level read quantification with Rsubread featureCounts (fan-in)
6. tmm_normalise - TMM normalization of read counts (edgeR)
7. pca          - Principal component analysis of samples
8. diffexpr     - Differential gene expression with DESeq2 (optional)

Usage:
    ./workflow_generator.py --sample-file samples.tsv --ref-genome ref.fasta \\
                            --ref-ann ref.gff --output workflow.yml

    # With differential expression:
    ./workflow_generator.py --sample-file samples.tsv --ref-genome ref.fasta \\
                            --ref-ann ref.gff --contrast-table contrasts.tsv \\
                            --output workflow.yml

    # A Slurm cluster:
    ./workflow_generator.py ... -e compute --site-style slurm \\
                            --queue cpu --project my_lab

Sites: the workflow names no scheduler. Jobs state cores, memory and a
wall-clock runtime; custom_sites.py writes sites.yml for an HTCondor pool by
default, or for Slurm via --site-style slurm --queue --project.
"""

import argparse
import csv
import logging
import os
import sys
from pathlib import Path

from Pegasus.api import *

# Site-catalog handling. Keeps this generator free of scheduler details.
sys.path.insert(0, str(Path(__file__).parent.resolve()))
from custom_sites import (  # noqa: E402
    HOSTED_SITE, STYLES, ensure_sites_yml, hosted_catalog, is_batch_site, parse_profile,
    parse_tag_profile, worker_package_url,
)

# Execution site when -e is not given: an HTCondor pool, unless ~/.pegasusrc
# names a hosted catalog (pegasushub pegasus-site-catalogs), whose one site is
# HOSTED_SITE ("compute").
DEFAULT_SITE = "condorpool"

# The container's OS, for the Pegasus worker package staged into it. Must
# run on the image's base (Apptainer/RNASeq_Container.def:
# mambaorg/micromamba:1.5-jammy, Ubuntu 22.04, glibc 2.35), not the submit
# host. Pegasus publishes no ubuntu_22 package, so use the oldest-glibc one,
# rhel_8 (glibc 2.28).
CONTAINER_PLATFORM = "x86_64_rhel_8"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Per-tool resource configuration — the only resource statements the workflow
# makes, portable on any site. Cores follow the Nextflow conf/base.config
# labels; memory is scaled down from its 24 GB (process_medium) and 46 GB
# (process_high) so every job fits a 16 GB node, which is ample for bacterial
# genomes. runtime is a wall-clock budget in seconds; batch sites (Slurm)
# refuse or kill jobs without one, so keep it generous.
TOOL_CONFIGS = {
    "fastp":          {"memory": "8 GB",  "cores": 4, "runtime": 3600},      # process_medium
    "bwa_index":      {"memory": "8 GB",  "cores": 4, "runtime": 1800},      # process_medium
    "bwa_align":      {"memory": "14 GB", "cores": 8, "runtime": 4 * 3600},  # process_high
    "bam2bigwig":     {"memory": "8 GB",  "cores": 4, "runtime": 3600},      # process_medium
    "count_reads":    {"memory": "8 GB",  "cores": 4, "runtime": 3600},      # process_medium
    "tmm_normalise":  {"memory": "8 GB",  "cores": 4, "runtime": 1800},      # process_medium
    "pca":            {"memory": "8 GB",  "cores": 4, "runtime": 1800},      # process_medium
    "diffexpr":       {"memory": "14 GB", "cores": 8, "runtime": 3600},      # process_high
}


class RnaseqWorkflow:
    """Pegasus workflow for bacterial RNA-Seq analysis."""

    wf = None
    tc = None
    rc = None
    props = None

    dagfile = None
    wf_dir = None
    shared_scratch_dir = None
    local_storage_dir = None
    wf_name = "rnaseq_workflow"
    # Set by main() from the site in use (see the site-catalog block there).
    worker_package_url = None
    bind_workflow_dir = False

    def __init__(self, args, dagfile="workflow.yml"):
        self.dagfile = dagfile
        self.wf_dir = str(Path(__file__).parent.resolve())
        self.shared_scratch_dir = os.path.join(self.wf_dir, "scratch")
        self.local_storage_dir = os.path.join(self.wf_dir, "output")
        self.args = args
        self.samples = self._parse_sample_file(args.sample_file, args.data_dir)
        self.contrasts = (
            self._parse_contrast_table(args.contrast_table)
            if args.contrast_table
            else []
        )

    def _parse_sample_file(self, sample_file, data_dir):
        """Parse sample TSV file and return list of sample dicts."""
        samples = []
        with open(sample_file, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if not row.get("sample"):
                    continue
                sample = {
                    "sample": row["sample"].strip(),
                    "file1": row["file1"].strip(),
                    "file2": row.get("file2", "").strip(),
                    "group": row.get("group", "").strip(),
                    "rep_no": row.get("rep_no", "").strip(),
                    "paired": row.get("paired", "0").strip(),
                    "strandedness": row.get("strandedness", "unstranded").strip(),
                }
                # Resolve file paths with data_dir
                if data_dir:
                    sample["file1_path"] = os.path.join(data_dir, sample["file1"])
                    if sample["file2"] and sample["paired"] == "1":
                        sample["file2_path"] = os.path.join(
                            data_dir, sample["file2"]
                        )
                    else:
                        sample["file2_path"] = ""
                else:
                    sample["file1_path"] = sample["file1"]
                    sample["file2_path"] = sample.get("file2", "")

                samples.append(sample)
        return samples

    def _parse_contrast_table(self, contrast_file):
        """Parse contrast table for differential expression."""
        contrasts = []
        with open(contrast_file, newline="") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader)  # skip header
            for row in reader:
                if len(row) >= 2:
                    contrasts.append(
                        {"group1": row[0].strip(), "group2": row[1].strip()}
                    )
        return contrasts

    def _generate_sample_metadata(self):
        """Generate sample_metadata.tsv for downstream R scripts."""
        metadata_path = os.path.join(self.wf_dir, "sample_metadata.tsv")
        with open(metadata_path, "w") as f:
            f.write(
                "sample\tfile1\tfile2\tgroup\trep_no\tpaired\tstrandedness\n"
            )
            for s in self.samples:
                file1 = s["file1"]
                file2 = s["file2"] if s["paired"] == "1" else ""
                f.write(
                    f"{s['sample']}\t{file1}\t{file2}\t{s['group']}\t"
                    f"{s['rep_no']}\t{s['paired']}\t{s['strandedness']}\n"
                )
        return metadata_path

    def write(self):
        """Write all catalogs and workflow to files."""
        self.props.write()
        self.rc.write()
        self.tc.write()
        self.wf.write(file=self.dagfile)

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------
    def create_pegasus_properties(self, sites_yml="sites.yml",
                                  bypass_input_staging=False):
        self.props = Properties()
        self.props["pegasus.transfer.threads"] = "16"
        # Symlink rather than copy when an input already sits on the
        # execution site. A no-op otherwise, so always on.
        self.props["pegasus.transfer.links"] = "true"
        if bypass_input_staging:
            # Jobs read inputs (notably the .sif) straight from the submit
            # host's paths. Only valid where workers share a filesystem with
            # the submit host (a Slurm cluster) — never on a condor pool
            # staging over HTCondor file transfer.
            self.props["pegasus.transfer.bypass.input.staging"] = "true"
        if self.worker_package_url:
            # Stage the container-compatible worker package named in the TC
            # (pegasus::worker) and never download one from inside a job:
            # the image may lack curl/wget, workers may lack internet, and the
            # submit host's kickstart may need a newer glibc than the image.
            self.props["pegasus.transfer.worker.package"] = "true"
            self.props["pegasus.transfer.worker.package.strict"] = "false"
            self.props["pegasus.transfer.worker.package.autodownload"] = "false"
        if os.path.isfile(sites_yml):
            # Lets pegasus-plan find sites.yml from any directory.
            self.props["pegasus.catalog.site"] = "YAML"
            self.props["pegasus.catalog.site.file"] = os.path.abspath(sites_yml)

    # ------------------------------------------------------------------
    # Transformation Catalog
    # ------------------------------------------------------------------
    def create_transformation_catalog(
        self,
        container_sif="Apptainer/RNASeq_Container.sif",
    ):
        """Container and transformations. Nothing here names a site: the
        scripts and the .sif live on "local" (the submit host) and are
        shipped to wherever the job runs."""
        self.tc = TransformationCatalog()

        # A local Apptainer .sif built with `apptainer build`. Pegasus stages
        # the file like any other input, so image_site is the site where the
        # .sif physically lives (the submit host = "local").
        sif_path = (
            container_sif
            if os.path.isabs(container_sif)
            else os.path.join(self.wf_dir, container_sif)
        )
        if not os.path.exists(sif_path):
            logger.warning(
                "Apptainer image not found at %s — build it first with: "
                "apptainer build %s Apptainer/RNASeq_Container.def",
                sif_path, sif_path)
        container = Container(
            "rnaseq_container",
            container_type=Container.SINGULARITY,
            image="file://" + sif_path,
            image_site="local",
        )
        if self.bind_workflow_dir:
            # Batch sites stage inputs as symlinks into the workflow
            # directory, and PegasusLite starts containers with --no-home, so
            # without this bind every job fails with kickstart "Unable to
            # execute the specified binary" (exit 127). Never on a condor
            # pool: the directory does not exist on its workers.
            container.add_pegasus_profile(
                container_arguments=f"--bind {self.wf_dir}")

        transformations = []
        for tool_name, config in TOOL_CONFIGS.items():
            tx = Transformation(
                tool_name,
                site="local",
                pfn=os.path.join(self.wf_dir, f"bin/{tool_name}.py"),
                is_stageable=True,
                container=container,
            ).add_pegasus_profile(
                memory=config["memory"],
                cores=config.get("cores", 1),
                runtime=config["runtime"],
            )
            transformations.append(tx)

        if self.worker_package_url:
            transformations.append(
                Transformation(
                    "worker",
                    namespace="pegasus",
                    site="local",
                    pfn=self.worker_package_url,
                    is_stageable=True,
                    arch=Arch.X86_64,
                    os_type=OS.LINUX,
                )
            )

        self.tc.add_containers(container)
        self.tc.add_transformations(*transformations)

    # ------------------------------------------------------------------
    # Replica Catalog
    # ------------------------------------------------------------------
    def create_replica_catalog(self):
        self.rc = ReplicaCatalog()

        # Reference genome FASTA
        self.rc.add_replica(
            "local",
            os.path.basename(self.args.ref_genome),
            "file://" + os.path.abspath(self.args.ref_genome),
        )

        # Reference annotation GFF
        self.rc.add_replica(
            "local",
            os.path.basename(self.args.ref_ann),
            "file://" + os.path.abspath(self.args.ref_ann),
        )

        # Sample metadata (generated at workflow creation time)
        metadata_path = self._generate_sample_metadata()
        self.rc.add_replica(
            "local",
            "sample_metadata.tsv",
            "file://" + os.path.abspath(metadata_path),
        )

        # R scripts (support files called by Python wrappers)
        r_scripts = [
            "count_reads.R",
            "TMM_normalise_counts.R",
            "pca.R",
            "diffexpr.R",
        ]
        for script in r_scripts:
            script_path = os.path.join(self.wf_dir, "bin", script)
            if os.path.exists(script_path):
                self.rc.add_replica(
                    "local", script, "file://" + script_path
                )

        # Contrast table (if provided)
        if self.args.contrast_table:
            self.rc.add_replica(
                "local",
                "contrast_table.tsv",
                "file://" + os.path.abspath(self.args.contrast_table),
            )

        # FASTQ files
        for s in self.samples:
            file1_path = os.path.abspath(s["file1_path"])
            if os.path.exists(file1_path):
                self.rc.add_replica(
                    "local",
                    os.path.basename(s["file1"]),
                    "file://" + file1_path,
                )
            else:
                logger.warning(f"FASTQ file not found: {file1_path}")

            if s["paired"] == "1" and s["file2_path"]:
                file2_path = os.path.abspath(s["file2_path"])
                if os.path.exists(file2_path):
                    self.rc.add_replica(
                        "local",
                        os.path.basename(s["file2"]),
                        "file://" + file2_path,
                    )
                else:
                    logger.warning(f"FASTQ file not found: {file2_path}")

    # ------------------------------------------------------------------
    # Workflow DAG
    # ------------------------------------------------------------------
    def create_workflow(self):
        self.wf = Workflow(self.wf_name, infer_dependencies=True)

        # Reference file objects
        ref_genome_file = File(os.path.basename(self.args.ref_genome))
        ref_ann_file = File(os.path.basename(self.args.ref_ann))
        metadata_file = File("sample_metadata.tsv")

        # R script file objects
        count_reads_r = File("count_reads.R")
        tmm_normalise_r = File("TMM_normalise_counts.R")
        pca_r = File("pca.R")
        diffexpr_r = File("diffexpr.R")

        # ---- BWA Index (single job) ----
        bwa_idx_files = []
        for ext in ["amb", "ann", "bwt", "pac", "sa"]:
            bwa_idx_files.append(File(f"ref_idx.{ext}"))

        bwa_index_job = (
            Job("bwa_index", _id="bwa_index", node_label="bwa_index")
            .add_args("--reference", ref_genome_file.lfn)
            .add_inputs(ref_genome_file)
            .add_outputs(
                *bwa_idx_files, stage_out=False, register_replica=False
            )
        )
        self.wf.add_jobs(bwa_index_job)

        # Track files for fan-in steps
        all_bam_files = []
        all_bai_files = []
        all_count_files = []

        # ---- Per-sample pipeline ----
        for s in self.samples:
            sample_id = s["sample"]
            is_paired = s["paired"] == "1"

            # Input FASTQ files
            read1_file = File(os.path.basename(s["file1"]))
            read2_file = (
                File(os.path.basename(s["file2"])) if is_paired else None
            )

            # -- FASTP (QC + trimming) --
            if is_paired:
                trimmed1 = File(f"{sample_id}_1_trimmed.fq.gz")
                trimmed2 = File(f"{sample_id}_2_trimmed.fq.gz")
            else:
                trimmed1 = File(f"{sample_id}_trimmed.fq.gz")
                trimmed2 = None
            fastp_html = File(f"fastp_qc/{sample_id}_fastp.html")

            fastp_args = [
                "--read1", read1_file.lfn,
                "--out1", trimmed1.lfn,
                "--html", fastp_html.lfn,
                "--threads", str(TOOL_CONFIGS['fastp']['cores']),
            ]
            fastp_inputs = [read1_file]

            if is_paired:
                fastp_args += [
                    "--read2", read2_file.lfn,
                    "--out2", trimmed2.lfn,
                ]
                fastp_inputs.append(read2_file)

            fastp_job = (
                Job(
                    "fastp",
                    _id=f"fastp_{sample_id}",
                    node_label=f"fastp_{sample_id}",
                )
                .add_args(*fastp_args)
                .add_inputs(*fastp_inputs)
                .add_outputs(
                    fastp_html, stage_out=True, register_replica=False
                )
                .add_outputs(
                    trimmed1, stage_out=False, register_replica=False
                )
                .add_pegasus_profiles(label=sample_id)
            )
            if is_paired:
                fastp_job.add_outputs(
                    trimmed2, stage_out=False, register_replica=False
                )
            self.wf.add_jobs(fastp_job)

            # -- BWA ALIGN --
            bam_file = File(f"bwa_aln/{sample_id}.bam")
            bai_file = File(f"bwa_aln/{sample_id}.bam.bai")
            counts_file = File(f"{sample_id}.counts")
            all_bam_files.append(bam_file)
            all_bai_files.append(bai_file)
            all_count_files.append(counts_file)

            bwa_args = [
                "--read1", trimmed1.lfn,
                "--output-bam", bam_file.lfn,
                "--output-bai", bai_file.lfn,
                "--output-counts", counts_file.lfn,
                "--threads", str(TOOL_CONFIGS['bwa_align']['cores']),
            ]
            bwa_inputs = [trimmed1] + bwa_idx_files
            if is_paired:
                bwa_args += ["--read2", trimmed2.lfn]
                bwa_inputs.append(trimmed2)

            bwa_job = (
                Job(
                    "bwa_align",
                    _id=f"bwa_{sample_id}",
                    node_label=f"bwa_{sample_id}",
                )
                .add_args(*bwa_args)
                .add_inputs(*bwa_inputs)
                .add_outputs(
                    bam_file, stage_out=True, register_replica=False
                )
                .add_outputs(
                    bai_file, stage_out=True, register_replica=False
                )
                .add_outputs(
                    counts_file, stage_out=False, register_replica=False
                )
                .add_pegasus_profiles(label=sample_id)
            )
            self.wf.add_jobs(bwa_job)

            # -- BAM2BIGWIG (parallel branch from BAM) --
            bw_file = File(f"bigwig/{sample_id}.bw")
            bw_job = (
                Job(
                    "bam2bigwig",
                    _id=f"bw_{sample_id}",
                    node_label=f"bw_{sample_id}",
                )
                .add_args(
                    "--input-bam", bam_file.lfn,
                    "--output", bw_file.lfn,
                    "--threads", str(TOOL_CONFIGS['bam2bigwig']['cores']),
                )
                .add_inputs(bam_file, bai_file)
                .add_outputs(
                    bw_file, stage_out=True, register_replica=False
                )
                .add_pegasus_profiles(label=sample_id)
            )
            self.wf.add_jobs(bw_job)

        # ---- COUNT_READS (fan-in: all BAMs → gene counts) ----
        # Hierarchical output files (used for both organized output and
        # inter-job data flow; downstream wrappers symlink to flat names)
        gene_counts_out = File("read_counts/gene_counts.tsv")
        gene_counts_pc_out = File("read_counts/gene_counts_pc.tsv")
        ref_gene_df_out = File("read_counts/ref_gene_df.tsv")
        counts_summary_out = File("read_counts/counts_summary.tsv")
        lib_comp_out = File("read_counts/library_composition.png")
        lib_comp_prop_out = File(
            "read_counts/library_composition_proportions.png"
        )

        # Build --bam and --bai args for symlinking in wrapper
        bam_args = [tok for f in all_bam_files for tok in ("--bam", f.lfn)]
        bai_args = [tok for f in all_bai_files for tok in ("--bai", f.lfn)]

        count_reads_job = (
            Job(
                "count_reads",
                _id="count_reads",
                node_label="count_reads",
            )
            .add_args(
                "--metadata", "sample_metadata.tsv",
                "--gff", ref_ann_file.lfn,
                "--threads", str(TOOL_CONFIGS['count_reads']['cores']),
                *bam_args, *bai_args,
            )
            .add_inputs(
                metadata_file,
                ref_ann_file,
                count_reads_r,
                *all_bam_files,
                *all_bai_files,
                *all_count_files,
            )
            # Hierarchical outputs (also used as inputs by downstream jobs)
            .add_outputs(
                gene_counts_out, gene_counts_pc_out,
                ref_gene_df_out, counts_summary_out,
                lib_comp_out, lib_comp_prop_out,
                stage_out=True, register_replica=False,
            )
        )
        self.wf.add_jobs(count_reads_job)

        # ---- TMM_NORMALISE ----
        cpm_counts_out = File("read_counts/cpm_counts.tsv")
        rpkm_counts_out = File("read_counts/rpkm_counts.tsv")

        tmm_job = (
            Job(
                "tmm_normalise",
                _id="tmm_normalise",
                node_label="tmm_normalise",
            )
            .add_args("--log-transform", "TRUE")
            .add_inputs(gene_counts_out, ref_gene_df_out, tmm_normalise_r)
            .add_outputs(
                cpm_counts_out, rpkm_counts_out,
                stage_out=True, register_replica=False,
            )
        )
        self.wf.add_jobs(tmm_job)

        # ---- PCA ----
        pca_rds = File("PCA_samples/pca.rds")
        pca_coords = File("PCA_samples/pca_coords.tsv")
        pca_plot = File("PCA_samples/pca_grouped.png")

        pca_job = (
            Job("pca", _id="pca_samples", node_label="pca_samples")
            .add_inputs(cpm_counts_out, metadata_file, pca_r)
            .add_outputs(
                pca_rds,
                pca_coords,
                pca_plot,
                stage_out=True,
                register_replica=False,
            )
        )
        self.wf.add_jobs(pca_job)

        # ---- DIFF_EXPRESSION (optional) ----
        if self.args.contrast_table and self.contrasts:
            contrast_file = File("contrast_table.tsv")

            diffexpr_outputs = []
            for c in self.contrasts:
                name = f"{c['group1']}_{c['group2']}"
                diffexpr_outputs.append(
                    File(f"diff_expr/DGE_{name}.tsv")
                )
                diffexpr_outputs.append(
                    File(f"diff_expr/volcano_plot_{name}.png")
                )

            diffexpr_job = (
                Job("diffexpr", _id="diffexpr", node_label="diffexpr")
                .add_args(
                    "--p-threshold", str(self.args.p_thresh),
                    "--l2fc-threshold", str(self.args.l2fc_thresh),
                )
                .add_inputs(
                    gene_counts_out,
                    metadata_file,
                    contrast_file,
                    diffexpr_r,
                )
            )
            for out_f in diffexpr_outputs:
                diffexpr_job.add_outputs(
                    out_f, stage_out=True, register_replica=False
                )
            self.wf.add_jobs(diffexpr_job)

        logger.info(
            f"Workflow created: {len(self.samples)} sample(s), "
            f"{len(self.contrasts)} contrast(s)"
        )


# ======================================================================
# main()
# ======================================================================
def main():
    parser = argparse.ArgumentParser(
        description="Pegasus RNA-Seq Workflow Generator "
        "(converted from chienlab-rnaseq Nextflow pipeline)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --sample-file samples.tsv --ref-genome ref.fasta --ref-ann ref.gff
  %(prog)s --sample-file samples.tsv --ref-genome ref.fasta --ref-ann ref.gff \\
           --contrast-table contrasts.tsv --data-dir /path/to/fastq/
  %(prog)s ... -e compute --site-style slurm --queue cpu --project my_lab
""",
    )

    # Standard Pegasus arguments. Every site option has a default, so a run
    # with only the inputs plans on an HTCondor pool, or on a hosted
    # catalog's "compute" when ~/.pegasusrc names one; the rest tailor
    # sites.yml for a batch cluster. dest "execution_site" is what Studio
    # matches to keep its "Where it runs" choice in sync.
    parser.add_argument(
        "-e",
        "--execution-site",
        "--execution-site-name",
        dest="execution_site",
        metavar="STR",
        type=str,
        default=None,
        help=f"Site to plan against (default: {HOSTED_SITE!r} when "
             "~/.pegasusrc names a hosted catalog, which call their site "
             f"that; otherwise {DEFAULT_SITE!r})",
    )
    parser.add_argument(
        "--site-style",
        choices=("auto",) + STYLES + ("none",),
        default="auto",
        help="auto (default): keep a sites.yml entry or hosted catalog if one "
             "exists, else add an HTCondor site. condor/slurm: (re)write this "
             "site's entry. none: leave sites.yml alone.",
    )
    parser.add_argument("--queue", metavar="PARTITION",
                        help="Batch partition/queue jobs submit to")
    parser.add_argument("--project", metavar="ACCOUNT",
                        help="Allocation/account charged on a batch site")
    parser.add_argument("--site-scratch", metavar="DIR",
                        help="Slurm: shared scratch visible to workers and "
                             "the submit host (default: ./work)")
    parser.add_argument("--site-profile", action="append", default=[],
                        type=parse_profile, metavar="NS:KEY=VALUE",
                        help="Extra profile on the execution site; repeatable")
    parser.add_argument("--tag-profile", action="append", default=[],
                        type=parse_tag_profile, metavar="TAG:NS:KEY=VALUE",
                        help="Profile for jobs carrying a tag; repeatable")
    parser.add_argument("--shared-filesystem", choices=("auto", "yes", "no"),
                        default="auto",
                        help="Let jobs read inputs directly from the submit "
                             "host. auto: on for Slurm sites, off for HTCondor.")
    parser.add_argument("--sites-yml", metavar="FILE", default="sites.yml",
                        help="Local site catalog (default: sites.yml)")
    parser.add_argument("-s", "--skip-sites-catalog", action="store_true",
                        help="Deprecated: same as --site-style none")
    parser.add_argument(
        "-o",
        "--output",
        metavar="STR",
        type=str,
        default="workflow.yml",
        help="Output file (default: workflow.yml)",
    )

    # Workflow-specific arguments
    parser.add_argument(
        "--sample-file",
        type=str,
        required=True,
        help="Sample metadata TSV file "
        "(columns: sample, file1, file2, group, rep_no, paired, strandedness)",
    )
    parser.add_argument(
        "--ref-genome",
        type=str,
        required=True,
        help="Reference genome FASTA file",
    )
    parser.add_argument(
        "--ref-ann",
        type=str,
        required=True,
        help="Reference annotation GFF file",
    )
    parser.add_argument(
        "--data-dir",
        type=str,
        default=None,
        help="Directory containing FASTQ files "
        "(prepended to file1/file2 in sample sheet)",
    )
    parser.add_argument(
        "--contrast-table",
        type=str,
        default=None,
        help="Contrast table TSV for differential expression (optional)",
    )
    parser.add_argument(
        "--p-thresh",
        type=float,
        default=0.05,
        help="Adjusted p-value threshold for DE (default: 0.05)",
    )
    parser.add_argument(
        "--l2fc-thresh",
        type=float,
        default=1.0,
        help="Log2 fold change threshold for DE (default: 1.0)",
    )
    parser.add_argument(
        "--container-sif",
        type=str,
        default="Apptainer/RNASeq_Container.sif",
        help="Path to the Apptainer .sif image, absolute or relative to the "
             "workflow directory (default: Apptainer/RNASeq_Container.sif)",
    )

    args = parser.parse_args()
    if args.execution_site is None:
        args.execution_site = HOSTED_SITE if hosted_catalog() else DEFAULT_SITE

    # Validation
    if not os.path.exists(args.sample_file):
        logger.error(f"Sample file not found: {args.sample_file}")
        sys.exit(1)
    if not os.path.exists(args.ref_genome):
        logger.error(f"Reference genome not found: {args.ref_genome}")
        sys.exit(1)
    if not os.path.exists(args.ref_ann):
        logger.error(f"Reference annotation not found: {args.ref_ann}")
        sys.exit(1)

    logger.info("=" * 70)
    logger.info("RNA-SEQ WORKFLOW GENERATOR")
    logger.info("=" * 70)
    logger.info(f"Sample file:    {args.sample_file}")
    logger.info(f"Reference:      {args.ref_genome}")
    logger.info(f"Annotation:     {args.ref_ann}")
    logger.info(f"Data directory: {args.data_dir or '(none)'}")
    logger.info(f"Contrast table: {args.contrast_table or '(none)'}")
    logger.info(f"Execution site: {args.execution_site}")
    logger.info(f"Output file:    {args.output}")
    logger.info("=" * 70)

    try:
        workflow = RnaseqWorkflow(args, dagfile=args.output)

        logger.info(f"Parsed {len(workflow.samples)} sample(s)")
        for s in workflow.samples:
            pe = "PE" if s["paired"] == "1" else "SE"
            logger.info(
                f"  {s['sample']} ({pe}, {s['strandedness']}, "
                f"group={s['group']})"
            )
        if workflow.contrasts:
            logger.info(f"Contrasts: {len(workflow.contrasts)}")
            for c in workflow.contrasts:
                logger.info(f"  {c['group1']} vs {c['group2']}")

        # --- Site catalog and the settings that depend on the site ---
        if args.skip_sites_catalog:
            args.site_style = "none"
        action, style = ensure_sites_yml(
            args.sites_yml, args.execution_site, workflow.wf_dir,
            style=args.site_style, queue=args.queue, project=args.project,
            scratch=args.site_scratch, profiles=args.site_profile,
            tag_profiles=args.tag_profile)
        hosted = hosted_catalog()
        logger.info(f"Site catalog: {args.sites_yml}: {action}"
                    + (f" (merged over hosted {hosted})" if hosted else ""))
        if (style is None and hosted and args.execution_site != "local"
                and args.execution_site != HOSTED_SITE):
            # Nothing was written for this site, so planning works only if
            # the hosted catalog happens to define it.
            logger.warning(
                f"{args.execution_site!r} is not defined in {args.sites_yml} "
                f"and hosted catalogs normally define only {HOSTED_SITE!r}: "
                f"pegasus-plan will fail unless {hosted} has it. Use "
                f"-e {HOSTED_SITE}, or --site-style condor/slurm to describe "
                f"{args.execution_site!r}.")
        if args.shared_filesystem == "auto":
            bypass = style is not None and style != "condor"
        else:
            bypass = args.shared_filesystem == "yes"
        workflow.bind_workflow_dir = bypass or is_batch_site(style)
        workflow.worker_package_url = worker_package_url(CONTAINER_PLATFORM)
        if not workflow.worker_package_url:
            logger.warning("pegasus-version not found: Pegasus will choose the "
                           "container's worker package itself (needs curl/wget "
                           "in the image and internet on the workers)")

        workflow.create_pegasus_properties(
            sites_yml=args.sites_yml, bypass_input_staging=bypass)
        workflow.create_transformation_catalog(
            container_sif=args.container_sif,
        )
        workflow.create_replica_catalog()
        workflow.create_workflow()
        workflow.write()

        logger.info(f"\nWorkflow written to {args.output}")
        logger.info(
            f"Submit: pegasus-plan --submit "
            f"-s {args.execution_site} -o local {args.output}"
        )

    except Exception as e:
        logger.error(f"Failed to generate workflow: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
